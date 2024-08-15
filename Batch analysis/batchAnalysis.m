function batchAnalysis(S, noiseFileIDs, maskSelections, Options)
% batchAnalysis processes SMLM data (for the OneFlowTraX app).
%
% Syntax:
%   batchAnalysis(S)
%
% Input Arguments:
%   (Required)
%   S                  Structure that holds the batch analysis parameters. See |OneFlowTraX.mlapp| for its required
%                      contents.
%                      (1,1) struct
%
%   noiseFileIDs       Assigns each input file in S.files.inputFiles a corresponding noise file, with |noiseFileIDs|
%                      having the same length as S.files.inputFiles, and each entry is an integer corresponding to the
%                      row ID of the noise file in S.files.noiseFiles.
%                      (:,1) double
%
%   maskSelections     Logical list with the same length as S.files.inputFiles, indicating which input files will have a
%                      mask applied to them.
%                      (:,1) logical
%
% Name-Value Pair Input Arguments:
%   (Optional)
%   updateProgress     When set to "true", this creates and updates a progress dialog box to inform the user about the
%                      current progress of the batch analysis. Alternatively, the handle to an existing progress dialog
%                      box can be provided (default: false).
%                      (1,1) logical | ProgressDialog object
%
% Other required m-files: [mList,pList] = matlab.codetools.requiredFilesAndProducts('fname.m')
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
%
% Tested: MATLAB Version: 9.13.0.2166757 (R2022b)
%	      Microsoft Windows 10 Enterprise Version 10.0 (Build 19045)
%
% Author: Sven zur Oven-Krockhaus
%	      Institute of Physical and Theoretical Chemistry
%	      University of Tuebingen, Tuebingen, Germany
% E-mail: sven.zur-oven-krockhaus@uni-tuebingen.de
%
% GNU placeholder
%
% Initial release: 2023-03-01
% Last revision: 2024-08-15

%% Function argument validation
arguments
    S (1,1) struct
    noiseFileIDs (:,1) double
    maskSelections (:,1) logical
    Options.updateProgress (1,1) = false
end

% Check if a GPU can be used (for the localization process and the variance map calculation).
useGPU = canUseGPU();

if islogical(Options.updateProgress)
    if Options.updateProgress == true
        hFig = uifigure;
        progressDialog = uiprogressdlg(hFig, 'Title', 'Batch Processing', 'Indeterminate','on');
    else
        % When no progress dialog is desired, we just make this object a structure. The assignments in the main code
        % (like progressDialog.Message = '...') then just stores data as a field without further consequences in the
        % GUI. This removes the need for conditional >if< statements.
        progressDialog = struct;
        progressDialog.CancelRequested = false;
    end
elseif isa(Options.updateProgress,'matlab.ui.dialog.ProgressDialog')
    progressDialog = Options.updateProgress;
end

%% Main

% Parse all parameters currently set in |S|, keep only the necessary ones and write them to the output file's parameter
% sheet.
O = formatParamsOutput(S);
O = O(S.general.analysisSteps{:,:});
O = vertcat(O{:});
writecell(O, S.files.outputFile,'Sheet','parameters');

% Calculate and store the variance maps for all noise files.
progressDialog.Message = 'Calculating variance maps for all selected noise files ...';
varianceMapArray = cell(size(S.files.noiseFiles));
for iNoiseFile = 1:numel(varianceMapArray)
    varianceMapArray{iNoiseFile} = calculateVarianceMap(S.files.noiseFiles{iNoiseFile},'useGPU',useGPU);
end

% Set to progress dialog to show the progressing bar and make it cancelable.
progressDialog.Indeterminate = 0;
progressDialog.Cancelable = 1;

% Loop through each file.
nFiles = numel(S.files.inputFiles);
for iFile = 1:nFiles

    % This is the point where batch processing is stopped (just before a new file would be analyzed).
    checkCancel(progressDialog);
    if progressDialog.CancelRequested == true
        break
    end

    [~, thisFilename] = fileparts(S.files.inputFiles{iFile});

    %% Localization
    progressDialog.Message = sprintf('%s: Localization ...', thisFilename);
    progressDialog.Value = (iFile-1)/nFiles;
    checkCancel(progressDialog);

    if iFile == 1
        % For the first file, create the sheets 'metadata' and 'localization results' and write the column headers.
        writecell({'file ID','filename','date','noise file','mask file','width (pixels)','height (pixels)','N frames'},...
            S.files.outputFile,'Sheet','metadata');
        writecell({'file ID','Original number of localizations','NaN / 0 / inf','obvious outliers','fitting errors',...
        'Remaining after pre-filtering','Remaining after histogram filtering'}, ...
        S.files.outputFile,'Sheet','localization results');
    end

    % Prepare the localization parameters and select the variance map according to the user-specified ID for this input
    % file.
    L = S.localization;
    histogramLimits = table2array(L.limits)';
    L = rmfield(L, 'limits');
    L.useGPU = useGPU;
    L = reshape(namedargs2cell(L),2,[]);
    L{2,1} = varianceMapArray{noiseFileIDs(iFile)};

    % Check the TIFF file integrity.
    try
        % Try opening the TIFF file using the MATLAB Gateway to LibTIFF library routines. Even for big image stacks,
        % this will only read the first frame (so it is still reasonably fast). If everything is fine, close the
        % Tiff object again and try to read it with |readTIFFstack|.
        warning('off','imageio:tiffutils:libtiffWarning'); 
        t = Tiff(S.files.inputFiles{iFile});
        warning('on','imageio:tiffutils:libtiffWarning'); 
        close(t);
        [imageStack, Metadata] = readTIFFstack(S.files.inputFiles{iFile});
    catch
        % If any error was encountered, skip this file entirely. In the output file, this will result in empty rows and
        % columns where data for this file would normally appear.
        continue
    end
   
    progressDialog.Value = (iFile-1+0.2)/nFiles;
    checkCancel(progressDialog);

    % Check the data again.
    if size(imageStack,3) == 1
        % When only one frame is present, this might be file that was put in the file list by mistake (e.g. a mask
        % file). Skip this file entirely.
        continue
    end

    % Populate the 'metadata' sheet, writing the date, ...
    writecell(horzcat({iFile},S.files.inputFiles{iFile},Metadata.FileModDate),S.files.outputFile,...
        'Sheet','metadata','Range',strcat('A',num2str(iFile+1)));

    % ... noise filename, ...
    writecell(S.files.noiseFiles(noiseFileIDs(iFile)),S.files.outputFile,...
        'Sheet','metadata','Range',strcat('D',num2str(iFile+1)));

    % ... mask filename (if applicable), ...
    if maskSelections(iFile) == true
        writecell(S.files.maskFiles(iFile),S.files.outputFile,'Sheet','metadata','Range',strcat('E',num2str(iFile+1)));
    end

    % ... and the image dimensions.
    writecell(num2cell([Metadata.Width,Metadata.Height,Metadata.Frames]),S.files.outputFile,...
        'Sheet','metadata','Range',strcat('F',num2str(iFile+1)));

    % For EMCCD files, overwrite the variance map with a generated dummy variance map that fits the dimensions of the
    % input data.
    if numel(L{2,1}) == 1
        dummyVarianceMap = zeros(Metadata.Height,Metadata.Width,'single');
        L{2,1} = dummyVarianceMap;
    end

    % Perform the localization.
    [pointList, outlierInfo] = localizeSMLMdata(imageStack, L{:});
   
    % Filter out localizations according to the histogram limits set by the user.    
    pointList = pointList(~any(pointList(:,[4 5 7]) < histogramLimits(1,:) |...
        pointList(:,[4 5 7]) > histogramLimits(2,:),2),:);

    % Update the outlier list and write it to the sheet 'localization results'.
    outlierInfo(6) = size(pointList,1);
    writecell(num2cell([iFile, outlierInfo]), S.files.outputFile,'Sheet','localization results',...
        'Range',strcat('A',num2str(iFile+1)));

    % Skip the rest of the analysis when the file has less than 10 localizations.
    if size(pointList,1) < 10
        continue
    end

    if S.general.rawDataOutput.("raw data output")(1)
        % The user specified that the localization data should be stored in a text file.
        [pathname, filename] = fileparts(S.files.inputFiles{iFile});
        outputFilename = fullfile(pathname, strcat(filename, S.general.rawDataOutput.("file suffix")(1), '.dat'));

        % Format the point list as a table to create the correct column headers and write the text file.
        T = array2table(pointList, ...
            'VariableNames', {'frame','posX_nm','posY_nm','sdvPSF_nm','photons','background','locprec_nm'});
        writetable(T, outputFilename, "Delimiter", "tab");
    end

    if S.general.analysisSteps.doAnalysis(2)
        %% Tracking
        % Tracking was chosen by the user to be included in the batch analysis, so discriminate between using the track
        % centroids or the localizations as base data for the potentially subsequent cluster analysis.
        progressDialog.Message = sprintf('%s: Tracking ...', thisFilename);
        progressDialog.Value = (iFile-1+0.4)/nFiles;
        checkCancel(progressDialog);

        if iFile == 1
            % If this is the first file, create the spreadsheet 'tracking results' and write the column headers.
            writecell({'file ID','Original number of tracks','Remaining tracks after filtering', ...
                'Remaining tracks after filtering and masking'}, S.files.outputFile,'Sheet','tracking results');
        end

        if contains(S.cluster.input, 'track')
            % Build the tracks based on all localizations.
            tmpTracks = assembleTracks(pointList, 'maxGapClose', S.tracking.maxGapClose, ...
                'maxLinkDistance', S.tracking.maxLinkDistance, 'exportAllInfo', true);
            
            % Get the limits that determine which tracks will be excluded.
            trackLimits = table2array(S.tracking.limits)';

            % Calculate the number of localizations, the duration, the displacement and range for each track.
            trackFeatures = [...,
                cellfun('size', tmpTracks, 1), ...
                (cellfun(@(x) range(x(:,1)),tmpTracks) + 1) / S.MSDanalysis.frameRate, ...
                cellfun(@(x) norm(diff(x([1,end],2:3))), tmpTracks), ...
                cellfun(@(x) hypot(range(x(:,2)), range(x(:,3))), tmpTracks)];

            % Remove invalid tracks.
            isValid = all(trackFeatures >= trackLimits(1,:) & trackFeatures <= trackLimits(2,:), 2);
            trackArray = tmpTracks(isValid);

            % Summarize the track results and write them to the sheet 'tracking results'.
            writecell(num2cell([iFile, numel(tmpTracks), sum(isValid), sum(isValid)]), S.files.outputFile, ...
                'Sheet', 'tracking results', 'Range', strcat('A',num2str(iFile + 1)));

            if maskSelections(iFile) == true
                % If a mask polygon is to be applied, remove all localizations accordingly. The polygon is based on the
                % original image size in pixels, so multiply with the pixel size to apply it to the localization
                % coordinates.
                maskPolygon = S.tracking.maskPolygons{iFile};
                maskPolygon.Vertices = maskPolygon.Vertices * S.localization.pixelSize;
                TFin = isinterior(maskPolygon, pointList(:,2:3));
                pointList = pointList(TFin,:);

                % Repeat the track building, overwriting the corresponding variables and also update the entry in the
                % 'tracking results' sheet.
                tmpTracks = assembleTracks(pointList, 'maxGapClose', S.tracking.maxGapClose, ...
                    'maxLinkDistance', S.tracking.maxLinkDistance, 'exportAllInfo', true);
                trackFeatures = [...,
                    cellfun('size', tmpTracks, 1), ...
                    (cellfun(@(x) range(x(:,1)),tmpTracks) + 1) / S.MSDanalysis.frameRate, ...
                    cellfun(@(x) norm(diff(x([1,end],2:3))), tmpTracks), ...
                    cellfun(@(x) hypot(range(x(:,2)), range(x(:,3))), tmpTracks)];
                isValid = all(trackFeatures >= trackLimits(1,:) & trackFeatures <= trackLimits(2,:), 2);
                trackArray = tmpTracks(isValid);
                writecell({sum(isValid)}, S.files.outputFile, ...
                'Sheet', 'tracking results', 'Range', strcat('D',num2str(iFile + 1)));
            end

            % Use the track centroids as points for the cluster analysis (NASTIC and segNASTIC will use |trackArray|
            % anyway).
            pointList = double(convertTracksToCentroids(trackArray));

        else
            % Use the individual localizations as points (mask files do not apply here).
            pointList = double(pointList(:,2:3));
        end

        if S.general.rawDataOutput.("raw data output")(2)
            % The user specified that the tracks should be stored in a text file.
            [pathname, filename] = fileparts(S.files.inputFiles{iFile});
            outputFilename = fullfile(pathname, strcat(filename, S.general.rawDataOutput.("file suffix")(2), '.dat'));

            % Format the tracks as a table to create the correct column headers and write the text file.
            trackID = repelem(1:numel(trackArray),cellfun('size', trackArray, 1))';
            trackList = vertcat(trackArray{:});
            T = array2table([trackID,trackList(:,1:7)], 'VariableNames',...
                {'trackID','frame','posX_nm','posY_nm','PSFsize_nm','photons','background','locPrecision_nm'});
            writetable(T, outputFilename, "Delimiter", "tab");
        end
    else
        % The user chose to not include tracking in the batch analysis, so just store the individual localizations for
        % the potentially subsequent cluster analysis.
        pointList = double(pointList(:,2:3));
    end

    if S.general.analysisSteps.doAnalysis(3)
        %% MSD analysis
        progressDialog.Message = sprintf('%s: Mobility analysis ...', thisFilename);
        progressDialog.Value = (iFile-1+0.6)/nFiles;
        checkCancel(progressDialog);

         if iFile == 1
            % If this is the first file, create the corresponding spreadsheets and write their column headers.
            writecell({'file ID';'lag time (s)'}, S.files.outputFile,'Sheet','mean MSD');
            writecell({'file ID'}, S.files.outputFile,'Sheet','D lists');
            writecell({'file ID','N','D (µm²/s)'}, S.files.outputFile,'Sheet','D fit (1 pop)');
            writecell({'file ID','N','slower D (µm²/s)','faster D (µm²/s)','fit quality'}, ...
                S.files.outputFile,'Sheet','D fit (2 pop)');
            writecell({'file ID','N','via split value','via component fit','fit quality'}, ...
                S.files.outputFile,'Sheet','slower fraction');
        end

        % Calculate the MSD.
        [MSDresultsArray, meanMSDresults] = ...
            calculateMSD(trackArray, S.MSDanalysis.frameRate, "minTrackLength", S.MSDanalysis.minTrackLength);

        % Fit all individual MSD curves.
        indivMSDfits = fitMSD(MSDresultsArray, S.MSDanalysis.frameRate, 'linFitPoints', S.MSDanalysis.linFitPoints);

        % Remove all data that do not conform with the chosen adj. R² value, have negative diffusion coefficients, were
        % tagged as invalid (row of zeroes) by >fitMSD< or have NaN or Inf values.
        isValid = indivMSDfits(:,1) > 0 & indivMSDfits(:,2) >= S.MSDanalysis.adjR2 & ...
            ~any(isnan(indivMSDfits), 2) & ~any(isinf(indivMSDfits), 2);
        Dlist = indivMSDfits(isValid, 1);
        logDlist = log10(Dlist);

        % Fit the logD data as it if they were normally distributed (strictly speaking, they should follow a lognormal
        % distribution), but this is how it is done across literature.
        statOptions = statset('normfit');
        statOptions.MaxFunEvals = 500;
        pd = fitdist(logDlist, 'Normal', 'Options', statOptions);

        % Try to fit the logD data with a Gaussian mixture model, assuming two populations.
        [pdComp1, pdComp2, fractionComp1, errorCode] = fit2pops(logDlist);

        % Extract the results.
        D1pop = 10^pd.mu;
        slowFractionSplit = sum(logDlist <= S.MSDanalysis.splitValue) / numel(logDlist);
        if isempty(pdComp1) || isempty(pdComp2)
            D2pop1 = 0;
            D2pop2 = 0;
            slowFractionComp = 0;
        else
            D2pop1 = 10.^pdComp1.mu;
            D2pop2 = 10.^pdComp2.mu;
            slowFractionComp = fractionComp1;
        end
        fitScoreList = ["ok","not converged","failed"];
        fitScore = fitScoreList(errorCode + 1);

        % Summarize the mobility analysis results.
        
        % Mean MSD results:
        % (Over)write the lag times to the first column of the sheet. The longest list of lag times will therefore be
        % represented here, while the MSD columns for the individual files can be shorter.
        writecell(num2cell(meanMSDresults(:,1)),S.files.outputFile,'Sheet','mean MSD','Range','A3');

        % The data is then written in column pairs of the MSD and the corresponding SEM* (the asterisk notes that these
        % are SEM values derived from weighted averages).
        writecell(vertcat(num2cell([iFile, iFile]),{'MSD (µm²)','SEM* (µm²)'}),S.files.outputFile,'Sheet','mean MSD',...
            'Range', strcat(matlab.io.spreadsheet.internal.columnLetter(2*iFile),'1'));
        SEM = meanMSDresults(:,3)./sqrt(meanMSDresults(:,4));
        writecell(num2cell([meanMSDresults(:,2),SEM]),S.files.outputFile,'Sheet','mean MSD',...
            'Range', strcat(matlab.io.spreadsheet.internal.columnLetter(2*iFile),'3'));
        
        % Raw D values:
        writecell({iFile; 'D (µm²/s)'}, S.files.outputFile,'Sheet','D lists',...
            'Range',strcat(matlab.io.spreadsheet.internal.columnLetter(iFile+1),'1'));
        writecell(num2cell(Dlist), S.files.outputFile,'Sheet','D lists',...
            'Range',strcat(matlab.io.spreadsheet.internal.columnLetter(iFile+1),'3'));
       
        % Results of single population analysis.
        writecell(num2cell([iFile, numel(logDlist), D1pop]), S.files.outputFile, 'Sheet', 'D fit (1 pop)',...
            'Range', strcat('A',num2str(iFile+1)));
              
        % Results of two populations analysis.
        writecell(horzcat(num2cell([iFile, numel(logDlist), D2pop1, D2pop2]),{fitScore}), ...
            S.files.outputFile, 'Sheet', 'D fit (2 pop)','Range', strcat('A',num2str(iFile+1)));

        % Results of slower/faster fraction analysis (two populations).
        writecell(horzcat(num2cell([iFile, numel(logDlist), slowFractionSplit, slowFractionComp]),{fitScore}), ...
            S.files.outputFile, 'Sheet', 'slower fraction','Range', strcat('A',num2str(iFile+1)));
    end

    %% Cluster analysis
    if S.general.analysisSteps.doAnalysis(4)
        progressDialog.Message = sprintf('%s: Cluster analysis ...', thisFilename);
        progressDialog.Value = (iFile-1+0.8)/nFiles;
        checkCancel(progressDialog);

        if iFile == 1
            % If this is the first file, create the spreadsheet 'cluster results' and write the column headers.
            writecell({'file ID'},S.files.outputFile,'Sheet','cluster results');
        end

        switch S.cluster.algorithm
            case 'Voronoi'
                V = S.cluster.Voronoi;

                % Create the Voronoi diagram and trim the Voronoi cells.
                if V.trimEdges == false
                    trimMethod = 'removeInf';
                else
                    trimMethod = V.trimMethod;
                end
                [DT, voronoiVertices, voronoiCells] = getTrimmedVoronoi(pointList, 'type', trimMethod, ...
                    'boundary', V.trimBoundary);

                % Calculate the area for all Voronoi cells.
                voronoiAreaList = calculateVoronoiArea(voronoiVertices, voronoiCells);

                if strcmp(V.method, 'density')
                    % Calculate the normalized first rank density for the value list.
                    valueList = calculateNormFirstRankDensity(DT, voronoiAreaList);

                    % Get the alpha value.
                    splitValue = V.alpha;

                elseif strcmp(V.method, 'area')
                    % The area list will be used as the value list.
                    valueList = voronoiAreaList;

                    % Get the area split value.
                    splitValue = V.areaThreshold;
                end

                % Group the Voronoi cells to clusters.
                clusterData = buildVoronoiClusters(DT, voronoiVertices, voronoiCells, V.method, valueList, splitValue);

                % Remove all clusters with less points than specified.
                clusterData(clusterData(:,3) < S.cluster.minN,:) = [];

            case 'DBSCAN'
                D = S.cluster.DBSCAN;
                clusterData = buildDBSCANClusters(pointList, D.epsilon, D.minPts);

                % Remove all clusters with less points than specified.
                clusterData(clusterData(:,3) < S.cluster.minN,:) = [];

            case {'NASTIC', 'segNASTIC'}
                radiusFactor = S.cluster.NASTIC.r;

                % Convert the tracks to axis-aligned bounding boxes and assign them to groups.
                AABB = convertTracksToAABB(trackArray, radiusFactor);
                overlapGroups = groupAABBoverlaps(AABB, inf);

                % Remove groups with too few tracks and store the remaining overlap groups. This is done here to save
                % time, especially for the segNASTIC algorithm.
                overlapGroups(cellfun('size', overlapGroups, 2) < S.cluster.minN) = [];

                % Discriminate between the two NASTIC types.
                if contains(S.cluster.algorithm, 'seg')
                    % The overlap threshold determines often one track segment has to overlap with segments of other
                    % tracks to be considered as part of a cluster.
                    clusterData = ...
                        buildSegNASTICclusters(trackArray, overlapGroups, S.cluster.NASTIC.nSegmentOverlaps);
                else
                    clusterData = buildNASTICclusters(trackArray, overlapGroups);
                end

                % Remove all clusters with less points than specified.
                clusterData(clusterData(:,3) < S.cluster.minN,:) = [];
        end

        % The data is written to the sheet 'cluster results' with five columns per file (area, diameter, N, cluster
        % centroid x position, cluster centroid y position). Internally, the cluster areas are in nm², so divide by 1E6
        % to get µm².
        clusterData(:,1) = clusterData(:,1)/1E6;
        
        if isempty(clusterData)
            % If no clusters are found at all, just set the cluster data to all zeros.
            clusterData = zeros(1,5);
        end

        writecell(vertcat(repmat({iFile},1,5), {'area (µm²)','diameter (nm)','N','x position (nm)','y position (nm)'}), ...
            S.files.outputFile, 'Sheet', 'cluster results',...
            'Range', strcat(matlab.io.spreadsheet.internal.columnLetter(5*iFile-3),'1'));
        writecell(num2cell(clusterData),S.files.outputFile,'Sheet','cluster results',...
            'Range', strcat(matlab.io.spreadsheet.internal.columnLetter(5*iFile-3),'3'));
    end
end

if isa(progressDialog,'matlab.ui.dialog.ProgressDialog')
    % Close the progress dialog (if it was not just a dummy structure).
    close(progressDialog);
end

end

%% Subfunctions
function O = formatParamsOutput(S)
O = cell(4,1);

% Localization.
O{1} = cell(13,2);
O{1}{1} = 'Localization';
O{1}(2:7,1) = {'offset (ADU)';'conversion (e-/ADU)';'pixel size (nm)';'filter size (pixel)';...
    'cut-off (photons)';'PSF ROI size (pixel)'};
L = struct2cell(S.localization);
O{1}(2:7,2) = L(2:7);

% Localization filter.
O{1}{9} = 'Localization filter';
O{1}(10:12,1) = {'PSF size (nm)';'photons';'localization precision (nm)'};
histogramLimits = table2cell(S.localization.limits);
O{1}(10:12,2) = cellstr("[" + join(string(histogramLimits),', ') + "]");

% Tracking.
O{2} = cell(10,2);
O{2}{1} = 'Tracking';
O{2}(2:3,1) = {'max. linking distance (nm)';'max. gap closing (frames)'};
T = struct2cell(S.tracking);
O{2}(2:3,2) = T(1:2);

% Tracking filter.
O{2}{5} = 'Tracking filter';
O{2}(6:9,1) = {'number of localizations';'duration (s)';'displacement (nm)';'range (nm)'};
trackingLimits = table2cell(S.tracking.limits);
O{2}(6:9,2) = cellstr("[" + join(string(trackingLimits),', ') + "]");

% MSD analysis.
O{3} = cell(8,2);
O{3}{1} = 'MSD analysis';
O{3}(2:5,1) = {'frame rate (Hz)';'min. track length (frames)';'min. adj. R²';'N points for linear fits'};
M = struct2cell(S.MSDanalysis);
O{3}(2:4,2) = M(1:3);
O{3}(5,2) = {M{4}(end)};
O{3}(6:7,1) = {'first point excluded','immobile/mobile split value (log10 of µm²/s)'};
if isscalar(S.MSDanalysis.linFitPoints)
    O{3}{6,2} = 'no';
else
    O{3}{6,2} = 'yes';
end
O{3}(7,2) = M(5);

% Cluster analysis.
if contains(S.cluster.input, "track")
    clusterInputString = "tracks";
else
    clusterInputString = "localizations";
end

switch S.cluster.algorithm
    case 'Voronoi'
        V = S.cluster.Voronoi;
        O{4} = cell(8,1);
        O{4}{1} = sprintf('Clustering (Voronoi, %s)', V.method);
        O{4}(2:3,1) = {'input'; sprintf('at least N %s per cluster', clusterInputString)};
        O{4}(5:7,1) = {'trim edges';'trim method';'trim boundary'};
        O{4}{2,2} = S.cluster.input;
        O{4}{3,2} = S.cluster.minN;

        if strcmp(V.method, 'density')
            O{4}(4,:) = {'alpha', V.alpha};
        else
            O{4}(4,:) = {'threshold (µm²)', V.areaThreshold};
        end

        logicalString = {'no';'yes'};
        O{4}(5,2) = logicalString(V.trimEdges + 1);
        O{4}{6,2} = V.trimMethod;

        if isnumeric(V.trimBoundary)
            O{4}(7,2) = char("[" + join(string(V.trimBoundary), ', ') + "]");
        else
            choiceString = {'convex hull', 'axis-aligned bounding box'};
            O{4}(7,2) = choiceString(strcmp(V.trimBoundary, 'AABB') + 1);
        end

    case 'DBSCAN'
        D = S.cluster.DBSCAN;
        O{4} = cell(6,1);
        O{4}{1} = 'Clustering (DBSCAN)';
        O{4}(2:5,1) = {'input'; sprintf('at least N %s per cluster', clusterInputString); ...
            'minimum points'; [char(949), ' (nm)']};
        O{4}{2,2} = S.cluster.input;
        O{4}{3,2} = S.cluster.minN;
        O{4}{4,2} = D.minPts;
        O{4}{5,2} = D.epsilon;

    case 'NASTIC'
        O{4} = cell(4,1);
        O{4}{1} = 'Clustering (NASTIC)';
        O{4}(2:3,1) = {'at least N tracks per cluster';'radius factor'};
        O{4}{2,2} = S.cluster.minN;
        O{4}{3,2} = S.cluster.NASTIC.r;

    case 'segNASTIC'
        O{4} = cell(5,1);
        O{4}{1} = 'Clustering (segNASTIC)';
        O{4}(2:4,1) = {'at least N tracks per cluster';'radius factor';'N segment overlaps'};
        O{4}{2,2} = S.cluster.minN;
        O{4}{3,2} = S.cluster.NASTIC.r;
        O{4}{4,2} = S.cluster.NASTIC.nSegmentOverlaps;
end

end

function checkCancel(progressDialog)
    if progressDialog.CancelRequested == true
        % Changes the title of the progress dialog to assure the user that the batch processing is cancelled after the
        % analysis of the current file is completed.
        currentFilename = extractBefore(progressDialog.Message,':');
        progressDialog.Title = sprintf('Stopping after processing %s', currentFilename);
    end
end
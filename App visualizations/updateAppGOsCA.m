function V = updateAppGOsCA(GO, hAxes, V)
% updateAppGOsCA updates graphic objects for data visualization in the Cluster Analysis tab of the OneFlowTraX app.
%
% Syntax:
%   V = updateAppGOsCA(GO, hAxes, V)
%
% Input Arguments:
%   (Required)
%   GO                 Array of graphic objects that can be passed from an app to this function (see notes below for
%                      specifics). Generally, they should contain all objects that might change in |hAxes| or that are
%                      directly influenced by calculated values in |hAxes| (like colorbars or data limit edit fields and
%                      their labels that have their values or visibility changed).
%                      (:,1) graphic objects array
%
%   hAxes              Handle to the axes (or array of axes) populated by the objects in |GO|.
%
%   V                  Structure holding all data and information to update |GO|. V(1) holds the old, and V(2) the new
%                      parameters.
%                      (2,1) struct
%
% Output Arguments:
%   V                  Parameter structure, reset to the new values (V(1) = V(2)).
%                      (2,1) struct
%
% Other required m-files: [mList,pList] = matlab.codetools.requiredFilesAndProducts('fname.m')
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Many apps have axes with multiple layers of information (images, lines, text, points) that will be plotted in certain
% combinations, depending on what display options were chosen by the user. Processing this information often bloats the
% app code, so it is outsourced here. This function is very specific (order of graphic objects, fields of the parameter
% structure, handle names of the axes and UI components) and must be developed along the corresponding app.
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
% Initial release: 2023-03-17
% Last revision: 2023-03-21

%% Function argument validation
arguments
    GO (:,1)
    hAxes (:,1) {mustBeA(hAxes, 'matlab.graphics.axis.Axes')}
    V (2,1) struct
end

%% Main

% Compare the first (old) and second (possibly changed) elements of |VisParams| and extract the
% field names of the updated values.
hasChanged = cellfun(@(A,B) ~isequal(A,B), struct2cell(V(1)), struct2cell(V(2)));
fields = string(fieldnames(V));
taskList = fields(hasChanged);

% After the task list has been created, store the old values in |oldV|, while |V| will now only
% hold the new values.
oldV = V(1);
V = V(2);

% Store the visibility status of all objects.
isVisible = cellfun(@logical, get(GO, 'Visible'));

% Process the task list from start to finish until no more tasks are left. This can also be cut short by
% programmatically deleting all remaining tasks in the list. Equally, new tasks can be added at any position.
while ~isempty(taskList)
    % Extract the first task and remove it from the list.
    currentTask = taskList(1);
    taskList(1) = [];

    switch currentTask
        case {'pointList', 'centroidList', 'trackArray'}
        % The point list changes when a new file from the list was selected.
        isVisible(1:9) = false;      

        % Set the axes (back) to their maximum limits.
        hAxes.XLim = V.axesLimits(1:2);
        hAxes.YLim = V.axesLimits(3:4);
        
        % Remove 'pointList', 'centroidList' and 'trackArray' from the task list, as they will have equal outcomes.
        % Although changing the tracking parameters would not change anything when the localizations are currently
        % displayed, it's better to "clean the slate". Recalculate the visual data by adding 'baseDataStyle' to the
        % list, which will also trigger 'colorscheme'. Also add 'colorMap' so the colormap is updated.
        taskList = setdiff(taskList, ["pointList";"centroidList";"trackArray"], 'stable');
        taskList = unique(["baseDataStyle"; "colorMap"; taskList(:)], 'stable');

        case 'showClusterAlg'
            switch V.clusterAlg
                case 'none'
                    isVisible([2 3 6]) = [0 0 0];
                case 'Voronoi'
                    isVisible([2 3 6]) = [1 1 0];
                case 'DBSCAN'
                     isVisible([2 3 6]) = [0 0 1];
                case 'NASTIC'
                    %%%% later
            end

            % Put 'baseDataStyle' in front to plot the correct base data.
            taskList = unique(["baseDataStyle"; taskList(:)], 'stable');
            
        case 'voronoiPoints'
            % Show the Voronoi cells (points are processed with 'baseDataStyle').
            isVisible([2 3]) = [1 1];
     
            % Remove all Voronoi tasks from the list. As there is a new set of Voronoi cells, add
            % 'baseDataStyle', which will set the base data for the Voronoi cells (track centroids
            % or localizations) and also trigger 'colorScheme' to update the data (density or area)
            % that the color map can index into.
            taskList = setdiff(taskList, ["voronoiPoints";"voronoiVertices";"voronoiCells";...
                "voronoiValueList";"showClusterAlg"], 'stable');
            taskList = unique(["baseDataStyle"; taskList(:)], 'stable');

        case 'DBSCANpoints'
            % Hide the Voronoi cells, if present.
            isVisible([2 3]) = [0 0];

            % Remove all DBSCAN tasks from the list. As there is a new DBSCAN points list, add
            % 'baseDataStyle', which will set the base data DBSCAN (track centroids or
            % localizations) and also trigger 'colorScheme' to update the point values that the
            % color map can index into. Also remove 'showClusterAlg' as the algorithm-specific plot
            % is always switched on for new/recalculated data. 
            taskList = setdiff(taskList, ["DBSCANpoints";"DBSCANclusterIndexList";...
                "DBSCANisCore";"showClusterAlg"], 'stable');
            taskList = unique(["baseDataStyle"; taskList(:)], 'stable');

        case 'overlapGroupsNASTIC'
            % Hide the Voronoi cells, if present.
            isVisible([2 3]) = [0 0];

            % Add 'colorScheme' to the task list to also update the coloring of the tracks.
            taskList = unique([taskList(:); "colorScheme"], 'stable');

        case 'clusterBoundaryArray'
            if isempty(V.clusterBoundaryArray)
                isVisible([4 5]) = [0 0];
            else
                isVisible([4 5]) = [1 1];

                % Calculate the plot data.
                if strcmp(V.colorScheme, 'Voronoi (density)')
                    clusterVertices = V.voronoiPoints;
                elseif strcmp(V.colorScheme, 'Voronoi (area)')
                    clusterVertices = V.voronoiVertices;
                elseif contains(V.colorScheme, 'DBSCAN')
                    clusterVertices = V.DBSCANpoints;
                elseif contains(V.colorScheme, 'NASTIC')
                    trackArray = V.trackArray(cellfun("size",V.trackArray,1) >= 3);
                    clusterVertices = ...
                        cell2mat(cellfun(@(x) x(:,2:3), trackArray, 'UniformOutput', false));
                end

                [clusterFaces, clusterEdges] = preparePatchData(clusterVertices, V.clusterBoundaryArray);
                set(GO(4), {'Faces','Vertices'}, {clusterFaces, clusterVertices});
                set(GO(5), {'XData','YData'}, {clusterEdges(:,1), clusterEdges(:,2)});

                % Add 'clusterAreaColor' at the top of the task list to force the subsequent execution of the cluster
                % style update.
                taskList = unique(["clusterAreaColor"; taskList(:)], 'stable');
            end

        case {'clusterBoundariesLineWidth', 'clusterBoundariesColor', ...
                'clusterAreaIsSemiTransparent','clusterAreaColor'}
            set(GO(4), 'FaceColor', V.clusterAreaColor, ...
                'FaceAlpha', 1 - 0.3*V.clusterAreaIsSemiTransparent);
            set(GO(5), 'Color', V.clusterBoundariesColor,  'LineWidth', V.clusterBoundariesLineWidth);

            % Cluster area/line styles are set all at once. They can then all be removed from
            % the task list.
            taskList = setdiff(taskList, ...
                ["clusterBoundariesLineWidth"; "clusterBoundariesColor"; ...
                "clusterAreaIsSemiTransparent"; "clusterAreaColor"], 'stable');

        case 'baseDataStyle'
            switch V.baseDataStyle
                case {'localizations', 'track centroids'}
                    isVisible([1 6]) = [0 1];
                case 'tracks'
                    isVisible([1 6]) = [1 0];
                otherwise
                    % Empty.
            end

            % Add 'colorScheme' at the top of the task list to force its subsequent execution.
            taskList = unique(["colorScheme"; taskList(:)], 'stable');

        case 'colorScheme'
            % This case sets the data for the plots and also handles the values the plot colors are
            % based on. With these calculated color values (CData), the limits, values and tooltips
            % of the limit edit fields are set. The colorbar visibility and the plot color limits
            % are also handled here. 

            % By default, the colorbar is visible (exceptions are handled in the cases below).
            isVisible(9) = true;

            % Adjust the 'random' case string to contain the desired number of cycles.
            if contains(V.colorScheme, 'random')
                nCycleCols = str2double(extract(V.colorMap, digitsPattern));
                colorScheme = "random" + num2str(nCycleCols);
            else
                colorScheme = V.colorScheme;
            end

            % Calculate the plot data according to the chosen color scheme. Also adjust the limits
            % and default values of the lower/upper limit edit fields and the unit label. Populate
            % the tooltips for the edit fields.
            if contains(V.baseDataStyle, 'localizations')
                pointList = V.pointList;
            else
                % For 'track centroids' (and as a dummy for 'tracks').
                pointList = V.centroidList;
            end

            switch colorScheme
                case 'unicolored'
                    set(GO(6), {'XData','YData'}, {pointList(:,1), pointList(:,2)});
                    isVisible(7) = false;

                case {'random4', 'random6', 'random8', 'random10'}
                    % Handle this entirely in the 'colormap' case, as the data needs to be
                    % recalculated for each color set (consequently, the limits are also set there).
                    taskList = unique(["colorMap"; taskList(:)], 'stable');   

                case 'immobile/mobile (logD)'
                    [X,Y,Z,C,colorValues] = prepareTracksPlotData(V.trackArray, ...
                        'colorScheme', colorScheme, 'MSDparams', V.MSDparams);
                    set(GO(1), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    set(GO(6), {'XData','YData','CData','MarkerFaceColor'}, ...
                        {pointList(:,1), pointList(:,2), colorValues, 'flat'});
                    [LL, UL] = bounds(double(colorValues));
                    GO(10).Limits = [LL, UL];
                    [GO(10).Value, V.dataLL] = deal(LL);
                    GO(12).Limits = [LL, UL];

                    % Set the default threshold value to the median of the data and round it to 2 digits to the right of
                    % the decimal point.
                    medianThreshold = round(median(double(colorValues), 'omitnan'), 2);
                    [GO(12).Value, V.dataUL] = deal(medianThreshold);
                    hAxes.CLim = 1000*[-eps, +eps] + GO(12).Value;
                    unitStr = "log" + char([0x2081, 0x2080]) + " of µm²/s";
                    GO(12).ValueDisplayFormat = "%.2f " + unitStr;
                    GO(12).Tooltip = ...
                        sprintf('Minimum: %.2f %s\nMaximum: %.2f %s', LL, unitStr, UL, unitStr);

                case {'mobility (D)', 'mobility (logD)', 'duration', 'displacement' 'range'}
                    [X,Y,Z,C,colorValues] = prepareTracksPlotData(V.trackArray, ...
                        'colorScheme', colorScheme, 'MSDparams', V.MSDparams);
                    set(GO(1), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    set(GO(6), {'XData','YData','CData','MarkerFaceColor'}, ...
                        {pointList(:,1), pointList(:,2), colorValues, 'flat'});

                    schemeOptions = ["mobility (D)", "mobility (logD)", "duration", "displacement", "range"];
                    unitStrArray = ["µm²/s", "log" + char([0x2081, 0x2080]) + " of µm²/s", "seconds", "nm", "nm"];

                    % Set the display format (i.e. precision) here for the different scheme options. Use capital "E" for
                    % scientific notation.
                    numFormatArray = ["%0.1E", "%.2f", "%.2f", "%.0f", "%.0f"];
                    unitStr = unitStrArray(ismember(schemeOptions, colorScheme));
                    numFormat = numFormatArray(ismember(schemeOptions, colorScheme));

                    % Determine how to round the values in the edit fields (otherwise they would be displayed with
                    % maximum precision the first time the user clicks into a field).
                    if contains(numFormat,"E")
                        % Values with scientific notation are a bit special as they must still be entered as decimal
                        % numbers. Therefore, just choose the smallest number as a template.
                        roundingDigits = -ceil(log10(min(colorValues,[],"all"))) + 1;
                    else
                        roundingDigits = str2double(extract(numFormat,digitsPattern));
                    end

                    % Set the limits and update the edit fields.         
                    [LL, UL] = bounds(double(colorValues));
                    LL = round(LL,roundingDigits);
                    UL = round(UL,roundingDigits);
                    GO(10).Limits = [LL, UL];
                    [GO(10).Value, V.dataLL] = deal(LL);
                    GO(12).Limits = [LL, UL];
                    [GO(12).Value, V.dataUL] = deal(UL);
                    hAxes.CLim = [LL, UL];                 
                   
                    set(GO([10,12]), 'ValueDisplayFormat', numFormat + " " + unitStr);
                    GO(10).Tooltip = sprintf("Minimum: " + numFormat + " %s", LL, unitStr);
                    GO(12).Tooltip = sprintf("Maximum: " + numFormat + " %s", UL, unitStr);

                case {'Voronoi (density)', 'Voronoi (area)'}
                    % The correct seed points (based on the localizations or centroids, possibly
                    % altered by edge trimming) are stored in DT. Paint them black.
                    pointList = V.voronoiPoints;
                    set(GO(6), {'XData','YData','MarkerFaceColor'}, ...
                        {pointList(:,1), pointList(:,2), [0 0 0]});
                    
                    % Plot the Voronoi cells, setting their values with the stored list (density or
                    % area). Also plot the edges (in black).
                    [voronoiFaces, voronoiEdges] = ...
                        preparePatchData(V.voronoiVertices, V.voronoiCells);
                    set(GO(2), {'Faces','Vertices','FaceVertexCData'}, ...
                        {voronoiFaces, V.voronoiVertices, V.voronoiValueList});
                    set(GO(3), {'XData','YData','Color'}, ...
                        {voronoiEdges(:,1), voronoiEdges(:,2), [0 0 0]});

                    % Update the limits.
                    [LL, UL] = bounds(V.voronoiValueList);
                    GO(10).Limits = [LL, UL];
                    [GO(10).Value, V.dataLL] = deal(LL);
                    GO(12).Limits = [LL, UL];

                     % Set the upper limit to the median of the data.
                    [GO(12).Value, V.dataUL] = deal(median(V.voronoiValueList));
                    hAxes.CLim = [LL, GO(12).Value];

                    if contains(colorScheme, 'density')
                        set(GO([10 12]), 'ValueDisplayFormat', "%.2f");
                        GO(10).Tooltip = ...
                            sprintf("Minimum (normalized first-rank density):\n%.2f", LL);
                         GO(12).Tooltip = ...
                            sprintf("Maximum (normalized first-rank density):\n%.2f", UL);
                    else
                        set(GO([10 12]), 'ValueDisplayFormat', "%.2f");
                        GO(10).Tooltip = sprintf("Minimum: %.0f nm²", LL);
                        GO(12).Tooltip = sprintf("Maximum: %.0f nm²", UL);
                    end

                case {'DBSCAN (noise/cluster)', 'DBSCAN (noise/core/border)'}
                    colorScheme = extractBetween(string(colorScheme),"(",")");
                    colorValues = prepareDBSCANplotData(V.DBSCANclusterIndexList, ...
                        V.DBSCANisCore, colorScheme);
                    set(GO(6), {'XData','YData','CData','MarkerFaceColor'}, ...
                        {V.DBSCANpoints(:,1), V.DBSCANpoints(:,2), colorValues, 'flat'});
                    [LL, UL] = bounds(colorValues);
                    GO(10).Limits = [LL, UL];
                    [GO(10).Value, V.dataLL] = deal(LL);
                    GO(12).Limits = [LL, UL];
                    [GO(12).Value, V.dataUL] = deal(UL);
                    hAxes.CLim = [LL, UL];

                case {'NASTIC', 'segNASTIC'}
                    % The first/second color of the chosen color duo will be used to paint all
                    % tracks that are outside/inside clusters.
                    % Get the all track IDs that are stored in the overlap groups.
                    clusterTrackIDs = horzcat(V.overlapGroupsNASTIC{:}); 

                    % These track IDs refer to the set of tracks that consist of at least 3
                    % localizations, so count them to generate a binary mask of tracks that are
                    % inside clusters.
                    validTrackArray = V.trackArray(cellfun("size", V.trackArray, 1) >=3);
                    binaryMask = false(numel(validTrackArray), 1);
                    binaryMask(clusterTrackIDs) = true;

                    [X,Y,Z,C] = prepareTracksPlotData(validTrackArray, ...
                        'colorScheme', 'binary', 'binaryMask', binaryMask);
                    set(GO(1), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    GO(10).Limits = [1, 2];
                    [GO(10).Value, V.dataLL] = deal(1);
                    GO(12).Limits = [1, 2];
                    [GO(12).Value, V.dataUL] = deal(2);
                    hAxes.CLim = [1, 2];

                otherwise
                    % Empty.
            end

            % As the limits are (re)set here, delete any tasks that would (inadvertently) change the
            % limits afterwards.
            taskList = setdiff(taskList, ["dataLL";"dataUL"], 'stable');
           
        case 'colorMap'
            % This case changes the plot colors according to the chosen option (but does not
            % recalculate the CData). The exception are the 'cycle' colormaps that are handled here
            % and not in 'colorScheme'.

            if strcmp(V.colorScheme, "unicolored")
                GO(6).MarkerFaceColor = V.colors{V.colorMap, :};
        
            elseif contains(V.colorMap, "cycle")
                nCycleCols = str2double(extract(V.colorMap, digitsPattern));
                colorScheme = "random" + num2str(nCycleCols);
                if contains(V.colorScheme, 'DBSCAN')
                    colorValues = prepareDBSCANplotData(V.DBSCANclusterIndexList, ...
                        V.DBSCANisCore, colorScheme);
                    set(GO(6), {'XData','YData','CData','MarkerFaceColor'}, ...
                        {V.DBSCANpoints(:,1), V.DBSCANpoints(:,2), colorValues, 'flat'});
                    [LL, UL] = bounds(colorValues);
                    GO(10).Limits = [LL, UL];
                    [GO(10).Value, V.dataLL] = deal(LL);
                    GO(12).Limits = [LL, UL];
                    [GO(12).Value, V.dataUL] = deal(UL);
                    hAxes.CLim = [LL, UL];

                    % Add gray as the first color for the noise points.
                    hAxes.Colormap = [0.5 0.5 0.5; getSMLMcolorMap(V.colorMap)];

                else
                    [X,Y,Z,C] = prepareTracksPlotData(V.trackArray, 'colorScheme', colorScheme);
                    set(GO(1), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    [LL, UL] = bounds(GO(1).CData, "all");
                    GO(10).Limits = [LL, UL];
                    [GO(10).Value, V.dataLL] = deal(LL);
                    GO(12).Limits = [LL, UL];
                    [GO(12).Value, V.dataUL] = deal(UL);
                    hAxes.CLim = [LL, UL];
                    hAxes.Colormap = getSMLMcolorMap(V.colorMap);
                end
            else
                hAxes.Colormap = getSMLMcolorMap(V.colorMap);
            end
 
        case 'pointSize'
            if V.pointSize == 0
                isVisible(6) = false;
                % Just set to invisible and process the next task.
                continue
            end

            isVisible(6) = true;
            GO(6).SizeData = V.pointSize;


        case 'gLineWidth'
             if V.gLineWidth == 0
                 isVisible([1 3]) = false;
                 % Just set to invisible and process the next task.
                 continue
             else
                 isVisible([1 3]) = true;
             end
             
             if contains(V.colorScheme, "Voronoi")
                 GO(3).LineWidth = V.gLineWidth;
             elseif strcmp(V.baseDataStyle, "tracks")
                 GO(1).LineWidth = V.gLineWidth;
             end

        case 'scaleBar'
            if contains(V.scaleBar,"full")
                isVisible(7:8) = [true, true];
                [rectanglePos, labelStr, labelPos] = ...
                    placeScaleBar([hAxes.XLim, hAxes.YLim], -9);
                GO(7).Position = rectanglePos;
                set(GO(8),{'String','Position'},{labelStr,labelPos});

            elseif contains(V.scaleBar,"bar")
                isVisible(7:8) = [true, false];
                rectanglePos = placeScaleBar([hAxes.XLim, hAxes.YLim], -9);
                GO(7).Position = rectanglePos;

            else
                isVisible(7:8) = [false, false];
            end

            [GO(7).FaceColor, GO(8).Color] = deal(ones(1,3)*contains(V.scaleBar, "(w"));

        case 'isWhiteBG'
            hAxes.Color = ones(1,3)*V.isWhiteBG;

        case 'hideBaseData'
            switch V.baseDataStyle
                case {'localizations', 'track centroids'}
                    isVisible([1 6]) = [0, ~V.hideBaseData];
                case 'tracks'
                    isVisible([1 6]) = [~V.hideBaseData, 0];
                otherwise
                    % Empty.
            end

        case {'dataLL', 'dataUL'}
            % Force the limits to be in the correct order and overwrite them in V (to ensure the
            % same sorting order).
            [LL, UL] = bounds([V.dataLL, V.dataUL]);
            V.dataLL = LL;
            V.dataUL = UL;

            if LL == UL
                % Equal limits are not allowed for CLim, so add 1000*>eps< to the upper limit
                % (this will be imperceptible in the plot).
                UL = UL + 1000*eps;
            end

            if strcmp(V.colorScheme, 'immobile/mobile (logD)')
                LL = UL - 1000*eps;
                UL = UL + 1000*eps;
            end
            hAxes.CLim = [LL, UL];
    end
end

% All new parameters of |V| have been processed so that old and new values can now be set the same.
V = [V;V];

% Finally, update the visibility of all objects.
set(GO, {'Visible'}, num2cell(isVisible)); drawnow

end
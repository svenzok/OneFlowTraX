function V = updateAppGOsTracksVis(GO, hAxes, V)
% updateAppGOsTracksVis updates graphic objects to visualize tracks in the Tracking tab of the OneFlowTraX app.
%
% Syntax:
%   V = updateAppGOsTracksVis(GO, hAxes, V)
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
% Initial release: 2023-07-05
% Last revision: 2023-07-05

%% Function argument validation
arguments
    GO (:,1)
    hAxes (:,1) {mustBeA(hAxes, 'matlab.graphics.axis.Axes')}
    V (2,1) struct
end

%% Main

% Compare the first (old) and second (possibly changed) elements of |V| and extract the
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
        case 'isSingleFileSelected'
            if V.isSingleFileSelected
                % When a single file is newly selected, set all corresponding plots to visible.
                isVisible(1) = false;
                isVisible(2:13) = true;

                % Add tasks that control the visibility of selected elements in the GUI. 
                taskList = unique(["showTIFF"; "scaleBar"; taskList(:)], 'stable');

            else
                % Hide all plot data, just display the text message.
                isVisible(1) = true;
                isVisible(2:13) = false;

                % Remove all other tasks to return to the calling app.
                taskList = string.empty;
            end

        case 'trackArray'
        % The track array changes when a new file from the list is selected or when new parameters lead to an updated
        % track array of the same file.
        
        % Set the axes of the left plot (back) to their maximum limits.
        hAxes(1).XLim = V.axesLimits(1:2);
        hAxes(1).YLim = V.axesLimits(3:4);

        % The initial, 10x magnified section shown in the right plot originates from center of the left plot.
        M = 10;
        maxX = V.axesLimits(2);
        maxY = V.axesLimits(4);
        hAxes(2).XLim = [maxX*(1-1/M)/2, maxX*(1+1/M)/2];
        hAxes(2).YLim = [maxY*(1-1/M)/2, maxY*(1+1/M)/2];
        
        % Add 'changedZoom' to the task list. Also recalculate the visual data by adding 'colorScheme' to the list and
        % add 'colorMap' so the colormap is updated.
        taskList = unique(["changedZoom"; "colorScheme"; "colorMap"; taskList(:)], 'stable');

        case 'colorScheme'
            % This case sets the data for the plots and also handles the values the plot colors are based on. With these
            % calculated color values (CData), the limits, values and tooltips of the limit edit fields are set. The
            % colorbar visibility and the plot color limits are also handled here.

            % By default, the colorbar is visible (exceptions are handled in the cases below).
            isVisible(13) = true;

            % Adjust the 'random' case string to contain the desired number of cycles.
            if contains(V.colorScheme, 'random')
                nCycleCols = str2double(extract(V.colorMap, digitsPattern));
                colorScheme = "random" + num2str(nCycleCols);
            else
                colorScheme = V.colorScheme;
            end

            switch colorScheme
                case {'random4', 'random6', 'random8', 'random10'}
                    % Handle this entirely in the 'colormap' case, as the data needs to be recalculated for each color
                    % set (consequently, the limits are also set there).
                    taskList = unique(["colorMap"; taskList(:)], 'stable');   

                case 'immobile/mobile (logD)'
                    [X,Y,Z,C,colorValues] = prepareTracksPlotData(V.trackArray, ...
                        'colorScheme', colorScheme, 'MSDparams', V.MSDparams);
                    set(GO(4), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    set(GO(5), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    [LL, UL] = bounds(double(colorValues));
                    GO(14).Limits = [LL, UL];
                    [GO(14).Value, V.dataLL] = deal(LL);
                    GO(16).Limits = [LL, UL];

                    % Set the default threshold value to the median of the data and round it to 2 digits to the right of
                    % the decimal point.
                    medianThreshold = round(median(double(colorValues), 'omitnan'), 2);
                    [GO(16).Value, V.dataUL] = deal(medianThreshold);
                    hAxes(1).CLim = 1000*[-eps, +eps] + GO(16).Value;
                    hAxes(2).CLim = 1000*[-eps, +eps] + GO(16).Value;
                    unitStr = "log" + char([0x2081, 0x2080]) + " of µm²/s";
                    GO(16).ValueDisplayFormat = "%.2f " + unitStr;
                    GO(16).Tooltip = sprintf('Minimum: %.2f %s\nMaximum: %.2f %s', LL, unitStr, UL, unitStr);

                case {'mobility (D)', 'mobility (logD)', 'duration', 'displacement' 'range'}
                    [X,Y,Z,C,colorValues] = prepareTracksPlotData(V.trackArray, ...
                        'colorScheme', colorScheme, 'MSDparams', V.MSDparams);
                    set(GO(4), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                    set(GO(5), {'XData','YData','ZData','CData'}, {X,Y,Z,C});

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
                    GO(14).Limits = [LL, UL];
                    [GO(14).Value, V.dataLL] = deal(LL);
                    GO(16).Limits = [LL, UL];
                    [GO(16).Value, V.dataUL] = deal(UL);
                    set(hAxes, 'CLim', [LL, UL]);
                    set(GO([14,16]), 'ValueDisplayFormat', numFormat + " " + unitStr);
                    GO(14).Tooltip = sprintf("Minimum: " + numFormat + " %s", LL, unitStr);
                    GO(16).Tooltip = sprintf("Maximum: " + numFormat + " %s", UL, unitStr);

                otherwise
                    % Empty.
            end

            % As the limits are (re)set here, delete any tasks that would (inadvertently) change the
            % limits afterwards.
            taskList = setdiff(taskList, ["dataLL";"dataUL"], 'stable');
           
        case 'colorMap'
            % This case changes the plot colors according to the chosen option (but does not recalculate the CData). The
            % exception are the 'cycle' colormaps that are handled here and not in 'colorScheme'.
            if contains(V.colorMap, "cycle")
                nCycleCols = str2double(extract(V.colorMap, digitsPattern));
                colorScheme = "random" + num2str(nCycleCols);
                [X,Y,Z,C] = prepareTracksPlotData(V.trackArray, 'colorScheme', colorScheme);
                set(GO(4), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                set(GO(5), {'XData','YData','ZData','CData'}, {X,Y,Z,C});
                [LL, UL] = bounds(GO(4).CData, "all");
                GO(14).Limits = [LL, UL];
                [GO(14).Value, V.dataLL] = deal(LL);
                GO(16).Limits = [LL, UL];
                [GO(16).Value, V.dataUL] = deal(UL);
                set(hAxes, 'CLim', [LL, UL]);
            end
            set(hAxes, 'Colormap', getSMLMcolorMap(V.colorMap));

        case 'lineWidth'
                 GO(5).LineWidth = V.lineWidth;
                 GO(4).LineWidth = V.lineWidth / 3;

        case 'scaleBar'
            if contains(V.scaleBar,"full")
                isVisible(8:11) = true;
                [rectanglePos, labelStr, labelPos] = placeScaleBar([hAxes(1).XLim, hAxes(1).YLim], -9);  
                GO(8).Position = rectanglePos;
                set(GO(10),{'String','Position'},{labelStr,labelPos});
                [rectanglePos, labelStr, labelPos] = placeScaleBar([hAxes(2).XLim, hAxes(2).YLim], -9);  
                GO(9).Position = rectanglePos;
                set(GO(11),{'String','Position'},{labelStr,labelPos});

            elseif contains(V.scaleBar,"bar")
                isVisible(8:11) = [true, true, false, false];
                GO(8).Position = placeScaleBar([hAxes(1).XLim, hAxes(1).YLim], -9);
                GO(9).Position = placeScaleBar([hAxes(2).XLim, hAxes(2).YLim], -9);

            else
                isVisible(8:11) = false;
            end

            [GO(8).FaceColor, GO(9).FaceColor, GO(10).Color, GO(11).Color] = deal(ones(1,3)*contains(V.scaleBar, "(w"));

        case 'isWhiteBG'
            % Set the background of the axes to the chosen color, use the opposite color for the text message and the
            % ROI square in the foreground.
            colorBG = ones(1,3)*V.isWhiteBG;
            colorFG = abs(colorBG-1);
            set(hAxes, 'Color', colorBG);
            GO(1).Color = colorFG;
            GO(12).EdgeColor = colorFG;

            % Add 'showUnmasked' to the task list to adjust the mask if necessary.
            taskList = unique(["showUnmasked"; taskList(:)], 'stable');

        case {'dataLL', 'dataUL'}
            % Force the limits to be in the correct order and overwrite them in V (to ensure the same sorting order).
            [LL, UL] = bounds([V.dataLL, V.dataUL]);
            V.dataLL = LL;
            V.dataUL = UL;

            if LL == UL
                % Equal limits are not allowed for CLim, so add 1000*>eps< to the upper limit (this will be
                % imperceptible in the plots).
                UL = UL + 1000*eps;
            end

            if strcmp(V.colorScheme, 'immobile/mobile (logD)')
                LL = UL - 1000*eps;
                UL = UL + 1000*eps;
            end
            set(hAxes, 'CLim', [LL, UL]);

        case 'changedZoom'
            % Changing the position or zoom of the right plot also has to change the position of the ROI square in the
            % left plot.
            GO(12).Position = [hAxes(2).XLim(1), hAxes(2).YLim(1) , diff(hAxes(2).XLim), diff(hAxes(2).YLim)];

            % Reset the 'changedZoom' info.
            V.changedZoom = false;

            % Add 'scaleBar' to the taskList to update it.
            taskList = unique(["scaleBar"; taskList(:)], 'stable');

        case 'showTIFF'
            isVisible(2:3) = V.showTIFF;

            % Match the stored maximum intensity projection to the axes limits. Image pixel coordinates are centered in
            % the middle of pixels, so this must be corrected to match the other data.
            [pixY, pixX] = size(V.MIP,[1 2]);
            xRange = [V.axesLimits(1) + 0.5*V.axesLimits(2)/pixX, V.axesLimits(2) - 0.5*V.axesLimits(2)/pixX];
            yRange = [V.axesLimits(3) + 0.5*V.axesLimits(4)/pixY, V.axesLimits(4) - 0.5*V.axesLimits(4)/pixY];
            set(GO(2),'XData', xRange, 'YData', yRange, 'CData', V.MIP);
            set(GO(3),'XData', xRange, 'YData', yRange, 'CData', V.MIP);

            % This option changes the mask that is applied to the images, so add that case to the task list.
            taskList = unique(["showUnmasked"; taskList(:)], 'stable');

        case 'showUnmasked'
            isVisible(6:7) = ~V.showUnmasked;

            % The masking images in both axes should just cover up the tracks that are outside the mask. This is done
            % here by overlaying them with background color images with accordingly set alpha values (or the MIP if the
            % TIFF underlay is active). Use the same range as the underlay TIFF files to correctly display the pixels.
            if V.showTIFF
                overlayMask = V.MIP;
            else
                overlayMask = V.MIP*0 + V.isWhiteBG*255;
            end
            set(GO(6),'XData',GO(2).XData,'YData',GO(2).YData,'CData',overlayMask,'AlphaData',~V.maskIm);
            set(GO(7),'XData',GO(3).XData,'YData',GO(3).YData,'CData',overlayMask,'AlphaData',~V.maskIm);
    end
end

% All new parameters of |V| have been processed so that old and new values can now be set the same.
V = [V;V];

% Finally, update the visibility of all objects.
set(GO, {'Visible'}, num2cell(isVisible)); drawnow

end
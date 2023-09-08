function [X,Y,Z,C,colorValues] = prepareTracksPlotData(trackArray, Options)
% prepareTracksPlotData converts particle tracks to patch object data for plotting.
%
% Syntax:
%   [X,Y,Z,C] = prepareTracksPlotData(trackArray)
%   [X,Y,Z,C] = prepareTracksPlotData(trackArray, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   trackArray         Particle trajectories as a cell array. Each cell (track) has at least three
%                      columns, in the order [frame, x-coordinate, y-coordinate] and at least two
%                      rows (localizations). The coordinates must be given in nm. Additional
%                      requirements depend on the |colorScheme| options (see below).
%                      (:,1) cell
%
% Name-Value Pair Input Arguments:
%   (Optional)
%   BoolSelect         Logical array to select which tracks will be plotted. Must have the same
%                      length as |trackArray| (default: plot all tracks).
%                      (:,1) logical
%
%   colorScheme        Defines which property of the individual particle track determines its color.
%                      'binary': each track is assigned a value of either 1 or 2 (needs the
%                                |binaryMask| logical vector).
%                      'randomX': each track gets a random value from 1 to X, with X = 4,6,8 or 10.
%                      'N': number of points the track consists of.
%                      'duration': duration (incl. gaps) of the track in seconds. Needs the 
%                                  |MSDparams| structure for the frame rate.
%                      'PSF size': average PSF size (PSF sizes in nm for each localization must be
%                                  provided in column 4 of each track).
%                      'photons': average number of photons (column 5).
%                      'precision': average CRLB precision (column 7, in nm).
%                      'displacement': distance between the first and last point in nm.
%                      'range': range estimation by encompassing all points of a track with a
%                               rectangle and calculating the length of its diagonal in nm.
%                      'mobility (D)': diffusion coefficients.
%                      'mobility (logD)': logarithmized diffusion coefficients.
%                      'immobile/mobile (logD)': Tracks are colored based on whether their
%                                                respective (logarithmized) diffusion coefficient is
%                                                below or above a threshold value (for a colormap
%                                                with only two colors ).
%                      The last three mobility plot options need at least three localizations per
%                      track and the |MSDparams| structure (default: 'random6').
%                      (1,:) char | (1,1) string
%                      {'random4', 'random6', 'random8', 'random10', 'N', 'duration', ...
%                       'PSF size', 'photons', 'precision', 'displacement', 'range', ...
%                       'mobility (D)', 'mobility (logD)', 'immobile/mobile'}
%
%   binaryMask         Logical vector for the plot option 'binary', set to true for all tracks that
%                      will have the value 2. If this input is missing, all tracks will be assigned
%                      the value 1 (default: all false, same length as |trackArray|).
%                      (:,1) logical
%
%   MSDparams          Structure that is required for all mobility plot options in |colorScheme|, 
%                      with fields:
%                      .frameRate: camera acquisition frame rate in Hz
%                      .minTrackLength: minimum number of localizations
%                      .adjR2: minimum goodness of fit
%                      .linFitPoints: datapoints for the linear fit
%                      (1,1) struct
%
% Output Arguments:
%   [X,Y,Z,C]          Vertex coordinate and color data for a patch object with the syntax:
%                      patch(X,Y,Z,C,'FaceColor','none','EdgeColor','flat')
%                      (1,:) double, all outputs
%
%   colorValues        Color data for each track (e.g., when the track centroids are plotted).
%                      (:,1) double
%
% Other required files: linearRegression.m, calculateMSD.m, fitMSD.m
% Subfunctions: mustBeTrackArray, checkOptionsCompatibility
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.4
%
% Notes:
% As there can be thousands of tracks for one plot, the use of >plot< in a loop is unreasonable, as
% this would create far to many graphics objects to handle.
% A fast method is to concatenate all data with NaN rows as separators, which only produces one
% graphics object, e.g. a single line with a gap between individual trajectories. With >plot< or
% >line<, such an object could only exist in one color.
% However, >patch< can be used by setting the z-coordinates to zero and defining a color array
% that paints each vertex (given by the x- and y-coordinates) in a certain color. The edges between
% the vertices are then equal to the particle tracks and will assume the color assigned to their
% respective vertices. Each color is represented by a real number in the color array. The smallest
% and highest values correspond to the first and last row of the chosen colormap, respectively. The
% intermediate values are mapped linearly to the intermediate rows in the colormap.
% Example 1: if the vertices are colored according to vector C = [1 6 6 6 3] and the "gray"
% colormap, then the first vertex would be black, the next three white and the last one gray. C can
% be constructed with >repelem< to provide the same color for one trajectory.
% Example 2: with NaN values as separators, the vector
% C = [5 5 5 5 NaN 15 15 15 NaN 10 10 10 10 10 10] with a colormap ranging from black to green would
% color the respective vertices black (5 5 5 5), green (15 15 15) and dark green (6x 10). The
% connecting edges assume the bounding vertex colors, and will have gaps (as they ignore NaN
% values). In the plot (vertices are not colored), this will look like three separately colored
% trajectories with 3, 2 and 5 lines.
%
% Tested: MATLAB Version 9.11.0.1769968 (R2021b),
%	      Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% Author: Sven zur Oven-Krockhaus
%	      Institute of Physical and Theoretical Chemistry
%	      University of Tuebingen, Tuebingen, Germany
% E-mail: sven.zur-oven-krockhaus@uni-tuebingen.de
%
% GNU placeholder
%
% Initial release: 2021-12-03
% Last revision: 2023-07-06

%% Function argument validation
arguments
    trackArray (:,1) cell {mustBeTrackArray}
    Options.BoolSelect (:,1) logical = true(numel(trackArray),1)
    Options.colorScheme (1,:) char {mustBeMember(Options.colorScheme, {'binary', ...
       'random4', 'random6', 'random8', 'random10', 'N', 'duration', 'PSF size', 'photons', ...
       'precision', 'displacement', 'range', 'mobility (D)', 'mobility (logD)', ...
       'immobile/mobile (logD)'})} = 'random6'
    Options.binaryMask logical = false(numel(trackArray) ,1)
    Options.MSDparams struct = struct.empty
end

% Check compatibility of chosen options.
checkOptionsCompatibility(trackArray, Options);

%% Main

% Select a set of tracks if that option was specified (otherwise all tracks are selected).
trackArray = trackArray(Options.BoolSelect);

if contains(Options.colorScheme, 'random')
    nCycleCols = sscanf(Options.colorScheme, 'random%u');
    Options.colorScheme = 'random';
end

% Process according to the chosen color scheme.
switch Options.colorScheme
    case 'binary'
        % Assign each track either 1 or 2 depending on the information in |binaryMask|.
        colorValues = ones(numel(trackArray), 1) + Options.binaryMask;

    case 'random'
        % Randomly assign values from 1:|nCycleCols| to each track.
        colorValues = randi(nCycleCols,numel(trackArray),1);

    case 'N'
        % Assign the number of points each track consists of.
        colorValues = cellfun('size', trackArray, 1);

    case 'duration'
        % Get the duration of each trajectory by calculating the range of each time vector (this
        % will include all gaps). Divide by the frame rate to get seconds.
        colorValues = (cellfun(@(x) range(x(:,1)),trackArray) + 1) / Options.MSDparams.frameRate;
    
    case 'PSF size'
        % extract the average PSF size for each track
        colorValues = cellfun(@(x) mean(x(:,4)),trackArray);

    case 'photons'
        % extract the average photon number for each track
        colorValues = cellfun(@(x) mean(x(:,5)),trackArray);

    case 'precision'
        % extract the average CRLB localization precision for each
        % track
        colorValues = cellfun(@(x) mean(x(:,7)),trackArray);

    case 'displacement'
        % Calculate the length of the vector between the first and last coordinate for each
        % track.
        colorValues = cellfun(@(x) norm(diff(x([1,end],2:3))), trackArray);

    case 'range'
        % Get the range of each coordinate vector (x and y) for each track, representing the
        % lengths a and b of an encompassing rectangle. Its hypotenuse is a rough estimation of
        % the particle range.
        colorValues = cellfun(@(x)...
            hypot(range(x(:,2)), range(x(:,3))), trackArray);

    case {'mobility (D)','mobility (logD)','immobile/mobile (logD)'}
        % The apparent diffusion coefficients for each trajectory are calculated.
        MSDresultsArray = calculateMSD(trackArray, Options.MSDparams.frameRate, ...
            'minTrackLength', Options.MSDparams.minTrackLength);
        indivMSDfits = fitMSD(MSDresultsArray, Options.MSDparams.frameRate, ...
            'linFitPoints', Options.MSDparams.linFitPoints);

        % Find valid results and set the invalid tracks to NaN.
        isValid = indivMSDfits(:,1) > 0 & indivMSDfits(:,2) >= Options.MSDparams.adjR2;
        trackArray(~isValid) = {NaN(1,7)};

        % The color values are the individual diffusion coefficients (logarithmized for
        % 'mobility (logD)' and 'immobile/mobile (logD)').
        colorValues = NaN(size(trackArray));
        colorValues(isValid) = indivMSDfits(isValid, 1);
        if ~contains(Options.colorScheme,'(D)')
            colorValues = log10(colorValues);
        end
end


% Append an NaN row to each trajectory, using only the x- and y-coordinates.
tracksNaN = cellfun(@(x) vertcat(x(:,2:3), [NaN NaN]), trackArray, 'uni', 0);

% get the length of each trajectory (including NaN)
tracklengths = cellfun('size', tracksNaN, 1);

% Construct |C| by repeating |colorValues| exactly |tracklengths| times to cover the combined
% trajectory.
C = repelem(colorValues', tracklengths);

% Concatenate and extract coordinates for the function output.
allTracks = vertcat(tracksNaN{:});
X = double(allTracks(:,1)');
Y = double(allTracks(:,2)');
Z = zeros(size(X));

end

%% Subfunctions

function mustBeTrackArray(trackArray)
% argument validation function
if ~all(cellfun("size",trackArray,2) > 2)
    eid = 'Size:wrongDimensions';
    msg = ['Each entry in the track array must have at least three columns ' ...
        '[frame, x-coordinate, y-coordinate].'];
    throwAsCaller(MException(eid,msg))
elseif ~all(cellfun(@(x) all(x(:,1) == floor(x(:,1))), trackArray))
    eid = 'Input:wrongType';
    msg = 'The first column of each track must contain (integer) frame numbers.';
    throwAsCaller(MException(eid,msg))
end

end

function checkOptionsCompatibility(trackArray, Options)
% Checks compatibility of chosen options.

% Check if the chosen |colorScheme| has all required inputs.
switch Options.colorScheme
    case 'binary'
        if ~(isvector(Options.binaryMask) && length(Options.binaryMask) == length(trackArray))
            eid = 'Size:wrongDimensions';
            msg = 'The number of elements in binaryMask and trackArray do not match.';
            throwAsCaller(MException(eid,msg))
        end

    case 'PSF size'
        if ~all(cellfun("size",trackArray,2) > 3)
            eid = 'Size:wrongDimensions';
            msg = 'Entries in the track array are missing the PSF size information.';
            throwAsCaller(MException(eid,msg))
        end

    case 'photons'
        if ~all(cellfun("size",trackArray,2) > 4)
            eid = 'Size:wrongDimensions';
            msg = 'Entries in the track array are missing the photons information.';
            throwAsCaller(MException(eid,msg))
        end

    case 'precision'
        if ~all(cellfun("size",trackArray,2) > 6)
            eid = 'Size:wrongDimensions';
            msg = 'Entries in the track array are missing the CRBL precision information.';
            throwAsCaller(MException(eid,msg))
        end

    case 'duration'
        if isempty(Options.MSDparams)
            eid = 'Missing:parameters';
            msg = 'Structure array containing the frame rate for the duration plot is required.';
            throwAsCaller(MException(eid,msg))

        elseif (~isfield(Options.MSDparams,'frameRate') || ...
                isempty(Options.MSDparams.frameRate))
            eid = 'Missing:parameters';
            msg = 'Structure array must contain the frame rate';
            throwAsCaller(MException(eid,msg))
        end

    case {'mobility (D)','mobility (logD)','immobile/mobile (logD)'}
        if isempty(Options.MSDparams)
            eid = 'Missing:parameters';
            msg = ['Structure array containing additional parameters ' ...
                'for mobility plots is required.'];
            throwAsCaller(MException(eid,msg))

        elseif (~isfield(Options.MSDparams,'frameRate') || ...
                isempty(Options.MSDparams.frameRate)) || ...
                (~isfield(Options.MSDparams,'minTrackLength') || ...
                isempty(Options.MSDparams.minTrackLength)) || ...
                (~isfield(Options.MSDparams,'adjR2') || ...
                isempty(Options.MSDparams.adjR2)) || ...
                (~isfield(Options.MSDparams,'linFitPoints') || ...
                isempty(Options.MSDparams.linFitPoints))
            eid = 'Missing:parameters';
            msg = ['Structure array containing additional parameters ' ...
                'for mobility plots is incomplete.'];
            throwAsCaller(MException(eid,msg))
        end

end

end
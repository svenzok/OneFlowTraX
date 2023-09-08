function [MSDresultsArray, meanMSDresults] = calculateMSD(trackArray, frameRate, Options)
% calculateMSD calculates the MSD of trajectories.
%
% Syntax:
%   [MSDresultsArray, meanMSDresults] = calculateMSD(trackArray, frameRate)
%
% Input Arguments:
%   (Required)
%   trackArray         Particle trajectories as a cell array. Each cell (track) has at least three
%                      columns, in the order [frame, x-coordinate, y-coordinate]. Track coordinates
%                      must be given in nanometers. Tracks with less than three coordinates will be
%                      excluded from the calculation (see output arguments).
%                      (:,1) cell
%
%   frameRate          Camera acquisition frame rate in Hz.
%                      (1,1) double
%   
%   (Optional)
%   minTrackLength     Minimum length of a track (including gaps) to be included in the
%                      calculation of both the individual MSD curves and the mean MSD (default: 1).
%                      (1,1) double
%
% Output Arguments:
%   MSDresultsArray    Cell array, one cell for each track, with columns:
%                      1   time lag in s
%                      2   MSD in µm²
%                      3   sdv in µm²
%                      4   number of datapoints
%                      When a track has insufficient data (less than three coordinates or less
%                      persistent than specified by |minTrackLength|), the corresponding cell will be
%                      filled with a row of zeros to preserve the track index order.
%                      (:,1) cell
%
%   meanMSDresults     Weighted mean results over all valid (see above) MSD curves, with columns:
%                      1   time lag (s)
%                      2   weighted MSD mean for each time lag in µm² 
%                      3   weighted standard deviation in µm²
%                      4   degrees of freedom of the weighted mean
%                      (:,4) double
%
% Other required m-files: none
% Subfunctions: mustBeTrackArray
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.4
%
% Notes:
% Utilizes parts of the program code:
% Maxime Deforet (2013). Kehl, a fast (no loop) method to compute MSD
% (https://www.mathworks.com/matlabcentral/fileexchange/41858-kehl-a-fast-no-loop-method-to-compute-msd),
% MATLAB Central File Exchange. Retrieved December 11, 2021.
%
% The program code was partly rewritten for better integration and faster execution.
%
% This program also calculates the weighted mean as mentioned in the program package msdanalyzer:
% Jean-Yves Tinevez (2014). msdanalyzer (https://github.com/tinevez/msdanalyzer).
% Retrieved November 11, 2022.
%
% Its subfunction getMeanMSD.m was used as a guide to write new code for faster execution.
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
% Initial release: 2021-12-12
% Last revision: 2023-08-25

%% Function argument validation
arguments
    trackArray (:,1) cell {mustBeTrackArray}
    frameRate (1,1) double
    Options.minTrackLength (1,1) double = 1
end

%% Main

% Pre-allocate.
nTracks = numel(trackArray);
MSDresultsArray = cell(nTracks,1);

% Determine which tracks are valid for further processing.
isValid = ...
    (cellfun(@(x) range(x(:,1)), trackArray) + 1 >= Options.minTrackLength) & ...
    (cellfun('size',trackArray,1) >= 3);

% Fill the cells assigned to invalid tracks with rows of zeros and only loop through the valid
% tracks.
MSDresultsArray(~isValid) = {zeros(1,4)};

for iTrack = find(isValid')

    % Extract a single track.
    thisTrack = trackArray{iTrack};

    % Convert the coordinates to micrometers.
    thisTrack(:,2:3) = thisTrack(:,2:3) / 1000;

    % List all time lag combinations (frame distances).
    lagTimeList = internal.stats.pdistmex(thisTrack(:,1)','euc');

    % List all squared point distances.
    sqDistanceList = internal.stats.pdistmex(thisTrack(:,2:3)','euc').^2;

    % Sort the lag times and apply the same sorting order to the squared distances.
    [lagTimeList, sortID] = sort(lagTimeList);
    sqDistanceList = sqDistanceList(sortID);

    % The lag times list now features groups of equal time lags and the respective squared distances
    % in the distances list.
    % For the MSD, these groups have to be averaged. |startID| and |endID| indicate the start and
    % end index of each group. >circshift ~= 0< is a fast method to find the group boundaries.
    startID = find(lagTimeList - circshift(lagTimeList,1) ~= 0);
    endID = find(lagTimeList - circshift(lagTimeList,-1) ~= 0);

    % The number of elements in each group is easily constructed from the group indices.
    nElements = (endID - startID + 1);

    % >cumsum< provides cumulative distances (include 0 for easier indexing later).
    cumDistanceList = cumsum([0,sqDistanceList]);

    % By using the group indices, the sum of each group can be determined in this fast, vectorized
    % way. For the average, the resulting vector just has to be divided by the number of elements
    % calculated before.
    MSD = (cumDistanceList(endID+1) - cumDistanceList(startID))./nElements;

    % The standard deviation can be calculated in a similar fashion. To get (x-mean(x))^2, x has to
    % subtracted by the already calculated average (the MSD) of each group, so the group MSD has to
    % be replicated |nElements| times to be available for each element.
    sqDeviationList = (sqDistanceList - repelem(MSD,nElements)).^2;

    % >cumsum< now provides the cumulative squared deviations.
    cumSquaredDeviationList = cumsum([0, sqDeviationList]);

    % The sum of the squared deviations is determined by the group indices. Then the sdv is
    % calculated by dividing by |nElements| and getting the square root (here, |nElements| is chosen
    % instead of |nElements|-1 to keep the calculation in line with other MSD calculation software).
    sdv = sqrt((cumSquaredDeviationList(endID+1) - cumSquaredDeviationList(startID))./nElements);

    % Store the output. The time lag for each group (first column) can be found by indexing into the
    % first value of |lagTimeList|. Convert to seconds using |frameRate|.
    MSDresultsArray{iTrack} = vertcat(lagTimeList(startID)./frameRate, MSD, sdv, nElements)';
end

% Calculate the mean MSD curve if requested.
if nargout > 1

    % Concatenate all MSD curves.
    allMSD = vertcat(MSDresultsArray{:});

    % Use >unique< and >accumarray< to group equal delays together in a cell array. Each cell will
    % contain (as columns) the MSD and number of data points for each entry in |timeLag|.
    [timeLag, ~, ID] = unique(allMSD(:,1));
    allMSDarray = accumarray(ID, 1:size(allMSD,1), [], @(r){allMSD(r,[2 4])});

    % Average the MSDs in each cell using the number of data points as a weight.
    weightedMeanMSD = cellfun(@(x) sum(prod(x,2))/sum(x(:,2)), allMSDarray, 'uni', 0);

    % Calculating the standard deviation is a bit more complicated (see
    % https://en.wikipedia.org/wiki/Weighted_arithmetic_mean#Reliability_weights).
    % The variance for each element i is calculated by:
    % sum(Ni*(MSDi - weightedMeanMSD)^2) / (sum(Ni) - sum(Ni^2) / sum(Ni)),
    % with Ni (number of datapoints for each element), MSDi (MSD value for each element) and
    % weightedMeanMSD as calculated before. The standard deviation is simply the square root of this
    % calculated variance.
    weightedSdvMSD = cellfun(@(x,y)...
        sqrt(...
        sum(x(:,2).*(x(:,1)-y).^2)./...
        (sum(x(:,2))-sum(x(:,2).^2)/sum(x(:,2)))...
        ), ...
        allMSDarray, weightedMeanMSD, 'uni', 0);

    % The degrees of freedom can be calculated by sum(Ni)^2 / sum(Ni^2). Using weights, this can
    % result in fractional numbers.
    df = cellfun(@(x) sum(x(:,2)).^2./sum(x(:,2).^2), allMSDarray, 'uni', 0);

    % Construct the output (convert the cell arrays to column vectors).
    meanMSDresults = [timeLag, cell2mat(horzcat(weightedMeanMSD, weightedSdvMSD, df))];

    % Cut away the rows with NaNs (the first row and all rows with only one degree of freedom).
    meanMSDresults = meanMSDresults(~any(isnan(meanMSDresults),2),:);

end

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
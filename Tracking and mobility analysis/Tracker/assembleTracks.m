function trackArray = assembleTracks(locList, Options)
% assembleTracks builds tracks from SMLM localization data.
%
% Syntax:
%   trackArray = assembleTracks(locList, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   locList            Array containing the localization data, with at
%                      least three columns in the order:
%                      column 1: the image frame of the localization
%                      column 2: x coordinate in nm
%                      column 3: y coordinate in nm
%                      (:,:) single|double
%
% Name-Value Pair Input Arguments:
%   (Optional)
%   maxGapClose        Maximum gap length (in frames) for which the
%                      algorithm tries to assign a nearby localization to
%                      the same track. As an example, localizations were
%                      connected from frame 4 to 11. Frame 12 has no
%                      localization in the same vicinity. With a value of 3
%                      for |maxGapClose|, the algorithm will look ahead to
%                      frame 15 (3 gap frames: 12, 13, 14). If a nearby
%                      localization is found, it will still be assigned to
%                      the track and the connection process is continued
%                      for the subsequent frames. Otherwise, new
%                      localizations (frame 16 or onwards) will be assigned
%                      to new tracks (default: 3).
%                      (1,1) double
%
%   maxLinkDistance    Maximum distance in nanometers that the algorithm
%                      will search around a localization in the current
%                      frame to connect it to one in the previous frame.
%                      For gaps (see |maxGapClose|), this distance will be
%                      adjusted, assuming free diffusion (default: 150).
%                      (1,1) double
%
%   exportAllInfo      If true, each localization in |trackArray| will
%                      retain all info that was given in |locList|. If
%                      false, only [frame, x coordinate, y coordinate] will
%                      be exported (default: false).
%                      (1,1) logical
%
% Output Arguments:
%   trackArray         Particle trajectories as a cell array. Each cell
%                      (track) has at least three columns, in the order
%                      [frame, x coordinate, y coordinate] and at least two
%                      rows (localizations). For subsequent calculations,
%                      it might be necessary to also delete entries that
%                      have only two localizations. This can be done with
%                      the command:
%                      trackArray(cellfun('size',trackArray,1)<3)=[]
%                      (:,1) cell
%
% Other required m-files: matchpairsSlim
% Subfunctions: checkInputSize
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.3
%
% Notes:
%
% Utilizes concepts and parts of the program code:
%
% Jean-Yves Tinevez (2019). simpletracker - A simple particle tracking
% algorithm for MATLAB that can deal with gaps.
% (https://github.com/tinevez/simpletracker). Retrieved November 17, 2022.
%
% One main difference concerns the linking of localizations, which uses
% MATLAB's >matchpairs< and not the Hungarian algorithm as in Tinevez'
% code. 
% Gap closing does not use the nearest neighbor algorithm here, but again
% uses >matchpairs< to look for connections beyond gap frames. The maximum
% gap lengths are also defined differently - in Tinevez' code, it's the
% difference of frames, here, it's the actual number of gaps (or lost
% frames) between two detected localizations. Thus, 'MaxGapClosing' of 4 in
% >simpletracker< is the same as 'maxGapClose' of 3 in >assembleTracks<.
% Although the Hungarian algorithm is believed to be more exact, results
% are almost identical (~0.1% tracks are merged or separated differently).
% Subsequent calculations of, e.g., diffusion coefficients, show no
% significant deviations.
% However, >assembleTracks< works ~4x faster, tested with ~30000
% localizations over 5000 movie frames, resulting in ~3000 tracks.
%
% Tested: MATLAB Version 9.12.0.1884302 (R2022a),
%	      Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% Author: Sven zur Oven-Krockhaus
%	      Institute of Physical and Theoretical Chemistry
%	      University of Tuebingen, Tuebingen, Germany
% E-mail: sven.zur-oven-krockhaus@uni-tuebingen.de
%
% GNU placeholder
%
% Initial release: 2022-11-16
% Last revision: 2022-11-17

%% Function argument validation
arguments
    locList (:,:) {checkInputSize}
    Options.maxGapClose (1,1) double = 3
    Options.maxLinkDistance (1,1) double = 150
    Options.exportAllInfo (1,1) logical = false
end

%% Rearrange the input data

% Count the occurrences for each frame by using an integer range, add +1 to
% the end as the rightmost bin edge.
rowID = histcounts(locList(:,1), 1:locList(end,1) + 1)';

% Use >mat2cell< to group the rows into cells according to |rowID| (can
% contain 0 as index -> empty cell). Lose the first column, as the cell
% index now serves as the frame number.
locArray = mat2cell(locList(:,2:end), rowID);

% Replace empty entries with an Inf vector that is adjusted to the
% dimensions and class of the input data.
infVec = Inf(1, size(locArray{1}, 2), 'like', locArray{1});
locArray(cellfun('isempty', locArray)) = {infVec};

% Only the point coordinates will be needed here.
pointArray = cellfun(@(x) x(:,1:2),locArray,'uni',0);

%% Link frame-to-frame

% Pre-allocate a cell array, each cell will later contain a list of
% possible point pairings.
pairs = cell(numel(pointArray) - 1, 1);

% Pre-allocate a cell array that will collect all pairs that are found with
% different gap sizes. The first cell contains pairs that were matched
% frame by frame (zero gap), the second cell hold pairs with a frame gap of
% one, and so on. Consequently, this will be a cell array of cell arrays.
allPairs = cell(Options.maxGapClose + 1, 1);

% Compare every frame with the frame directly following it.
source = pointArray(1:end-1);
target = pointArray(2:end);

for p = 1:numel(source)

    if isinf(source{p}(1)) || isinf(target{p}(1))
        % Skip when either the source or target has no localizations.
        continue

    else
        % Calculate euclidean distances, directly use the mex function
        % (skipping the inquiries of >pdist2< for enhanced speed).
        D = internal.stats.pdist2mex(source{p}', target{p}', ...
            'euc', [], [], [], 1);
        
        % Mark all distances greater than the maximal linking distance with
        % Inf to prevent linking.
        D(D > Options.maxLinkDistance) = Inf;
        
        if all(isinf(D), 'all')
            % Skip over combinations that did not result in any viable
            % links.
            continue

        else
            % Link by solving linear assignment problem, using a modified
            % version of the MATLAB function >matchpairs<
            % (>matchpairsSlim<) that skips all unnecessary steps for
            % enhanced speed.
            pairs{p} = matchpairsSlim(D, 1E7);

        end
    end
end

% Fill the results into the |allPairs| array. 
allPairs{1} = pairs;

%% Deal with gaps

% When introducing frame gaps, for example, trying to match localizations
% from frames 5->7, 6->8, 7->9 etc., we have to omit already linked sources
% and targets. This is done by replacing all linked locations in |points|
% with Inf (instead of deleting them, which would mess up the correct
% indexing later on). This must be done separately for the sources and
% targets (for example, a source in frame 10 must still be available as a
% target for frame 8 if its localization is missing in frame 9). In
% contrast to the frame-by-frame connections, the targets are now shifted
% further, according to the currently analyzed frame gap.

if Options.maxGapClose > 0
    
    source = pointArray;
    target = pointArray;
    
    for gap = 1:Options.maxGapClose
        
        for p = 1:numel(pairs)
            if ~isempty(pairs{p})
                source{p}(pairs{p}(:,1),:) = Inf;
                target{p+gap}(pairs{p}(:,2),:) = Inf;
            end
        end
        
        % Pre-allocate a cell array for the pairs.
        pairs = cell(numel(source)-gap-1, 1);
        
        % With frame gaps, the particles have more time to diffuse, thus
        % the maximum linking distance must be adjusted. A linear relation
        % would only make sense for particles moving straight into one
        % direction, and leaving the distance as it is would only apply to
        % immobile particles. The assumption of a random walk is more
        % reasonable, increasing the linking distance by the square root of
        % the frame steps.
        % Edit: taken out for the moment, as this can have a significant
        % impact. Need to simulate tracks with gaps to clear this up.
        maxLinkDistanceGap = Options.maxLinkDistance; %* sqrt(gap + 1);
        
        for p = 1:numel(pairs)

            if all(isinf(source{p}),'all') || all(isinf(target{p+gap+1}),'all')
                % Skip, no localizations left in source or target.
                continue
                
            else
                % Calculate euclidean distances.
                D = internal.stats.pdist2mex(source{p}', target{p+gap+1}', ...
                    'euc', [], [], [], 1);
                
                % Mark all distances greater than the maximal linking
                % distance, also replace NaN with inf.
                D(D > maxLinkDistanceGap | isnan(D)) = Inf;
                
                if all(isinf(D),'all')
                    % Skip over combinations that did not result in any
                    % viable links.
                    continue

                else
                    % Link by solving linear assignment problem.
                    pairs{p} = matchpairsSlim(D, 1E7);

                end
            end
        end
        
        % Add the gap pairs to the |allPairs| array.
        allPairs{gap+1} = pairs;
        
    end
end

% For absolute indices that can later refer to |locList|, the number of
% localizations in each frame are added to the relative pair indices per
% frame.
nLocs = histcounts(locList(:,1), 1:locList(end,1) + 1)';
addN = cumsum(nLocs);
add2target = addN;
add2source = [0; addN(1:end-1)];

for gap = 1:numel(allPairs)
    for p = 1:numel(allPairs{gap})
        if ~isempty(allPairs{gap}{p})
            allPairs{gap}{p} = allPairs{gap}{p} + ...
                [add2source(p), add2target(p+gap-1)];
        end
    end
end

% Join all arrays and unpack.
allPairs = vertcat(allPairs{:});
pairList = cell2mat(allPairs);

%% Graph-based track building

% So far, each localization (source) has only been paired to one other
% localization (target) in another frame. The target can, however, appear
% in another pairing as a source for another target, and so on. The next
% step is to find these chains, i.e. connecting multiple localizations to
% tracks.
% As |pairList| contains the sources in the first column and the targets in
% the second column, the track assignment is easily done using >graph<.
source = pairList(:,1);
target = pairList(:,2);
G = graph(source, target);
refArray = conncomp(G, 'OutputForm', 'cell');

% The cell array |refArray| features one track per cell, but only as
% reference row numbers of |locList|. The following loop replaces the
% references with the corresponding info stored in |locList|.
trackArray = cell(numel(refArray), 1);

if Options.exportAllInfo == false
    % Lose all info except [frame, x coordinate, y coordinate].
    locList = locList(:,1:3);
else
    % Do nothing, export all info as given in |locList|.
end

for iTrack = 1:numel(refArray)
    thisRefs = refArray{iTrack};
    trackArray{iTrack} = locList(thisRefs,:);
end

% Remove entries with a single localization (they are remaining sources
% that did not find a target), which will also remove placeholders with Inf
% values that were previously introduced for frames without any
% localizations.
trackArray(cellfun('size', trackArray, 1) == 1) = [];

end

%% Subfunctions
function checkInputSize(locList)

if size(locList,2) < 3
    eid = 'Size:missingColumns';
    msg = ['The input list of localizations must have at least three ' ...
        'columns in the order: frame, x, y.'];
    throwAsCaller(MException(eid,msg))
elseif size(locList,1) < 2
    eid = 'Size:missingRows';
    msg = ['The input list of localizations must have at least two ' ...
        'localizations (rows)'];
    throwAsCaller(MException(eid,msg))
end

end
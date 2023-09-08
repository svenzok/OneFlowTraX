function [clusterData, boundaryVertexArray] = buildNASTICclusters(trackArray, overlapGroups)
% buildNASTICclusters groups tracks into clusters based on NASTIC.
%
% Syntax:
%   [clusterPointsArray, clusterIndexList, isCore] = buildNASTICclusters(trackArray, overlapGroups)
%
% Input Arguments:
%   (Required)
%   trackArray         Particle trajectories as a cell array. Each cell (track) has at least three
%                      columns, in the order [frame, x-coordinate, y-coordinate] and at least three
%                      rows (localizations).
%                      (:,1) cell
%
%   overlapGroups      Grouped indices as a cell array. Each cell indexes into |trackArray| to
%                      collect all tracks whose bounding boxes were overlapping. Cells with singular
%                      entries contain the non-overlapping remnants.
%                      (:,1) cell
%
% Output Arguments:
%   clusterData        Data calculated for each cluster as a Nx3 list, with columns
%                      1: area in nm² (all input data is assumed to be in nm),
%                      2: diameter in nm,
%                      3: number of points inside the cluster,
%                      4: x position of the cluster's centroid in nm,
%                      5: y position of the cluster's centroid in nm.
%                      (:,5) double
%
%   boundaryVertexArray   Indices for each cluster that provide the sorted border vertices for
%                         its contained tracks.
%                         (:,1) cell
%
% Other required m-files: [mList,pList] = matlab.codetools.requiredFilesAndProducts('fname.m')
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
%
% Independently coded implementation of the procedure described in:
%
% Wallis et al., Molecular Videogaming: Super-Resolved Trajectory-Based Nanoclustering Analysis
% Using Spatio-Temporal Indexing, bioRxiv 2021.09.08.459552,
% doi: https://doi.org/10.1101/2021.09.08.459552.
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
% Initial release: 2023-05-05
% Last revision: 2023-07-26

%% Function argument validation
arguments
    trackArray (:,1) cell {mustBeTracksArray}
    overlapGroups (:,1) cell
end

%% Main

% Pre-allocate the output.
nClusters = numel(overlapGroups);
clusterData = zeros(nClusters, 3);
boundaryVertexArray = cell(nClusters, 1);

% The cluster boundaries should be indexable using the track localizations. From the original
% |trackArray|, all individual locations can extracted with:
allLocs = cell2mat(cellfun(@(x) x(:,2:3), trackArray, 'UniformOutput', false));

% Loop through all clusters to find the localizations of the clustered tracks and also to calculate
% the cluster size, number of enclosed tracks and the cluster centroid position.
for iCluster = 1:numel(overlapGroups)
    thisClusterTracks = trackArray(overlapGroups{iCluster});
    clusterPoints = cell2mat(cellfun(@(x) x(:,2:3), thisClusterTracks, 'UniformOutput', false));

    % Write the number of tracks to the output list.
    clusterData(iCluster, 3) = numel(thisClusterTracks);
    clusterData(iCluster, 4:5) = mean(clusterPoints,1);

    % Create a convex hull around all track localizations in the cluster and calculate the area.
    [k, clusterData(iCluster, 1)] = convhull(double(clusterPoints));

    % Find the indices based on the whole set of points stored in |trackArray|. The output of
    % |convhull| doubles the start and end vertices, so the first one can be ignored.
    [~, LocB] = ismember(clusterPoints(k(2:end),:), allLocs, "rows");
    boundaryVertexArray{iCluster} = LocB';
end

% The diameter is calculated as the diameter of a circle with the same area as the convex hull.
clusterData(:, 2) = sqrt(4*clusterData(:,1)/pi);

end


%% Subfunctions
function mustBeTracksArray(tracks)
% argument validation function
if ~all(cellfun("size",tracks,2) > 2)
    eid = 'Size:wrongDimensions';
    msg = 'Each entry in tracks array must have at least three columns [frame, x-coordinate, y-coordinate].';
    throwAsCaller(MException(eid,msg))
elseif ~all(cellfun("size",tracks,1) >= 2)
    eid = 'Size:wrongDimensions';
    msg = 'Each track must have at least two localizations.';
    throwAsCaller(MException(eid,msg))
end

end
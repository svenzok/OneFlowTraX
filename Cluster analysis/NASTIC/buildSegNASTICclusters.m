function [clusterData, boundaryVertexArray] = buildSegNASTICclusters(trackArray, overlapGroups, overlapThreshold)
% buildSegNASTICclusters groups tracks into clusters based their overlapping segments.
%
% Syntax:
%   [clusterData, boundaryVertexArray] = segNASTIC(trackArray, overlapGroups, overlapThreshold)
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
%   overlapThreshold   Determines how often one track segment has to overlap with segments of
%                      other tracks to be considered as part of a cluster. Overlap here refers to
%                      the segment bounding boxes.
%                      (1,1) double
%
% Output Arguments:
%   clusterData        Data calculated for each cluster as a Nx3 list, with columns
%                      1: area in nm² (all input data is assumed to be in nm),
%                      2: diameter in nm,
%                      3: number of tracks assigned to the cluster (by their centroids).
%                      4: x position of the cluster's centroid in nm,
%                      5: y position of the cluster's centroid in nm.
%                      (:,3) double
%
%   boundaryVertexArray   Indices for each cluster that provide their sorted border vertices. The
%                         points these indices refer to can be generated directly from |trackArray|:
%                         cell2mat(cellfun(@(x) x(:,2:3), trackArray, 'UniformOutput', false));
%                         (:,1) cell
%
% Other required m-files: convertTracksToCentroids, groupAABBoverlaps
% Subfunctions: mustBeTracksArray
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
% This version differs from their procedure in that the regular NASTIC is done first to grab sets of
% clusters that are then fine-tuned with this program. Also, the temporal separation of clusters is
% not implemented. Moreover, 
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
    overlapThreshold (1,1) double
end

%% Main

% Get the track centroids.
trackCentroids = convertTracksToCentroids(trackArray);

% Pre-allocate an array of 2-D polygonal shapes that will hold the boundaries of the subclusters
% (for now, it has the same size as |overlapGroups|, subclusters can later be stored as columns).
subClusters = repmat(polyshape, numel(overlapGroups), 1);

% Loop through each cluster.
for iGroup = 1:numel(overlapGroups)

    % Extract all tracks that are in one cluster.
    thisClusterTracks = trackArray(overlapGroups{iGroup});

    % The number of segments per track are always one less than the contained localizations.
    nSegmentsPerTrack = cellfun('size',thisClusterTracks,1) - 1;

    % Create a block diagonal matrix to mark the segments that belong to the same track.
    blocks = cellfun(@(x) true(x), num2cell(nSegmentsPerTrack), 'UniformOutput', false);
    isFromSameTrack = blkdiag(blocks{:});

    % Pre-allocate and loop through each track.
    nTracks = numel(trackArray);
    segList = cell(nTracks, 1);

    for iTrack = 1:numel(thisClusterTracks)

        % Extract the coordinates of the track localizations. In the following, skipping one or
        % more frames will be ignored (as this is typically a seldom event in relation to the normal
        % frame-to-frame connections).
        thisCoords = thisClusterTracks{iTrack}(:,2:3);

        % Reshape the coordinates to get the segment coordinates (two points for each segment).
        tmp = repelem(thisCoords,2,1);
        tmp = tmp(2:end-1,:)';
        segList{iTrack} = reshape(tmp,2,2,[]);
    end

    segList = cat(3, segList{:});

    % Sort the coordinates to get the bounding boxes.
    AABB = sort(segList, 2);

    % The original publication quantifies the overlap as the number of times a segment from one
    % track overlaps with segments of other tracks. First, get all overlaps with a fast and
    % vectorized method. For a more detailed description of the code, see >groupAABBoverlaps<.
    x0 = squeeze(AABB(1,1,:));
    x1 = squeeze(AABB(1,2,:));
    isOverlappingX = x0' < x1;
    isOverlappingX = isOverlappingX.' & isOverlappingX;
    y0 = squeeze(AABB(2,1,:));
    y1 = squeeze(AABB(2,2,:));
    isOverlappingY = y0' < y1;
    isOverlappingY = isOverlappingY.' & isOverlappingY;
    isOverlappingXY = isOverlappingX & isOverlappingY;

    % Ignore all overlaps within the same track.
    isOverlappingXY = isOverlappingXY & ~isFromSameTrack;

    % The sum of each column shows how often one segment overlaps with another segment from any
    % track but its own.
    nSegOverlaps = sum(isOverlappingXY)';

    % Extract all bounding boxes for the segments that have overlaps above the chosen threshold and
    % group them (for >groupAABBoverlaps< the last row is the frame range which is set to NaN here).
    isValid = find(nSegOverlaps > overlapThreshold);
    validAABB = AABB(:,:, isValid);
    validAABB(3,:,:) = NaN;
    segGroups = groupAABBoverlaps(validAABB, inf);

    % Remove groups that consist of only one segment.
    segGroups(cellfun(@numel, segGroups) <= 1) = [];

    % Create the boundaries for each segment cluster.
    nSegGroups = numel(segGroups);
    for iSub = 1:nSegGroups
        % Extract the clustered segments from the original set.
        clusteredSegs = segList(:,:, isValid(segGroups{iSub}));
        clusteredSegs = reshape(clusteredSegs(:),2,[],1)';

        % Create a convex hull around the points and use >polyshape< to store it for later.
        k = convhull(double(clusteredSegs));
        subClusters(iGroup,iSub) = polyshape(clusteredSegs(k,:));
    end
end

% Rearrange the polyshape array (it might have multiple columns now) to a one column vector and
% delete all polygons without vertices.
subClusters = subClusters(:);
subClusters(arrayfun(@(x) isempty(x.Vertices), subClusters)) = [];

% Pre-allocate the output.
nClusters = numel(subClusters);
clusterData = zeros(nClusters, 3);

% Loop through each of the found clusters to calculate the cluster metrics.
for iCluster = 1:nClusters

    % Get the number of track centroids that are within the convex hull of the found cluster.
    clusterData(iCluster, 3) = sum(isinterior(subClusters(iCluster), trackCentroids));

    % The cluster area is given by the area of the convex hull.
    clusterData(iCluster, 1) = subClusters(iCluster).area;

    % For now, calculate the centroid of the cluster as the centroid of the convex hull (maybe this can be
    % changed to the centroid of the segment points later).
    [clusterCentroidX, clusterCentroidY] = subClusters(iCluster).centroid;
    clusterData(iCluster, 4:5) = [clusterCentroidX, clusterCentroidY];
end

% The diameter is calculated as the diameter of a circle with the same area as the convex hull.
clusterData(:, 2) = sqrt(4*clusterData(:,1)/pi);

% Delete all entries that contain less than two track centroids.
isLacking = clusterData(:,3) < 2;
clusterData(isLacking, :) = [];
subClusters(isLacking, :) = [];

% The cluster boundaries should be indexable using the track localizations. From the original
% |trackArray|, all individual locations can extracted with:
allLocs = cell2mat(cellfun(@(x) x(:,2:3), trackArray, 'UniformOutput', false));

% Pre-allocate the output.
nClusters = numel(subClusters);
boundaryVertexArray = cell(nClusters, 1);

% Loop through all clusters again to find the correct indices.
for iCluster = 1:numel(subClusters)
    thisVertices = subClusters(iCluster).Vertices;
    [~, LocB] = ismember(thisVertices, allLocs, "rows");
    boundaryVertexArray{iCluster} = LocB';
end

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
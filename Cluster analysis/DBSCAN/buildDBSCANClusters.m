function [clusterData, boundaryVertexArray, clusterIndexList, isCore] = buildDBSCANClusters(pointList, epsilon, minPts)
% buildDBSCANClusters groups points into clusters based on DBSCAN.
%
% Syntax:
%   [clusterPointsArray, clusterIndexList, isCore] = buildDBSCANClusters(pointList, epsilon, minPts)
%
% Input Arguments:
%   pointList          List of point coordinates in the format [x-coordinate, y-coordinate].
%                      (:,2) double
%
%   epsilon            Epsilon value to discriminate between clusters and noise.
%                      (1,1) double
%
%   minPts             Minimum number of points that are contained in a circle around one point with
%                      radius epsilon (including the central point), so that the central point is
%                      classified as a core point in a cluster. |minPts| is usually set to two times
%                      the dimensionality of the problem (= 4 for 2-D).
%                      (1,1) double
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
%   boundaryVertexArray   Indices for each cluster that provide the sorted border vertices for its
%                         contained points.
%                         (:,1) cell
%
%   clusterIndexList   Cluster assignments (first output of >dbscan<).
%                      (:,1) double
%
%   isCore             Denotes if a point is a core point (second output of >dbscan<).
%                      (:,1) logical
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.5
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
% Initial release: 2021-12-23
% Last revision: 2023-26-07

%% Function argument validation
arguments
   pointList (:,2) double
   epsilon (1,1) double
   minPts (1,1) double
end

%% Main

% DBSCAN algorithm. 
[clusterIndexList, isCore] = dbscan(pointList, epsilon, minPts);

% Mark the cluster points (everything that is not noise).
isCluster = clusterIndexList ~= -1;

% Group the point indices into cells, one for each cluster index. The found clusters will not be
% further entangled, preserving the original functionality of DBSCAN (for a possibly better cluster
% separation, the input parameters |epsilon| and |minPts| can be changed).
clusterPointsArray = accumarray(clusterIndexList(isCluster), ...
    find(isCluster), [], @(x) {x'});


% The number of points per cluster can already be filled into the output data list.
nClusters = numel(clusterPointsArray);
clusterData = zeros(nClusters,5);
clusterData(:,3) = cellfun(@numel, clusterPointsArray);

% Pre-allocate the boundary index array.
boundaryVertexArray = cell(nClusters, 1);

% There is no general rule how to determine the cluster shapes or sizes found by DBSCAN. For tracing
% the cluster borders, >boundary< usually delivers robust results. Unlike >convhull<, >boundary< can
% shrink towards the cluster interior to envelop the points (based on >alphaShape<). The default
% value for this function is 0.5.
for iCluster = 1:nClusters
    thisClusterPointIDs = clusterPointsArray{iCluster};
    thisClusterPoints = pointList(thisClusterPointIDs,:);

    [k, area] = boundary(thisClusterPoints);
    clusterData(iCluster, 1) = area;
    clusterData(iCluster, 4:5) = mean(thisClusterPoints,1);
    boundaryVertexArray{iCluster} = thisClusterPointIDs(k(2:end));

    % The diameter of the cluster can be estimated using several methods. For example, SR-Tesseler
    % applies a PCA to the underlying point cloud, while ClusterViSu/SharpViSu uses the diameter of
    % a circle that has the same area as the cluster. Here, we choose the SR-Tesseler method,
    % fitting a two-component normal distribution to the data to get the average FWHM as the
    % diameter. If needed, the diameter according to the ClusterViSu/SharpViSu method can be quickly
    % calculated from the area.

    if size(thisClusterPoints, 1) <= 3
        % A two-component distribution fit needs at least three data points. Just leave the
        % list entry empty when there are too few points in the cluster.
    else
        GMModel = fitgmdist(thisClusterPoints, 1);
        clusterData(iCluster, 2) = 2*sqrt(2*log(2)) * sqrt(trace(GMModel.Sigma)/2);
    end
end

end
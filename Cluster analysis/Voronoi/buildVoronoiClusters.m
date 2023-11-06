function [clusterData, boundaryVertexArray] = buildVoronoiClusters(DT, voronoiVertices, voronoiCells, method, valueList, splitValue)
% buildVoronoiClusters groups points into clusters according to a chosen threshold value (based on
%                      the Voronoi tessellation methods of SR-Tesseler and SharpViSu/ClusterViSu).
%
% Syntax:
%   [clusterData, boundaryVertexArray] = ...
%             buildVoronoiClusters(DT, voronoiVertices, voronoiCells, method, valueList, splitValue)
%
% Input Arguments:
%   DT                 Delaunay triangulation object, as created with >delaunayTriangulation<.
%                      (:,3) delaunayTriangulation
%
%   voronoiVertices    2-D Voronoi vertices as they would be produced by the first output of the
%                      functions >voronoin< or >voronoiDiagram<.
%                      (:,2) double
%
%   voronoiCells       2-D Voronoi cells as they would be produced by the second output of the
%                      functions >voronoin< or >voronoiDiagram<.
%                      (:,1) cell
%
%   method             Method to discriminate between clusters and noise, using either the Voronoi
%                      cell areas or their normalized first-rank densities.
%                      (1,:) char {'area', 'density'}
%
%   valueList          List of either Voronoi cell areas or normalized first-rank densities in the
%                      same order as the points they were calculated from.
%                      (:,1) double
%
%   splitValue         Threshold value to discriminate between clusters and noise.
%                      (1,1) double
%
% Output Arguments:
%   clusterData        Data calculated for each cluster as a Nx3 list, with columns
%                      1: area in nm² (all input data is assumed to be in nm),
%                      2: diameter in nm,
%                      3: number of points inside the cluster.
%                      4: x position of the cluster’s centroid in nm
%                      5: y position of the cluster's centroid in nm
%                      The calculation of the area and diameter are different depending on the
%                      chosen method.
%                      (:,5) double
%
%   boundaryVertexArray   Indices for each cluster that provide the sorted border vertices for
%                         either its Voronoi cells ('area' method) or its contained points
%                         ('density' method).
%                         (:,1) cell
%
% Other required m-files: calculateVoronoiArea
% Subfunctions: traceVoronoiClusters
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
% Initial release: 2021-12-06
% Last revision: 2023-07-26

%% Function argument validation
arguments
   DT (:,3) delaunayTriangulation
   voronoiVertices (:,2) double
   voronoiCells (:,1) cell
   method (1,:) char {mustBeMember(method,["area","density"])}
   valueList (:,1) double
   splitValue (1,1) double
end

%% Main

% Tag the Voronoi cells that comply with the chosen threshold value.
if strcmp(method, 'area')
    % Use the Voronoi cell area for segmentation.
    validCellIDs = valueList < splitValue;
elseif strcmp(method, 'density')
    % Use the normalized first-rank density for segmentation.
    validCellIDs = valueList > splitValue;
else
    % empty
end

% For each valid point, determine the distance to the nearest neighbor and calculate the median of
% all found distances. This value will be used later to possibly untangle merged clusters.
% XXXX maybe take the mean of the distances in each cluster.
[~,D] = knnsearch(DT.Points(validCellIDs,:), DT.Points(validCellIDs,:), 'K', 2);

if sum(D(:)) == 0 || all(isempty(D(:)))
    % If only one valid Voronoi cell is found, or if none are found at all, exit the program and return one cluster with
    % only one dummy point (it will be filtered out later).
    clusterData = zeros(1,5);
    boundaryVertexArray = {1};
    return
else
    medianPointDistance = median(D(:,2));
end

% Get an array of boundary vertices (tracing along the border edges of the corresponding Voronoi
% cells), using a subfunction.
boundaryVertexArray = traceVoronoiClusters(voronoiCells(validCellIDs));

% Collect the corresponding cluster points.
nClusters = numel(boundaryVertexArray);
clusterPointsArray = cell(nClusters,1);

% When the Voronoi cells were cut at the image boundary, some edges will be on the same line,
% for which >polyshape< will issue a warning that can be safely silenced here.
warning('off','MATLAB:polyshape:repairedBySimplify');

for iCluster = 1:nClusters
    boundaryVertexIDs = boundaryVertexArray{iCluster};

    % Create a polyshape object from the Voronoi cluster boundary vertices.
    pgon = polyshape(voronoiVertices(boundaryVertexIDs,:));

    % Collect the point IDs that are enclosed by the boundary.
    clusterPointIDs = find(isinterior(pgon, DT.Points));
    clusterPointsArray{iCluster} = clusterPointIDs;
end

% Re-enable the warning.
warning('on','MATLAB:polyshape:repairedBySimplify');

% Try to untangle each cluster, using >alphaShape< with default parameters (so no points will be
% lost) to see if the cluster can be separated into several regions. Described simply, this
% algorithm tries to push a ball with radius alpha between the points to "scoop out" polygonal
% regions that tightly enclose the points. The default value for alpha chooses the smallest
% value that does not leave single points behind.

for iCluster = 1:nClusters
    clusterPointIDs = clusterPointsArray{iCluster};
    S = alphaShape(DT.Points(clusterPointIDs,:));

    % Compare the calculated default alpha value for |S| with the previously calculated median
    % distance and choose the higher one (to better conserve the underlying data, double the median
    % distance).
    S.Alpha = max(S.Alpha, 2*medianPointDistance);

    % When this still results in multiple regions, separate the points into subclusters.
    nRegions = S.numRegions;
    subPointIDs = cell(nRegions,1);
    if nRegions ~=1
        for iRegion = 1:nRegions
            % Get the point IDs from each region.
            isInside = inShape(S, DT.Points(clusterPointIDs,:), iRegion);
            subPointIDs{iRegion} = clusterPointIDs(isInside);
        end

        % The >alphaShape< algorithm might sometimes put points into new regions even if the Voronoi
        % cells are tightly packed (but the points inside are distributed in a way that might
        % suggest an empty region). The following is based on a conservative approach that
        % only separates the clusters for very clear cases.

        if nRegions == 2 && any(cellfun(@numel,subPointIDs) < 5)
            % Only two regions are found and one of them (or both) contain less than five points.
            % This is an accidental separation and will be ignored (the original cluster is
            % preserved).
            subPointIDs = {clusterPointIDs};

        elseif nRegions > 2 && sum(cellfun(@numel,subPointIDs) < 5) == 1
            % Multiple larger regions are found, but one contains less than five points.
            disp("This happens") %XXXXX

        else
            % Do nothing, all found regions are valid and their separated points will be processed
            % accordingly in the following.
        end

        % Trace the new boundary around each subcluster.
        for iSub = 1:numel(subPointIDs)
            % Overwrite the original arrays holding the cluster points and their outer Voronoi
            % borders, using additional columns.
            pointIDs = subPointIDs{iSub};
            thisVoronoiCells = voronoiCells(pointIDs);
            subBoundaryVertexArray = traceVoronoiClusters(thisVoronoiCells);
            clusterPointsArray{iCluster, iSub} = pointIDs;
            boundaryVertexArray(iCluster, iSub) = subBoundaryVertexArray;
        end

        % Rearrange the arrays (they might have multiple columns now) to one-column arrays.
        clusterPointsArray = clusterPointsArray(:);
        clusterPointsArray(cellfun(@isempty,clusterPointsArray)) = [];
        boundaryVertexArray = boundaryVertexArray(:);
        boundaryVertexArray(cellfun(@isempty,boundaryVertexArray)) = [];
    end
end

% The number of points per cluster and the centroid positions of the clusters can already be filled into the output data
% list.
nClusters = numel(clusterPointsArray);
clusterData = zeros(nClusters,5);
clusterData(:,3) = cellfun(@numel, clusterPointsArray);
clusterData(:,4:5) = cell2mat(cellfun(@(x) mean(DT.Points(x,:),1), clusterPointsArray, 'UniformOutput', false));

% The cluster sizes and diameters must be calculated for the chosen method.
switch method
    case 'area'
        % The area according to the method used in ClusterViSu/SharpViSu is given by the area of the
        % clustered Voronoi cells. We can use >calculateVoronoiArea< to calculate the area,
        % putting in |boundaryVertexArray| as the Voronoi cells argument (each entry traces around a
        % cluster instead of a single Voronoi cell).
        clusterData(:,1) = calculateVoronoiArea(voronoiVertices, boundaryVertexArray);

        % The diameter is defined here as the diameter of a circle with the same area as the
        % corresponding cluster area.
        clusterData(:,2) = sqrt(4*clusterData(:,1)/pi);

        % For visualization purposes, the output boundary array is |boundaryVertexArray|
        % that can index into the Voronoi vertices.

    case 'density'
        % SR-Tesseler estimates the cluster area by drawing a polygon around the cluster points. The
        % publication does not describe exactly how this is done. Creating a convex hull or an alpha
        % shape would be fastest, but to be most accurate, we can use |boundaryVertexArray| that
        % contains the vertices of the border Voronoi cells in sorted order to go around the
        % contained points step by step.

        % Extract all edges of the thresholded Voronoi cells, connecting the last vertex of each
        % cell back to the first. This both reduces the number edges that need to be compared and
        % simplifies the assigment (any cluster border edge can only match one Voronoi cell).
        allEdgesArray = cellfun(@(x) [x',circshift(x,-1)'], voronoiCells(validCellIDs), ...
            "UniformOutput", false);
        allEdges = cell2mat(allEdgesArray);

        % Pre-allocate.
        pointBoundaryVertexArray = cell(size(boundaryVertexArray));

        % XXXXX When clusters are separated, edges and their assignment to only one point sometimes
        % fail.

        for iCluster = 1:nClusters

            % Use the sorted boundary vertices to create sorted edges. Sort the columns so that the
            % lower index comes first.
            boundaryVertexIDs = boundaryVertexArray{iCluster};
            sortedOuterEdges = [boundaryVertexIDs', circshift(boundaryVertexIDs,-1)'];
            sortedOuterEdges = sort(sortedOuterEdges,2);

            % Find which edge of |allEdges| is the respective edge on the Voronoi border.
            sortedAllEdges = sort(allEdges,2);
            [~, edgeIDs] = ismember(sortedOuterEdges, sortedAllEdges, 'rows');

            % Create an index vector for |allEdges| to track back which edges belong to which
            % original Voronoi cell.
            nEdgesPerCell = cellfun(@numel, voronoiCells(validCellIDs));
            edgeIDList = repelem(find(validCellIDs), nEdgesPerCell);

            % The IDs of the Voronoi cells can now be extracted - in the order of the edge that runs
            % around the cluster. Use >unique< while preserving the order as one Voronoi cell might
            % have several border edges.
            pointBoundaryVertexArray{iCluster} = unique(edgeIDList(edgeIDs), 'stable')';
        end

        % Calculate the area, using |pointBoundaryVertexArray| as the indexing input for the points
        % stored in |DT.Points|.
        clusterData(:,1) = calculateVoronoiArea(DT.Points, pointBoundaryVertexArray);

        % According to SR-Tessleler, the length and width of each cluster is analyzed by PCA,
        % determining sigma along the principal components and converting them to equate the FWHM of
        % a Gaussian distribution. The diameter then is defined as the average of length and width.
        % Here, in a similar fashion, the points are fitted with a bivariate normal distribution,
        % calculating the FWHM for each component and returning their average.
        for iCluster = 1:nClusters
            thisClusterPoints = DT.Points(clusterPointsArray{iCluster},:);
            if size(thisClusterPoints, 1) <= 3
                % A two-component distribution fit needs at least three data points. Just leave the
                % list entry empty when there are too few points in the cluster.
            else
                GMModel = fitgmdist(thisClusterPoints, 1);
                clusterData(iCluster, 2) = 2*sqrt(2*log(2)) * sqrt(trace(GMModel.Sigma)/2);
            end
        end

        % For visualization purposes, the output boundary array will be |pointBoundaryVertexArray|
        % that can index into the point coordinates (e.g. |DT.Points|).
        boundaryVertexArray = pointBoundaryVertexArray;
end


end

% Subfunctions
function boundaryVertexArray = traceVoronoiClusters(voronoiCells)

% Extract all edges of the input Voronoi cells, connecting the last vertex of each cell back to the
% first.
allEdgesArray = cellfun(@(x) [x',circshift(x,-1)'], voronoiCells, "UniformOutput", false);
allEdges = cell2mat(allEdgesArray);

% Voronoi cells of the same cluster will be connected by their edges. They can be grouped quickly
% using >graph<.
G = graph(allEdges(:,1), allEdges(:,2));

% Extract the clusters by separating the connected graph components. This will also produce
% components with single points (which is normal) that can be deleted.
clusterCellArray = conncomp(G, "OutputForm", "cell")';
clusterCellArray(cellfun(@numel, clusterCellArray) == 1) = [];

% Loop through each cluster to find the outer Voronoi cell borders.
nClusters = numel(clusterCellArray);
boundaryVertexArray = cell(nClusters,1);

for iCluster = 1:nClusters
    thisClusterVertexIDs = clusterCellArray{iCluster};

    % Extract all edges belonging to the current cluster.
    isClusterEdge = ismember(allEdges(:,1), thisClusterVertexIDs) | ...
        ismember(allEdges(:,2), thisClusterVertexIDs);
    clusterEdges = allEdges(isClusterEdge,:);

    % Sort each row to make sure all segments have their vertices in the same order and then
    % sort all rows in ascending order.
    sortedClusterEdges = sortrows(sort(clusterEdges,2));
       
    % All interior edges are shared by two Voronoi cells, so the border edges are the only ones that
    % are unique.
    k = find([true; any(diff(sortedClusterEdges, 1, 1), 2); true]);
    outerEdges  = sortedClusterEdges(k(diff(k) == 1), :);

    % Use >graph< and >allcycles< to quickly connect these edges and put them into the correct
    % order.
    C = graph(outerEdges(:,1), outerEdges(:,2));
    cycles = allcycles(C);
    if numel(cycles) > 1
        % If multiple cycles exist, the smaller ones represent interior edges that form holes. These
        % will be ignored as the actual boundary is given by the largest cycle.
        [~, boundaryCycleID] = max(cellfun(@numel, cycles));
        cycles = cycles(boundaryCycleID);
    end
    boundaryVertexArray{iCluster} = cycles{:};
end

end
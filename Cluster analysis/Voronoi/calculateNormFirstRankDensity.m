function normFirstRankDensityList = calculateNormFirstRankDensity(DT, voronoiAreaList)
% calculateNormFirstRankDensity calculates the normalized first-rank density for Voronoi cells.
%
% Syntax:
%   normFirstRankDensityList = calculateNormFirstRankDensity(DT, voronoiAreaList)
%
% Input Arguments:
%   (Required)
%   DT                 Delaunay triangulation object, as created with >delaunayTriangulation<.
%                      (:,3) delaunayTriangulation
%
%   voronoiAreaList    Area list of the sample data.
%                      (:,1) double
%
% Output Arguments:
%   normFirstRankDensityList   Normalized first-rank density list of the sample data.
%                              (:,1) double
%
% Other required m-files: none
% Subfunctions: mustHaveMatchingPointNumbers
% Additional required MATLAB products: none
%
%
% Coded along the procedure described in the publication:
%
% Levet, F. et al.: SR-Tesseler: a method to segment and quantify localization-based
% super-resolution microscopy data. Nat Methods 12, 1065–1071 (2015).
% https://doi.org/10.1038/nmeth.3579
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
% Last revision: 2023-04-03

%% Function argument validation
arguments
    DT (:,3) delaunayTriangulation
    voronoiAreaList (:,1) double {mustHaveMatchingPointNumbers(DT,voronoiAreaList)}
end

%% Main

% Find the Delaunay triangles that are connected to each vertex.
indexDT = vertexAttachments(DT);

% Get their vertex IDs (= coordinate IDs = Voronoi cell IDs) and use >unique< to collect shared
% vertices only once. This will include the original "seed" vertex.
vertexIDs = cellfun(@(x) unique(DT.ConnectivityList(x,:)), indexDT, "UniformOutput", false);

% Calculate the first-rank density for each seed:
% The Voronoi cell of point P directly touches other Voronoi cells around it (these are those of the
% first rank). Take the number of all Voronoi cells of the first rank (including the one of point P)
% and divide by their collective area.
firstRankDensityList = cellfun(@(x) numel(x)/sum(voronoiAreaList(x)), vertexIDs);

% Divide by the average density for the normalized first-rank density.
meanDensity = 1/mean(voronoiAreaList, 'omitnan');
normFirstRankDensityList = firstRankDensityList / meanDensity;

end

%% Subfunctions
function mustHaveMatchingPointNumbers(DT,valueList)
% argument validation function
if ~isequal(length(DT.Points),length(valueList))
    eid = 'Length:notEqual';
    msg = 'Number of points in DT does not match the valueList.';
    throwAsCaller(MException(eid,msg))
end

end
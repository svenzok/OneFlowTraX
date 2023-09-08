function [DT, voronoiVertices, voronoiCells] = getTrimmedVoronoi(pointList, Options)
% getTrimmedVoronoi returns the Voronoi vertices and regions, trimmed at a boundary.
%
% Syntax:
%   [voronoiVertices, voronoiCells] = getTrimmedVoronoi(pointList, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   pointList          List of point coordinates in the format [x-coordinate, y-coordinate].
%                      (:,2) double
%
%   (Optional)
%   type               Specifies the trim type. All options will remove all Voronoi cells with at
%                      least one vertex placed at [Inf, Inf]. The corresponding row in |DT.Points|
%                      will also be deleted to keep the correct order. The default option
%                      'removeInf' will do that and nothing else. Alternatively, 'remove' will also
%                      remove all Voronoi cells that have at least one vertex outside a specified
%                      boundary, while 'intersect' intersects affected Voronoi cells with the
%                      boundary and replaces them with their trimmed versions
%                      (default: 'removeInf').
%                      (1,:) char {'removeInf', 'remove', 'intersect'}
%
%   boundary           Specifies the boundary for the 'type' options 'remove' and 'intersect'. It
%                      can either be specified as 'convhull' to use the convex hull of |DT| (i.e., a
%                      convex hull around the points stored in |DT|), as 'AABB' to use an axis-
%                      aligned bounding box around those points, or as a 1x4 vector specifying
%                      the boundaries in the format [xmin, xmax, ymin, ymax] (default: 'AABB').
%                      (1,:) char {'convhull', 'AABB'} | (1,4) double
%
% Output Arguments:
%   DT                 Delaunay triangulation object, with all invalid points removed.
%                      (:,3) delaunayTriangulation
%
%   voronoiVertices    2-D Voronoi vertices as they would be produced by the first output of the
%                      functions >voronoin< or >voronoiDiagram<. Depending on the chosen options,
%                      vertices are added from the intersection of Voronoi cells with the boundary.
%                      (:,2) double
%
%   voronoiCells       2-D Voronoi cells as they would be produced by the second output of the
%                      functions >voronoin< or >voronoiDiagram<. Voronoi cells that would extend
%                      into infinity are removed. Depending on the chosen options, Voronoi cells
%                      that intersect with the boundary have updated indices that now also index
%                      into the added vertices of |voronoiVertices|.
%                      (:,1) cell
%
% Other required m-files: none
% Subfunctions: checkBoundaryOption
% Additional required MATLAB products: none
%
% Notes:
% 
% This program is used to adjust the perimeter regions when constructing Voronoi diagrams. Instead
% of just disregarding all Voronoi cells that extend beyond the borders (e.g., extents of the point
% cloud, coordinate limits of the image the points were calculated for), the user has also the
% option to intersect affected Voronoi cells with the border, "rescuing" some of those cells for
% calculations in certain situations (they would otherwise be unusable because of their
% disproportionately large areas).
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
% Initial release: 2021-12-20
% Last revision: 2023-04-06

%% Function argument validation
arguments
    pointList (:,2) double
    Options.type {mustBeMember(Options.type,["removeInf","remove","intersect"])} = 'removeInf'
    Options.boundary {checkBoundaryOption} = 'AABB'
end

%% Main

% Create the Voronoi vertices and cells.
DT = delaunayTriangulation(pointList);
[voronoiVertices, voronoiCells] = voronoiDiagram(DT);

% Remove all Voronoi cells when they reference to [Inf, Inf], which is the standard first entry in
% |voronoiVertices|. To keep the correct references to |DT.Points| and |DT.ConnectivityList|, the
% cell contents are first replaced with 1 and removed later.
hasInf = cellfun(@(x) any(x==1, 2), voronoiCells);
voronoiCells(hasInf) = {1};

% When the 'type' option is 'removeInf', remove all invalid Voronoi cells and the corresponding
% points from |DT| and let the function end here.
if strcmp(Options.type, 'removeInf')
    DT.Points(find(hasInf),:) = [];
    voronoiCells(hasInf) = [];
    return
end

% Create the boundary as a polyshape object. Use either a convex hull, an axis-aligned bounding box
% around all points or a user-specified rectangle.
if isstring(Options.boundary) || ischar(Options.boundary)
    Options.boundary = string(Options.boundary);
    if Options.boundary == "convhull"
        C = convexHull(DT);
        boundaryPolygon = polyshape(DT.Points(C(2:end),:));
    elseif Options.boundary == "AABB"
        [xmin, xmax] = bounds(DT.Points(:,1));
        [ymin, ymax] = bounds(DT.Points(:,2));
        boundaryPolygon = polyshape([xmin xmax, xmax, xmin], [ymin, ymin, ymax, ymax]);
    end
else
    boundaryPolygon = polyshape(Options.boundary([1 3; 2 3; 2 4; 1 4]));
end 

% Get all Voronoi cell vertices that lie outside the boundary.
outIDs = find(~isinterior(boundaryPolygon, voronoiVertices));

% Get the IDs of affected Voronoi cells that contain at least one vertex ID that is in |outIDs|,
% but exclude those Voronoi cells that had vertices at [Inf, Inf] as they were already set to 1.
outCellIDs = find(cellfun(@(x) any(ismember(x,outIDs)), voronoiCells) & ~hasInf);

if strcmp(Options.type, 'remove')
    % Replace the affected Voronoi cells with 1.
    voronoiCells(outCellIDs) = {1};

elseif strcmp(Options.type, 'intersect')
    % Pre-allocate a cell array that will hold the new vertices for the affected cells.
    newVertices = cell(size(outCellIDs));

    % Cut the affected cells at the boundary.
    for iCell = outCellIDs'
        % Extract coordinates and remove Inf values.
        thisCell = voronoiVertices(voronoiCells{iCell},:);
        thisCell(any(isinf(thisCell),2),:) = [];

        % Get the intersection using polyshapes.
        newPolyshape = intersect(boundaryPolygon, polyshape(thisCell));

        % Extract the new coordinates and store them.
        newVertices{iCell} = newPolyshape.Vertices;
    end

    % Collect all new coordinates and remove duplicates that are already in the original vertex list.
    allNewVertices = vertcat(newVertices{:});
    isDuplicate = ismembertol(allNewVertices, voronoiVertices, 'ByRows', true);
    allNewVertices(isDuplicate,:) = [];

    % Append the new vertices to the original list.
    voronoiVertices = [voronoiVertices; allNewVertices];

    % For each of the reshaped voronoi cells, assign the vertex coordinates to the corresponding
    % indices.
    for iCell = outCellIDs'
        [~, idx] = ismembertol(newVertices{iCell}, voronoiVertices, 'ByRows', true);
        voronoiCells{iCell} = idx';
    end

end

% Remove all invalid Voronoi cells and the corresponding points from |DT|.
invalidIDs = find(cellfun(@(x) all(x==1,2), voronoiCells));
voronoiCells(invalidIDs) = [];
DT.Points(invalidIDs, :) = [];


end

%% Subfunctions
function checkBoundaryOption(boundary)
% argument validation function
if isstring(boundary) || ischar(boundary)
    if string(boundary) == "convhull" || string(boundary) == "AABB"
        % These are the expected string options, everything is fine.
    else
        eid = 'boundary:unknownInput';
        msg = 'The boundary must either be specified as ''convhull'', ''AABB'' or as a 1x4 numerical vector';
        throwAsCaller(MException(eid,msg))
    end
elseif isnumeric(boundary) && length(boundary) == 4
    % This is the expected numerical input, everything is fine.
else
    eid = 'boundary:unknownInput';
    msg = 'The boundary must either be specified as ''convhull'', ''AABB'' or as a 1x4 numerical vector';
    throwAsCaller(MException(eid,msg))
end

end
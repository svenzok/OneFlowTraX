function [patchFaces, patchEdges] = preparePatchData(vertexList, connectivityArray)
% preparePatchData converts connected and grouped vertices to patch object inputs for plotting.
%
% Syntax:
%   [patchFaces, patchEdges] = preparePatchData(vertexList, connectivityArray)
%
% Input Arguments:
%   (Required)
%   vertexList         List of all vertex coordinates that will be indexed into with
%                      |connectivityArray|.
%                      (:,2) double
%
%   connectivityArray   Each cell contains a list of indices showing which vertices in
%                       |vertexList| must be connected in which order to create one polygon. The
%                       last vertex in each cell will be connected the the first vertex in each cell
%                       to create closed polygons.
%                       (:,1) cell
%
% Output Arguments:
%   patchFaces         The Faces property values for a patch object that can be used together with
%                      |vertexList| to plot all polygons in |connectivityArray| at once.
%                      Example syntax:
%                      patch('Faces', patchFaces, 'Vertices', vertexList)
%                      (:,:) double
%
%   patchEdges         XData and YData property values for line object, one column for each. This
%                      output can be used to color the patches in a different style (transparency,
%                      interpolated color etc.), which would otherwise be the same for both the
%                      lines and faces of the patch object.
%                      Example syntax (the patch object's 'EdgeColor' must be set to 'none'):
%                      line(patchEdges(:,1), patchEdges(:,2))
%                      (:,2) double
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% The input array |connectivityArray| holds the information how the vertices in |vertexList|
% must be connected to form faces. >patch< can plot multiple faces at once when the Faces property
% value is a matrix (one row per face). As not all polygons will have the same number of vertices,
% padding with NaN is allowed here. The edges can be plotted as a separate line object to allow
% for a different colors or transparency. Putting all lines into one object, separated by NaN, is
% the fastest option for plotting.
% The polygons of a patch object can be colored individually by setting 'FaceColor' to 'flat' and
% 'FaceVertexCData' to values (e.g., their area) that are interpreted by a colormap.
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
% Last revision: 2023-03-28

%% Function argument validation
arguments
    vertexList (:,2) double
    connectivityArray (:,1) cell
end

%% Main

% Determine the maximum number of vertices used in |connectivityArray|.
maxVertices = max(cellfun('size', connectivityArray, 2));

% Pad the vertex list in each cell with NaN and convert to a matrix.
patchFaces = cellfun(@(x) horzcat(x, NaN(1, maxVertices - numel(x))), ...
    connectivityArray, "UniformOutput", false);
patchFaces = cell2mat(patchFaces);

% For the edges, duplicate the first vertex to ensure that the lines completely enclose each patch 
% and add NaN as a separator.
patchEdges = cellfun(@(x) ...
    [vertexList(x,:); ...
     vertexList(x(1),:); ...
     NaN(1,2)], ...
    connectivityArray, "UniformOutput", false);
patchEdges = cell2mat(patchEdges);

end
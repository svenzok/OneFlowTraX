function voronoiAreaList = calculateVoronoiArea(voronoiVertices, voronoiCells)
%calculateVoronoiArea computes the area of Voronoi cells. 
%
% Syntax:
%   voronoiAreaList = calculateVoronoiArea(voronoiVertices, voronoiCells)
%
% Input Arguments:
%   (Required)
%   voronoiVertices    2-D Voronoi vertices as produced by the first output of the functions
%                      >voronoin< or >voronoiDiagram<.
%                      (:,2) double
%
%   voronoiCells       2-D Voronoi cells as produced by the second output of the functions
%                      >voronoin< or >voronoiDiagram<.
%                      (:,1) cell
%
% Output Arguments:
%   voronoiAreaList    Area vector in the same order as |voronoiCells|
%                      (:,1) double
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Low-level version of >polyarea<, but uses >circshift< instead of determining the vertex number,
% making it faster than >polyarea< by a factor of ~ 3.5
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
% Initial release: 2021-03-12
% Last revision: 2022-03-06

%% Function argument validation
arguments
    voronoiVertices (:,2) double
    voronoiCells (:,1) cell
end

%% Main

voronoiAreaList = zeros(numel(voronoiCells),1);
for iCell = 1:numel(voronoiCells)
    x = voronoiVertices(voronoiCells{iCell},1);
    y = voronoiVertices(voronoiCells{iCell},2);
    voronoiAreaList(iCell) = ...
        abs(sum((y + circshift(y,1)) / 2 .* (x - circshift(x,1))));
end

end
function [distanceList, intersectionPointList] = calculatePointsToLineDistance2D(pointList, lineVector)
% calculatePointsToLineDistance2D determines the shortest distance from each point to a line in two dimensions. 
%
% Syntax:
%   [distanceList, intersectionPointList] = ...
% calculatePointsToLineDistance2D(pointsX, pointsY, lineX, lineY)
%
% Input Arguments:
%   (Required)
%   pointList          List of point coordinates in the format [x-coordinate, y-coordinate].
%                      (:,2) double
%
%   lineVector         Line vector in the format [x1 x2; y1 y2].
%                      (2,2) double
%
% Output Arguments:
%   distanceList       Distances from each point to the line.
%                      (:,1) double
%
%   intersectionPointList   Projection coordinates of the points onto the line (intersection of the
%                           line tangent passing through each respective point).
%                           (:,2) double
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Coded along coded along the procedure described in:
% Bourke, P. (1988, October). Minimum Distance between a Point and a Line.
% Points, lines, and planes. http://paulbourke.net/geometry/pointlineplane/
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
% Initial release: 2021-12-22
% Last revision: 2023-03-07

%% Function argument validation
arguments
    pointList (:,2) double
    lineVector (2,2) double
end

%% Main

% The line is defined by points P1 (x1,y1) and P2 (x2,y2), and the query
% points are given by P3 (x3, y3).
x1 = lineVector(1,1);
x2 = lineVector(1,2);
y1 = lineVector(2,1);
y2 = lineVector(2,2);
x3 = pointList(:,1);
y3 = pointList(:,2);

% The line can also be expressed as the equation:
% P = P1 + u * (P2-P1).
% As point P3 (x3, y3) is closest to the line at the tangent to the line which passes through P3,
% the dot product of the tangent and line is zero:
% (P3 - P) dot (P2 - P1) = 0.
% Substituting P with the line equation and solving for u gives:
u = ((x3 - x1)*(x2 - x1) + (y3 - y1)*(y2 - y1)) / ...
    norm([x1 x2; y1 y2])^2;
% (MATLAB's implicit expansion feature automatically replicates the line coordinates for each point
% in P3.)

% The intersection points (x, y) for each tangent are given by:
x = x1 + u*(x2 - x1);
y = y1 + u*(y2 - y1);
intersectionPointList = [x, y];

% The corresponding distances then are calculated with:
distanceList = vecnorm([x, y]' - [x3, y3]')';

end
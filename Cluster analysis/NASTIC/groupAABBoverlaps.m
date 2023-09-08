function overlapGroups = groupAABBoverlaps(AABB, timeWindow)
%groupAABBoverlaps finds overlapping axis-aligned bounding boxes, grouping their indices.
%
% Syntax:
%   overlapGroups = groupAABBoverlaps(AABB)
%   overlapGroups = groupAABBoverlaps(AABB, timeWindow)
%
% Input Arguments:
%   (Required)
%   AABB               Axis-aligned bounding boxes as a 3-D array with one page for each track and
%                      [x0,x1; y0,y1; z0,z1] as respective box coordinates.
%                      (3,2,:) double
%
%   (Optional)
%   timeWindow         Temporal "thickness" of the bounding boxes. This extends them in z by the
%                      specified number of frames, preceding AND following the original temporal
%                      extent of each trajectory (default: all bounding boxes cover the entire frame
%                      range, thus calculating only their spatial overlap).
%                      (1,1) double
%
% Output Arguments:
%   overlapGroups      Grouped indices as a cell array. Each cell holds a self-contained collection
%                      of overlapping bounding boxes (indices = pages of AABB). Cells with singular
%                      entries contain the non-overlapping remnants.
%                      (:,1) cell
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Implementation of the time window along the procedure described in:
%
% Wallis et al., Molecular Videogaming: Super-Resolved Trajectory-Based
% Nanoclustering Analysus Using Spatio-Temporal Indexing,
% bioRxiv 2021.09.08.459552,
% doi: https://doi.org/10.1101/2021.09.08.459552.
%
% Vectorization in MATLAB is very effective here. Instead of checking the overlap of all AABBs with
% each other in a loop, the use of array comparisons and the consecutive graph-based analysis to
% find groups is quite fast.
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
% Initial release: 2021-12-03
% Last revision: 2023-05-03

%% Function argument validation
arguments
    AABB (3,2,:) double
    timeWindow (1,1) = inf
end

%% Analysis of each dimension
% Two lines A0-A1 and B0-B1 in one-dimensional space overlap when A0 <= B1 and B0 <= A1. Using
% vectorization, this can be calculated by comparing all x0 coordinates with all x1 coordinates.
x0 = squeeze(AABB(1,1,:));
x1 = squeeze(AABB(1,2,:));
isOverlappingX = x0' <= x1;

% Using the example from above, this results in a logical square matrix, with the lower left
% triangle now containing the A0 <= B1 info and the upper right triangle containing the B0 <= A1
% info, with the main diagonal as symmetry axis for the line indices. To overlap, both conditions
% must be true, so the logical matrix is compared with its transpose.
isOverlappingX = isOverlappingX.' & isOverlappingX;

% The diagonal will then be irrelevant (self-comparison), but both triangles now hold the overlap
% information. Further calculations will thus only use one triangle.

% This is also done for the other dimensions:
y0 = squeeze(AABB(2,1,:));
y1 = squeeze(AABB(2,2,:));
isOverlappingY = y0' <= y1;
isOverlappingY = isOverlappingY.' & isOverlappingY;

% Calculate the temporal overlap only if necessary.
if ~isinf(timeWindow)
    z0 = squeeze(AABB(3,1,:));
    z1 = squeeze(AABB(3,2,:));
    
    % Enlarge the time range by |timeWindow|. This will extend some ranges beyond the experimental
    % data acquisition, but fortunately this is irrelevant for the logical comparison to find
    % overlapping ranges.
    z0 = z0 - round(timeWindow/2);
    z1 = z1 + round(timeWindow/2);
    
    isOverlappingZ = z0' <= z1;
    isOverlappingZ = isOverlappingZ.' & isOverlappingZ;
else
    % Set z to all true (overlap is only determined by the x and y coordinates)
    isOverlappingZ = true(size(isOverlappingX));
end

%% Grouping

% For cuboids in 3D space to overlap, they must overlap in each dimension
isOverlappingXYZ = isOverlappingX & isOverlappingY & isOverlappingZ;

% To find the indices of overlapping cuboids, the >graph< function offers a convenient solution,
% using |isOverlappingXYZ| as an adjacency matrix (using only one of the triangles and disregarding
% the main diagonal)
G = graph(isOverlappingXYZ, 'upper', 'omitselfloops');

% Using the graph object functionality, groups are easily found. Each overlapping cuboid is a 'node'
% and 'connected nodes' are not only directly overlapping cuboids, but also those via intermediates.
% As an example: A only overlaps with B and B overlaps with A and C. Then A, B and C are part of the
% same cluster (A is connected to C via B). With >conncomp<, connected nodes can be returned as a
% cell array.
overlapGroups = conncomp(G,'OutputForm','cell')';

end
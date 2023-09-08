function AABB = convertTracksToAABB(trackArray, radiusFactor)
%convertTracksToAABB converts particle trajectories to axis-aligned bounding boxes.
%
% Syntax:
%   AABB = convertTracks2AABB(tracks)
%   AABB = convertTracks2AABB(tracks, radiusFactor)
%
% Input Arguments:
%   (Required)
%   trackArray         Particle trajectories as a cell array. Each cell (track) has at least three
%                      columns, in the order [frame, x-coordinate, y-coordinate] and at least three
%                      rows (localizations).
%                      (:,1) cell  
%
%   (Optional)
%   radiusFactor       Factor to approximate the radius of the convex hull for the idealized
%                      bounding boxes (default: 1.3).
%                      (1,1) double
%
% Output Arguments:
%   AABB               Axis-aligned bounding boxes as a 3-D array with one page for each track and
%                      [x0,x1; y0,y1; z0,z1] as respective box coordinates.
%                      (3,2,:) double
%
% Other required m-files: none
% Subfunctions: mustBeTracksArray
% Additional required MATLAB products: none
%
% Notes:
% Coded independently along the procedure described in the publication:
%
% Wallis et al., Molecular Videogaming: Super-Resolved Trajectory-Based
% Nanoclustering Analysus Using Spatio-Temporal Indexing,
% bioRxiv 2021.09.08.459552,
% doi: https://doi.org/10.1101/2021.09.08.459552.
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
    trackArray (:,1) cell {mustBeTracksArray}
    radiusFactor (1,1) double = 1.3
end

%% Main
centroidArray = zeros(numel(trackArray),2);
frameRange = zeros(numel(trackArray),2);
areaArray = zeros(numel(trackArray),1);

for iTrack = 1:numel(trackArray)
    % Determine the centroid coordinates.
    centroidArray(iTrack,:) = mean(trackArray{iTrack}(:,2:3),1);
    
    % Get the frame range of each track.
    [frameRange(iTrack,1), frameRange(iTrack,2)] = bounds(trackArray{iTrack}(:,1),1);
    
    % Calculate the area of the convex hull.
    [~, areaArray(iTrack)] = convhull(double(trackArray{iTrack}(:,2:3)));
end

% The circle radii are based on the calculated areas, multiplied with the user-defined radius factor.
radii = sqrt(areaArray/pi) * radiusFactor;

% The start/end coordinates of the bounding squares can now be calculated by subtracting/adding the
% radii from/to the centroid coordinates.
AABB = zeros(3, 2, numel(trackArray));
AABB(1,1,:) = centroidArray(:,1) - radii;
AABB(1,2,:) = centroidArray(:,1) + radii;
AABB(2,1,:) = centroidArray(:,2) - radii;
AABB(2,2,:) = centroidArray(:,2) + radii;

% The third dimension comprises the frame ranges.
AABB(3,:,:) = permute(frameRange, [3 2 1]);

end

%% Subfunctions
function mustBeTracksArray(tracks)
% argument validation function
if ~all(cellfun("size",tracks,2) > 2)
    eid = 'Size:wrongDimensions';
    msg = 'Each entry in tracks array must have at least three columns [frame, x-coordinate, y-coordinate].';
    throwAsCaller(MException(eid,msg))
elseif ~all(cellfun("size",tracks,1) > 2)
    eid = 'Size:wrongDimensions';
    msg = 'Each track must have at least three localizations.';
    throwAsCaller(MException(eid,msg))
end

end
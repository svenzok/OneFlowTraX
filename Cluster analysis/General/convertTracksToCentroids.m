function centroidList = convertTracksToCentroids(trackArray)
%convertTracksToCentroids converts particle trajectories to their respective centroid coordinates.
%
% Syntax:
%   centroidList = convertTracksToCentroids(trackArray)
%
% Input Arguments:
%   (Required)
%   trackArray         Particle trajectories as a cell array. Each cell (track) has at least three
%                      columns, in the order [frame, x-coordinate, y-coordinate]. Only the
%                      coordinates are used here, but they have to be at these positions.
%                      (:,1) cell  
%
% Output Arguments:
%   centroidList       List of centroid coordinates in the format [x-coordinate, y-coordinate],
%                      rows corresponding to the cells in |trackArray|.
%                      (:,2) double
%
% Other required m-files: none
% Subfunctions: mustBeTracksArray
% Additional required MATLAB products: none
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
% Last revision: 2023-05-19

%% Function argument validation
arguments
    trackArray (:,1) cell {mustBeTracksArray}
end

%% Main
centroidList = cell2mat(cellfun(@(x) [mean(x(:,2)), mean(x(:,3))], ...
    trackArray, "UniformOutput", false));

end

%% Subfunctions
function mustBeTracksArray(tracks)
% argument validation function
if ~all(cellfun("size",tracks,2) > 2)
    eid = 'Size:wrongDimensions';
    msg = 'Each entry in tracks array must have x and y data in columns 2 and 3.';
    throwAsCaller(MException(eid,msg))
end

end
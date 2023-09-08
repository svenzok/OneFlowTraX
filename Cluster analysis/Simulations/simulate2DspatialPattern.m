function pointList = simulate2DspatialPattern(limits, nPoints, Options)
%simulate2DspatialPattern simulates several types of 2-D spatial patterns
%
% Syntax:
%   pointList = simulate2DspatialPattern(limits, nPoints)
%   pointList = simulate2DspatialPattern(limits, nPoints, ...
%                'model', 'dispersed', 'radiusDispersed', 50)
%
% Input Arguments:
%   (Required)
%   limits             Area limits in the format
%                      [xmin, xmax, ymin, ymax].
%                      The coordinates of the simulated points will have
%                      the same length unit.
%                      (1,4) double
%
%   nPoints            Number of simulated points.
%                      (1,1) double
%
% Name-Value Pair Input Arguments (Options):
%   (Optional)
%   model              Type of spatial pattern, either complete spatial
%                      randomness (CSR, default), clustered or dispersed.
%                      (1,:) char {'CSR', 'clustered', 'dispersed'}
%                      
%   nClusters          Total number of clusters (default: 100).
%                      (1,1) double
%
%   distClusters       Distribution of clusters, either random (default) or
%                      dispersed (prevents overlapping)
%                      (1,:) char {'random', 'dispersed'}
%
%   distInClusters     Distribution type of points inside clusters, either
%                      they are randomly placed inside a circular disk
%                      (default) or follow a 2-D normal distribution.
%                      (1,:) char {'disk', 'normal'}
%
%   diameterClusters   Diameter of each circular cluster, or FWHM if the
%                      clusters are modeled as 2-D normal distribution
%                      (default: 100)
%                      (1,1) double
%
%   fractInClusters    Fraction (0-1) of points that will be placed into
%                      clusters (default: 0.8). The remaining points will
%                      be modeled as noise.
%                      (1,1) double
%
%   noiseOutClusters   Determines how to model the noise outside the
%                      clusters, either with CSR (default) or dispersed.
%                      (1,:) char {'CSR', 'dispersed'}
%
%   radiusDispersed    Minimal distance of two points (or cluster centers)
%                      in a dispersed spatial pattern, also for dispersed
%                      noise in the clustered pattern (default: 20).
%                      (1,1) double
%
% Output Arguments:
%   pointList          List of simulated point coordinates in the format
%                      [x-coordinate, y-coordinate].
%                      (:,2) double
%
% Other required m-files: none
% Subfunctions: placeRandomPoints, placeClusteredPoints, placeDispersedPoints
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.2
%
% Notes:
% The placement of points or clusters in a dispersed pattern may not always
% be possible due to spatial constraints. In this case, an error message
% will be generated.
% A long input list of name-value pairs may be cumbersome for calling a
% function, an alternative would be passing a structure array with its
% fields and values as the name-value pairs. To make it work with the
% arguments block, use
% tmp = reshape(namedargs2cell(OptionsStructure),2,[])
% and pass tmp{:} to the function, as {:} generates a comma-separated list
% that is readable as name-value pairs.
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
% Last revision: 2022-12-01

%% Function argument validation
arguments
    limits (1,4) double
    nPoints (1,1) double
    Options.model {mustBeMember(Options.model,["CSR","clustered","dispersed"])} = 'CSR'
    Options.nClusters (1,1) double = 100;
    Options.distClusters {mustBeMember(Options.distClusters,["random","dispersed"])} = 'random'
    Options.distInClusters {mustBeMember(Options.distInClusters,["disk","normal"])} = 'disk'
    Options.diameterClusters (1,1) double = 100;
    Options.fractInClusters (1,1) double = 0.8;
    Options.noiseOutClusters {mustBeMember(Options.noiseOutClusters,["CSR","dispersed"])} = 'CSR'
    Options.radiusDispersed (1,1) double = 20;
end

%% Main

switch Options.model
    case 'CSR'
        pointList = placeRandomPoints(nPoints, limits);

    case 'clustered'
        % points in each cluster
        nPerCluster = round(Options.fractInClusters * ...
                            nPoints / Options.nClusters);

        % pre-allocate cell array (one cell for each cluster) and
        % generate the points for each of them
        clusterPointArray = cell(Options.nClusters,1);
        for iCluster = 1:Options.nClusters
            clusterPointArray{iCluster} = placeClusteredPoints(nPerCluster, ...
                Options.distInClusters, Options.diameterClusters);
        end
        
        % generate either random or dispersed cluster origins
        if strcmp(Options.distClusters,'random')
            clusterOriginList = placeRandomPoints(Options.nClusters, limits);

        elseif strcmp(Options.distClusters,'dispersed')
            % Place the centers of disk-like clusters at least two
            % diameters apart. Place the 'softer', normal-distributed
            % clusters at least four times their FWHM apart.
            if strcmp(Options.distInClusters,'disk')
                radiusDispersedCluster = 2 * Options.diameterClusters;
            elseif strcmp(Options.distInClusters,'normal')
                radiusDispersedCluster = 4 * Options.diameterClusters;
            else
                % empty
            
            end

            clusterOriginList = placeDispersedPoints(Options.nClusters, ...
                limits, radiusDispersedCluster);
        else
            % empty
        
        end
        
        % Unpack the cluster points.
        clusterPointList = cell2mat(clusterPointArray);

        % Shift the cluster points to their respective positions.
        clusterPointList = clusterPointList + ...
            repelem(clusterOriginList, nPerCluster, 1);

        % Add random or dispersed points as noise.
        nPointsNoise = nPoints - size(clusterPointList,1);
      
        if strcmp(Options.noiseOutClusters,'CSR')
             % Noise points might also be placed into the cluster regions,
             % this is usually negligible.
             noisePointList = placeRandomPoints(nPointsNoise, limits);

        elseif strcmp(Options.noiseOutClusters,'dispersed')
            % Dispersed noise points take into account the existing cluster
            % points.
             noisePointList = placeDispersedPoints(nPointsNoise, limits, ...
                 Options.radiusDispersed, clusterPointList);

        else
            % empty

        end

        % Delete cluster coordinates beyond the limits. This may result in
        % a smaller number of output points than originally specified by
        % |nPoints|.
        isOutsideX = (clusterPointList(:,1) < limits(1) | ...
                      clusterPointList(:,1) > limits(2));
        isOutsideY = (clusterPointList(:,2) < limits(3) | ...
                      clusterPointList(:,2) > limits(4));
        clusterPointList(isOutsideX | isOutsideY,:) = [];

        % Combine clustered points and noise
        pointList = [clusterPointList; noisePointList];

    case 'dispersed'
        pointList = placeDispersedPoints(nPoints, limits, ...
            Options.radiusDispersed);

    otherwise
        % empty

end

end

%% Subfunctions
function pointList = placeRandomPoints(nPoints,limits)
pointList =...
    [rand(nPoints,1) * (limits(2) - limits(1)) + limits(1), ...
     rand(nPoints,1) * (limits(4) - limits(3)) + limits(3)];
end

function pointList = placeClusteredPoints(nPerCluster,distInClusters,diameterClusters)
% placeClusteredPoints generates a cluster of points around position [0,0]

switch distInClusters
    case "disk"
            r = rand(nPerCluster,1);
            theta = 2*pi*rand(nPerCluster,1);
            x = sqrt(r).*cos(theta);
            y = sqrt(r).*sin(theta);
            pointList = diameterClusters/2 * [x y];

    case "normal"
        % Here, |diameterClusters| is the FWHM.
        covariance = (diameterClusters/(2*sqrt(2*log(2))))^2;
        pointList = mvnrnd([0 0], [covariance covariance], nPerCluster);

    otherwise

end
end

function pointList = placeDispersedPoints(nPoints,limits,radiusDispersed,existingPoints)
% placeDispersedPoints tries to put points into a pre-defined area while
%                      keeping a minimal distance between any of them.

% squared exclusion distance for calculations
sqaredDistance  = radiusDispersed ^ 2;

% adjust if there are existing points
if nargin == 4
    x = [existingPoints(:,1); zeros(nPoints,1)];
    y = [existingPoints(:,2); zeros(nPoints,1)];

    nPoints = nPoints + size(existingPoints,1);
    nValid  = size(existingPoints,1);

else
    x = zeros(nPoints, 1);
    y = zeros(nPoints, 1);
    nValid  = 0;

end

% count loops as an exit condition when no more points can be placed into
% the area defined by |limits|
loopCounter = 1;

% Try to find points that are sufficiently separated from all other points
% that have been stored so far. Stop when the number of valid points equals
% the target number or when the loop counter exceeds one million.
while nValid < nPoints && loopCounter < 1E6
    xNew = rand * (limits(2) - limits(1)) + limits(1);
    yNew = rand * (limits(4) - limits(3)) + limits(3);
    if all(((x(1:nValid) - xNew).^2 +...
            (y(1:nValid) - yNew).^2) >= sqaredDistance) == true
        % new random point ist valid, so append it
        nValid    = nValid + 1;
        x(nValid) = xNew;
        y(nValid) = yNew;
    end
    loopCounter = loopCounter + 1;
end

if nValid ~= nPoints
    error(['Error when placing points in a dispersed pattern.',...
           '\nOnly %d of %d points possible at this density.\n'],...
           nValid, nPoints);
end

pointList = [x(1:nValid), y(1:nValid)];

end
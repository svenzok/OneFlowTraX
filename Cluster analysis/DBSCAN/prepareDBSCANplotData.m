function colorValues = prepareDBSCANplotData(clusterIndexList, isCore, colorScheme)
% prepareDBSCANplotData sets color values (for colormaps) for points categorized by DBSCAN.
%
% Syntax:
%   colorValues = prepareDBSCANplotData(classificationList, colorScheme)
%
% Input Arguments:
%   (Required)
%   clusterIndexList   Cluster assignments (first output of >dbscan<). 
%                      (:,1) double
%
%   isCore             Denotes if a point is a core point (second output of >dbscan<).
%                      (:,1) logical
%
%   colorScheme        Defines how to color the individual points:
%                      'randomX': all points in one cluster get a random value from 1 to X, with
%                                 X = 4,6,8 or 10. Noise points are set to 0.
%                      'noise/cluster': noise points are set to 0, cluster points to 1.
%                      'noise/core/border': |colorValues| is just a copy of the second column of
%                                           |classificationList|.
%                      (1,:) char | (1,1) string
%                      {'random4', 'random6', 'random8', 'random10', ...
%                       'noise/cluster', 'noise/core/border'}
%
% Output Arguments:
%   colorValues        Color data for each point.
%                      (:,1) double
%
% Other required files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
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
% Initial release: 2023-04-02
% Last revision: 2023-04-02

%% Function argument validation
arguments
    clusterIndexList (:,1) double
    isCore (:,1) logical
    colorScheme (1,:) char {mustBeMember(colorScheme, ...
       {'random4', 'random6', 'random8', 'random10', 'noise/cluster', 'noise/core/border'})}
end

%% Main

if contains(colorScheme, 'random')
    nCycleCols = sscanf(colorScheme, 'random%u');
    colorScheme = 'random';
end

% Process according to the chosen color scheme.
switch colorScheme
    case 'random'
        % Randomly assign values from 1:|nCycleCols| to cluster points.
        colorValues = zeros(size(clusterIndexList));
        for iCluster = 1:max(clusterIndexList)
            colorValues(clusterIndexList == iCluster) = randi(nCycleCols);
        end

    case 'noise/cluster'
        % Just discriminate between noise and cluster points.
        colorValues = double(clusterIndexList ~= -1);

    case 'noise/core/border'
        % Assign: noise = 0, core = 1 and border = 2.
        isNoise = clusterIndexList == -1;
        isBorder = ~isCore & ~isNoise;
        colorValues = zeros(size(isNoise));
        colorValues(isCore) = 1;
        colorValues(isBorder) = 2;
end

end
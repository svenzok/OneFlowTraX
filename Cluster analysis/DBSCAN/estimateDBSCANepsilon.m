function [epsilon, PlotData] = estimateDBSCANepsilon(pointList, minPts)
% estimateDBSCANepsilon estimates the epsilon value for DBSCAN.
%
% Syntax:
%   [epsilon, PlotData] = estimateDBSCANepsilon(pointList, minPts)
%
% Input Arguments:
%   (Required)
%   pointList          List of point coordinates in the format
%                      [x-coordinate, y-coordinate].
%                      (:,2) double
%
%   minPts             Minimum number of points that are contained in a
%                      circle around one point with radius epsilon
%                      (including the central point), so that the central
%                      point is classified as a core point in a cluster.
%                      minPts is usually set to two times the
%                      dimensionality of the problem (= 4 for 2-D).
%                      (1,1) double
%
% Output Arguments:
%   epsilon            Estimated epsilon value.
%                      (1,1) double
%
%   PlotData           Data for plotting rectangles marking the noise
%                      and cluster regions and the dividing line, as a
%                      structure with fields:
%                      .noiseRegion: Sub-fields .XData and .YData hold the
%                                    values for a patch object
%                      .clusterRegion: Sub-fields .XData and .YData hold
%                                      the values for a patch object
%                      .splitLine: Sub-fields .XData and .YData hold the
%                                  values for line object
%                      .kDistance: Sub.fields .XData and .YData hold the values for a line object
%                      (1,1) struct
%
% Other required m-files: calculatePointsToLineDistance2D
% Subfunctions: none
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.2
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
% Last revision: 2021-12-26

%% Function argument validation
arguments
    pointList (:,2) double
    minPts (1,1) double
end

%% Main

% Calculate the distances for the k nearest neighbors (kNN) of every point
% and use minPts as the value for k.
kNN = pdist2(pointList, pointList, 'euclidean', 'Smallest', minPts);

% sorted distances
dataY = sort(kNN(end,:));

% index
nPoints = length(kNN);
dataX = 1:nPoints;

% The "knee" of this plot (the largest curvature) is a good estimate for
% the threshold between cluster and noise. Construct a line from the first
% to last datapoint and collect the shortest distances for each datapoint
% to this line (i.e., each point is the endpoint of a vector that starts at
% the line perpendicular to the line).
lineVector = [dataX(1), dataX(end); dataY(1), dataY(end)];
distanceList = ...
    calculatePointsToLineDistance2D([dataX', dataY'], lineVector);

% The largest curvature is at the point that has the largest perpendicular
% distance to the line.
[~,splitValueID] = max(distanceList);
epsilon = dataY(splitValueID);

% Store data for plotting two squares, a line and the k-distance graph.
PlotData.clusterRegion.XData = dataX([1 nPoints nPoints 1]);
PlotData.clusterRegion.YData = dataY([1 1 splitValueID splitValueID]);
PlotData.noiseRegion.XData = PlotData.clusterRegion.XData;
PlotData.noiseRegion.YData = ...
    dataY([splitValueID splitValueID nPoints nPoints]);
PlotData.splitLine.XData = dataX([1 nPoints]);
PlotData.splitLine.YData = dataY([splitValueID splitValueID]);
PlotData.kDistance.XData = dataX;
PlotData.kDistance.YData = dataY;

end
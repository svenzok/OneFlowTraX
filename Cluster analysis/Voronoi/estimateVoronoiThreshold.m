function [splitValue, histDataSample, histDataMC] = estimateVoronoiThreshold(pointList, method, nIterations, Options)
% estimateVoronoiThreshold uses Monte Carlo simulations to estimate the threshold value separating
%                          complete spatial randomness and clustering using Voronoi data.
%
% Syntax:
%   splitValue = estimateVoronoiThreshold(pointList, method, nIterations)
%   [splitValue, histDataSample, histDataMC] = ...
%   estimateVoronoiThreshold(pointList, method, nIterations)
%
% Input Arguments:
%   (Required)
%   pointList          List of point coordinates in the format [x-coordinate, y-coordinate].
%                      (:,2) double
%
%   method             Method to discriminate between clusters and noise, using either the Voronoi
%                      cell areas or the normalized first-rank densities.
%                      (1,:) char {'area', 'density'}
%
%   nIterations        Number of Monte Carlo simulations
%                      (1,1) double
%
% Name-Value Pair Input Arguments:
%   (Optional)
%   updateHandle       Handle to a graphic object with a .Text property. This can be used to update
%                      the simulation progress to an external app.
%                      (1,1) matlab.ui.control.Label (other classes are possible)
%
% Output Arguments:
%   splitValue         The value at which the histograms of the sample data and simulated complete
%                      spatial randomness (CSR) data are intersecting, useful to find Voronoi cells
%                      belonging to a cluster (their area has to be smaller, their normalized 
%                      first-rank density has to be larger than this value).
%                      (1,1) double
%
%   histDataSample     Data points to plot the sample data histogram in the format
%                      [bin centers, counts].
%                      (:,2) double
%
%   histDataMC         Data points to plot the Monte Carlo simulation results for CSR in the format
%                      [bin centers, mean counts, min. counts, max. counts].
%
% Other required m-files: calculateVoronoiArea, calculateFirstRankDensity
% Subfunctions: findIntersection
% Additional required MATLAB products:
% - Statistics and Machine Learning Toolbox 12.2
%
% Notes:
%
% Coded along the procedure (and using the area histogram parameters) described in the publication:
%
% Andronov, L. et al. ClusterViSu, a method for clustering of protein complexes by Voronoi
% tessellation in super-resolution microscopy. Sci Rep 6, 24084 (2016).
% https://doi.org/10.1038/srep24084.
%
% The histogram parameters for estimating alpha (using the normalized first-rank density) were
% determined by simulations of data sets with different densities, cluster sizes and numbers.
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
% Initial release: 2021-12-07
% Last revision: 2023-05-30

%% Function argument validation
arguments
    pointList (:,2) double
    method (1,:) char {mustBeMember(method,["area","density"])}
    nIterations (1,1)
    
    Options.updateHandle (1,1) {mustHaveTextProperty(Options.updateHandle)} = [];
end

%% Main

% Calculate the Delaunay triangulation and convert to Voronoi vertices and cells, intersecting them
% with an axis-aligned bounding box to deal with out-of-limits vertices.
DT = delaunayTriangulation(pointList);
[DT, voronoiVertices, voronoiCells] = getTrimmedVoronoi(DT.Points, 'type', 'intersect', 'boundary', 'AABB');

% Get the data limits for subsequent simulations.
[xmin, xmax] = bounds(DT.Points(:,1));
[ymin, ymax] = bounds(DT.Points(:,2));
limits = [xmin, xmax, ymin, ymax];

 % Calculate the Voronoi cell areas of the sample data.
voronoiAreaList = calculateVoronoiArea(voronoiVertices, voronoiCells);

switch method
    case 'area'
        % Calculate the area histogram. Andronov et al. use 0 to 4 times the median value of the areas as the range for
        % the histograms, then use the Rice rule for the number of bins. We apply their method to calculate the bin
        % width, but extend the range (to include some cases where the intersection is beyond that original range).
        nBins = round(2 * sum(~isnan(voronoiAreaList)) ^ (1/3));
        areaLimit = 4 * median(voronoiAreaList, 'omitnan');
        binWidth = areaLimit / nBins;
        [counts, binEdges] = histcounts(voronoiAreaList, 0:binWidth:areaLimit*3);
        nBins = numel(binEdges)-1;

        % Start Voronoi Monte Carlo simulations (for complete spatial randomness).
        countsMC = zeros(nIterations, nBins);
        for iIteration = 1:nIterations
            % Simulate a CSR pattern with the same limits and number of points as the sample data.
            pointListMC = simulate2DspatialPattern(limits, length(voronoiAreaList));

            % Voronoi tesselation
            DT = delaunayTriangulation(pointListMC);
            [~, voronoiVertices, voronoiCells] = ...
                getTrimmedVoronoi(DT.Points, 'type', 'intersect', 'boundary', 'AABB');

            % calculation of areas and histogram
            voronoiAreaListMC = calculateVoronoiArea(voronoiVertices, voronoiCells);
            countsMC(iIteration, :) = histcounts(voronoiAreaListMC, binEdges);

            % Update a text message if needed.
            if ~isempty(Options.updateHandle)
                if nIterations - iIteration == 0
                    Options.updateHandle.Text = "Estimating averaged threshold value ..."; drawnow
                else
                    Options.updateHandle.Text = sprintf("%u simulations left. Please wait ...", ...
                    nIterations - iIteration); drawnow
                end
            end
        end

    case 'density'
        % Calculate the normalized first-rank densities of the sample.
        normFirstRankDensityList = calculateNormFirstRankDensity(DT, voronoiAreaList);

        % Calculate the density histogram. This estimation method is not used in the publication for the SR-Tesseler
        % software, but typical values for alpha in literature are set to around 2. Correspondingly, the histogram will
        % cover values between 1 (the average density) and a value of ten times the median value.
        nBins = 50;
        densityLimit = 10 * median(normFirstRankDensityList, 'omitnan');
        isInsideLimits = normFirstRankDensityList > 1 & normFirstRankDensityList <= densityLimit;
        [counts, binEdges] = histcounts(normFirstRankDensityList(isInsideLimits), nBins);

        % Voronoi Monte Carlo simulations (complete spatial randomness)
        countsMC = zeros(nIterations, numel(binEdges) - 1);
        for iIteration = 1:nIterations
            % simulate CSR pattern with the same limits and number of
            % points as the sample data
            pointListMC = simulate2DspatialPattern(limits, length(normFirstRankDensityList));

            % Voronoi tesselation
            DTMC = delaunayTriangulation(pointListMC);
            [vMC,cMC] = voronoiDiagram(DT);

            % calculation of densities and histogram
            voronoiAreaListMC = calculateVoronoiArea(vMC, cMC);
            normFirstRankDensityListMC = calculateNormFirstRankDensity(DTMC, voronoiAreaListMC);
            countsMC(iIteration,:) = histcounts(normFirstRankDensityListMC, binEdges);

            % Update a text message if needed.
            if ~isempty(Options.updateHandle)
                if nIterations - iIteration == 0
                    Options.updateHandle.Text = "Estimating averaged threshold value ..."; drawnow
                else
                    Options.updateHandle.Text = sprintf("%u simulations left. Please wait ...", ...
                    nIterations - iIteration); drawnow
                end
            end
        end

    otherwise
        % empty
end

meanCountsMC = mean(countsMC);
[minCountsMC, maxCountsMC] = bounds(countsMC);

% Shift and crop the edge vector so that the data now corresponds to the center of each bin (using
% >conv<);
binCenters = conv(binEdges, [0.5 0.5], 'valid');

% Find the intersection of the sample and CSR histograms to estimate the segmentation threshold.
splitValue = findIntersection(binCenters, counts, meanCountsMC);

% Export the histogram data for plotting.
histDataSample = [binCenters; counts].';
histDataMC = [binCenters; meanCountsMC; minCountsMC; maxCountsMC].';

end

%% Subfunctions
function mustHaveTextProperty(updateHandle)
if isgraphics(updateHandle) && isprop(updateHandle,'Text')
    % Graphic object has the .Text property, all is well.
else
    eid = 'TextProperty:notFound';
    msg = 'The handle to the graphic object must have a .Text property.';
    throwAsCaller(MException(eid,msg))
end
end


function xIntersection = findIntersection(xData,yData1,yData2)
% Finds the first position on the x-axis where two lines intersect.
%
% Adapted from Douglas Schwarz (2021). Fast and Robust Curve Intersections
% (https://www.mathworks.com/matlabcentral/fileexchange/11837-fast-and-robust-curve-intersections),
% MATLAB Central File Exchange. Retrieved December 6, 2021.
%
% Here, the matrix notation for linear solving is also used, but without
% error handling or speed improvements as the input vectors are usually
% short.

x = xData(:); y1 = yData1(:); y2 = yData2(:);
xy1 = [x y1]; dxy1 = diff(xy1);
xy2 = [x y2]; dxy2 = diff(xy2);

n = length(x) - 1;
T = zeros(4,n);
AA = zeros(4,4,n);
AA([1 2],3,:) = -1;
AA([3 4],4,:) = -1;
AA([1 3],1,:) = dxy1.';
AA([2 4],2,:) = dxy2.';
B = -[x x y1 y2].';

warning('off','MATLAB:singularMatrix')
for k = 1:n
    [L,U] = lu(AA(:,:,k));
    T(:,k) = U\(L\B(:,k));
end
warning('on','MATLAB:singularMatrix')

isInRange = (T(1,:) >= 0 & T(2,:) >= 0 & T(1,:) < 1 & T(2,:) < 1).';
x0 = T(3,isInRange).';

% set output to NaN if no intersections can be found
if isempty(x0)
    xIntersection = NaN;
else
    xIntersection = x0(1);
end

end
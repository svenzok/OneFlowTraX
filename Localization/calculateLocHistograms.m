function [histogramData, axesLimits] = calculateLocHistograms(data, dataFiltered)
% calculateLocHistograms calculates the localization histograms for the OneFlowTraX app.
%
% Syntax:
%   histogramData = calculateLocHistograms(data, limits)
%
% Input Arguments:
%   (Required)
%   data               Nx3 list of localization data, with columns PSF size, photons and localization precision.
%                      (:,3) double
%
%   dataFiltered       As |data|, but filtered with user-defined histogram limits.
%                      (:,3) double
%
% Output Arguments:
%   histogramData      Histogram data, one cell for each set (same order as in |data|). In each set, the first column
%                      are the centers of the histogram bins, the second column is the data for the original histogram
%                      and the third column is the data for the histogram when user-defined limits are applied.
%                      (1,3) cell
%
%   axesLimits         Axes limits for the individual histograms, in the row order xmin, xmax, ymax. The columns
%                      correspond to |data|.
%                      (3,3) double
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
%
% Tested: MATLAB Version: 9.13.0.2166757 (R2022b)
%	      Microsoft Windows 10 Enterprise Version 10.0 (Build 19045)
%
% Author: Sven zur Oven-Krockhaus
%	      Institute of Physical and Theoretical Chemistry
%	      University of Tuebingen, Tuebingen, Germany
% E-mail: sven.zur-oven-krockhaus@uni-tuebingen.de
%
% GNU placeholder
%
% Initial release: 2023-06-15
% Last revision: 2023-08-09

%% Function argument validation
arguments
    data (:,3) double
    dataFiltered (:,3) double
end

%% Main

% Pre-allocate the output.
histogramData = cell(1,3);
axesLimits = zeros(3);

for iHist = 1:3
    % Estimation of x-axis limits and binsize.
    [~,L,U] = isoutlier(data(:,iHist),'median','ThresholdFactor',2.5);

    % Determine the appropriate rounding precision.
    prc = round(log10(U))-1;

    % Calculate the axis limits and bin size.
    xmin = max(round(L,-prc),0);
    xmax = round(1.25*U,-prc);
    binsize = round(0.01*(xmax - xmin),-prc+1);

    while binsize == 0 || (xmax - xmin)/binsize < 30
        prc = prc - 1;
        binsize = round(0.01*(xmax - xmin),-prc+1);
    end

    % Histogram for the original data. Always set the left edge to 0.
    [N1, edges] = histcounts(data(:,iHist), 0:binsize:xmax);

    % Histogram for the filtered data.
    N2 = histcounts(dataFiltered(:,iHist), edges);

    % Combine with the histogram bin centers and store.
    binCenters = edges(1:end-1) + mean(diff(edges))/2;
    histogramData{iHist} = [binCenters', N1', N2'];

    % Store the axis limits (leave some room between the top of the y axis and the data).
    axesLimits(1,iHist) = edges(1);
    axesLimits(2,iHist) = edges(end);
    axesLimits(3,iHist) = 1.05*max(N1);
end

end
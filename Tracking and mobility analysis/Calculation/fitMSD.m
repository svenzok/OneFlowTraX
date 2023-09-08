function [indivMSDfits, fitCoeffs] = fitMSD(MSDresultsArray, frameRate, Options)
% fitMSD fits MSD curves to determine diffusion coefficients.
%
% Syntax:
%   [indivMSDfits, fitCoeffs] = fitMSD(MSDresultsArray, frameRate, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   MSDresultsArray    Cell array, one cell for each track, with columns:
%                      1   time lag (s)
%                      2   MSD (square micrometers)
%                      3   sdv (square micrometers)
%                      4   N (number of datapoints)
%                      For fitting only one MSD curve (e.g., to fit the results of a mean MSD
%                      calculation), provide the respective MSD data as a cell by enclosing the
%                      variable in {} brackets.
%                      (:,1) cell
%
%   frameRate          Camera acquisition frame rate in Hz.
%                      (1,1) double
%
% Name-Value Pair Input Arguments:
%   (Optional)
%   linFitPoints       Selects the points of the MSD curve that the linear regression will be
%                      applied to. If |linFitPoints| is a single number X, then the first X points
%                      of the MSD curve will be fitted. If |linFitPoints| is a 1x2 vector, then the
%                      fit will be done for an interval. For example, [2 4] will apply the linear
%                      regression to the data points 2, 3 and 4 of the MSD curve
%                      (default: [2 4]).
%                      (1,1)|(1,2) double
%
% Output Arguments:
%   indivMSDfits       Nx2 vector with columns D (µm²/s) and the adjusted R² value for each entry in
%                      |MSDresultsArray|. Invalid entries (the original tracks were too short or had
%                      gaps in the required fitting interval) are not fitted, indicated by a [0 0]
%                      pair in the corresponding row of |indivMSDfits|. This preserves the row order
%                      to connect a track with its corresponding fitting results.
%                      (:,2) double
%
%   fitCoeffs          Raw fit coefficients of the linear fits, with:
%                      y = fitCoeffs(1) + fitCoeffs(2)*x.
%                      (:,2) double
%
% Other required m-files: linearRegression
% Subfunctions: checkValidPoints
% Additional required MATLAB products: none
%
% Notes:
% This method only fits the first few points of the MSD curves to extract D, which is a common
% practice in the single particle tracking field. The calculated values are referred to differently
% in the literature (instantaneous, apparent, effective, short-term D or D2-4), but these terms have
% very different meanings, depending on the field of research. D2-4 is defined best, and is used as
% the default setting for this program. Here, D is only derived from the time lags 2-4, for which
% super- or subdiffusive behavior should not yet manifest. The first point is omitted, as this
% time interval shows complex behavior [1].
%
% [1] Kusumi, A. et al. Confined Lateral Diffusion of Membrane Receptors as Studied by Single
% Particle Tracking (Nanovid Microscopy). Effects of Calcium-induced Differentiation in Cultured
% Epithelial Cells. Biophys J 65, 2021-2040 (1993).
% https://doi.org/10.1016/s0006-3495(93)81253-0
% 
% Tested: MATLAB Version 9.12.0.1884302 (R2022a),
%	      Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% Author: Sven zur Oven-Krockhaus
%	      Institute of Physical and Theoretical Chemistry
%	      University of Tuebingen, Tuebingen, Germany
% E-mail: sven.zur-oven-krockhaus@uni-tuebingen.de
%
% GNU placeholder
%
% Initial release: 2022-10-27
% Last revision: 2023-03-23

%% Function argument validation
arguments
    MSDresultsArray (:,1) cell
    frameRate (1,1) double
    Options.linFitPoints (1,2) double {mustBePositive, checkValidPoints(Options.linFitPoints)} = [2 4]
end

%% Main

% Get the fit point indices. If only one number was specified, the arguments block will have
% expanded it to a vector, so start with 1.
if diff(Options.linFitPoints) == 0
    pointIDs = 1:Options.linFitPoints(2);
else
    pointIDs = Options.linFitPoints(1):Options.linFitPoints(2);
end

% Pre-allocate.
indivMSDfits = zeros(numel(MSDresultsArray), 2);
fitCoeffs = zeros(numel(MSDresultsArray), 2);

% We need to make sure that the specified time lags for the linear fit are present. It is sufficient
% to compare the n-th time point of each track (with n as the last point for the linear fit) with
% the expected time lag for an uninterrupted track at that position. A difference tolerance of 1E-5
% is used here to account for possible precision problems.
isValid = cellfun(@(x) abs(x(pointIDs(end)) - pointIDs(end) / frameRate) < 1E-5, MSDresultsArray);


% Loop through the valid MSD curves.
for iTrack = find(isValid')
    % The linear regression uses the lag time as X and the MSD as Y. Additionally, the number of
    % data points that were available for each X/Y pair (4th column of |MSDresultsArray|) is used as
    % weights. Note that division by 4 only produces correct diffusion coefficients if the original
    % tracks were acquired in 2D.
    [c, adjR2] = linearRegression(...
        MSDresultsArray{iTrack}(pointIDs,1),...
        MSDresultsArray{iTrack}(pointIDs,2),...
        MSDresultsArray{iTrack}(pointIDs,4));
    indivMSDfits(iTrack,:) = [c(2)/4, adjR2];
    fitCoeffs(iTrack,:) = c;
end

end

%% Subfunctions
function checkValidPoints(linFitPoints)

if diff(linFitPoints) == 0
    % User specified a 1x2 vector with the same numbers or a scalar (that was expanded to a 1x2
    % vector in the arguments block)
    if linFitPoints(2) < 3
        eid = 'linFitPoints:tooFew';
        msg = 'Linear regression is required to have at least three data points.';
        throwAsCaller(MException(eid,msg))
    else
        % Do nothing.
    end

else
    % User specified a 1x2 vector with different numbers or the default value is used.
    if diff(linFitPoints) < 0
        eid = 'linFitPoints:wrongOrder';
        msg = 'Second value of linFitPoints must be higher than the first';
        throwAsCaller(MException(eid,msg))

    elseif diff(linFitPoints) < 2
        eid = 'linFitPoints:tooFew';
        msg = 'Linear regression is required to have at least three data points.';
        throwAsCaller(MException(eid,msg))

    else
        % Do nothing.
    end
end

end
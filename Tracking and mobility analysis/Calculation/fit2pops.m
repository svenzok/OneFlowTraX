function [pdComp1, pdComp2, fractionComp1, errorCode] = fit2pops(logDlist)
% fit2pops attempts a two-component distribution fit with some error handling (OneFlowTrax app).
%
% Syntax:
%   [pdComp1, pdComp2] = fit2pops(logDlist)
%
% Input Arguments:
%   (Required)
%   logDlist           List of logarithmized diffusion coefficient values (usually between -5 and 1).
%                      (:,1) double
%
% Output Arguments:
%   pdComp1            Normal distribution object for the slower component.
%                      (1,1) probability distribution object
%
%   pdComp1            Normal distribution object for the faster component.
%                      (1,1) probability distribution object
%
%   fractionComp1      Fraction of the slower component.
%                      (1,1) double
%
%   errorCode          0 when all worked well, 1 when warnings were issued, 2 when failed.
%                      (1,1) double
%
% Other required m-files: [mList,pList] = matlab.codetools.requiredFilesAndProducts('fname.m')
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
% Initial release: 2023-03-01
% Last revision: 2023-08-08

%% Function argument validation
arguments
    logDlist (:,1) double
end

%% Main

% Disable warnings for non-converging fits.
warning('off','stats:gmdistribution:FailedToConverge');
warning('off','stats:gmdistribution:FailedToConvergeReps');
warning('off','stats:gmdistribution:IllCondCov');

% The covariances (sigmas) are always estimated to be 0.1 (if needed later).
s(:,:,1) = 0.1; s(:,:,2) = 0.1;

% Create an Options structure.
statOptions = statset('normfit');
statOptions.MaxFunEvals = 500;
statOptions.MaxIter = 1000;

% Use try/catch blocks to catch unexpected fitting fails.
try
    % First, let the internal algorithm use 5 different start parameter sets to find the best fit.
    GMModel = fitgmdist(logDlist, 2, 'CovarianceType', 'diagonal', 'Options', statOptions, 'Replicates', 5);
    errorCode = 0;
catch
    errorCode = 2;
end

if errorCode == 2
    % If the previous method failed, try again by setting the start parameters to commonly encountered values in
    % sptPALM experiments.
    StartP = struct('mu', [-3, -2], 'Sigma', s, 'ComponentProportion', [0.5, 0.5]);

    try
        GMModel = fitgmdist(logDlist, 2, 'Start', StartP, 'CovarianceType', 'diagonal', 'Options', statOptions);
        errorCode = 0;
    catch
        errorCode = 2;
    end
end

if errorCode == 2
    try
        % If the previous method failed, try again by estimating start parameters based on a smoothed histogram of the
        % data.
        [N, edges] = histcounts(logDlist,'BinWidth',0.1);
        Nsmooth = sgolayfilt(N,3,5);

        % Use the data gradient to find the locations of the two highest peaks.
        Ngrad = abs(gradient(Nsmooth));
        minHeight = max(Ngrad)/5;
        [~, idx] = findpeaks(Ngrad,'NPeaks',2,'MinPeakHeight',minHeight,'SortStr','descend');

        % Use these locations as starting values for the distribution fit.
        mu = edges(idx)';

        % Get the peak heights in the original histogram at these locations to guess initial component proportions.
        ph = Nsmooth(idx);
        pcomp = ph/sum(ph);
        StartP = struct('mu', mu, 'Sigma', s, 'ComponentProportion', pcomp);

        % Fit the distribution
        GMModel = fitgmdist(logDlist, 2, 'Start', StartP, 'CovarianceType', 'diagonal', 'Options', statOptions);
        errorCode = 0;
    catch
        errorCode = 2;
    end
end

if errorCode == 0
    % Check if the fit converged and adjust the error code accordingly.
    errorCode = 1 - GMModel.Converged;

elseif errorCode == 2
    % Pass empty outputs and return.
    [pdComp1, pdComp2, fractionComp1] = deal([]);
    return
end

% Sort the component parameters, so that the slower one comes first.
[mu, sortID] = sort(GMModel.mu);
sigma = sqrt(squeeze(GMModel.Sigma));
sigma = sigma(sortID);
compFractions = GMModel.ComponentProportion(sortID);
fractionComp1 = compFractions(1);

% Create normal distribution objects for both components.
pdComp1 = makedist('normal', mu(1), sigma(1));
pdComp2 = makedist('normal', mu(2), sigma(2));

% Re-enable warnings for non-converging fits.
warning('on','stats:gmdistribution:FailedToConverge');
warning('on','stats:gmdistribution:FailedToConvergeReps');
warning('on','stats:gmdistribution:IllCondCov');

end
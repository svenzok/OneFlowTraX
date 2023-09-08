function [c, adjR2] = linearRegression(X, Y, W)
% linearRegression performs a fast, optionally weighted linear regression.
%
% Syntax:
%   [c, adjR2] = linearRegression(X, Y)
%   [c, adjR2] = linearRegression(X, Y, W)
%
% Input Arguments:
%   (Required)
%   X                  Independent variable.
%                      (:,1) double
%
%   Y                  Dependent variable, same length as |X|.
%                      (:,1) double
%
%   (Optional)
%   W                  Weight for each data point, same length as |X|. A
%                      valid input for these weights would be, for example,
%                      the number of observations for each data point
%                      (default: equal weights for all data points).
%                      (:,1) double
%
% Output Arguments:
%   c                  Regression coefficients for Y = c(1) + c(2)*X.
%                      (2,1) double
%
%   adjR2              Adjusted R² value for the linear regression.
%                      (1,1) double
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% To calculate many regressions in a loop that (for some of them) might
% result in badly scaled, (nearly) singular or rank-deficient matrices, the
% loop can be enclosed with
% warning('off')
% warning('on')
% to suppress warnings for enhanced speed.
%
% The calculation of the coefficient of determination was coded along the
% equations described in the publication:
%
% Willet, J.B. and Singer, J.D. Another Cautionary Note about R²: Its Use
% in Weighted Least-Squares Regression Analysis. Am Stat 42, 236-238
% (1988).
% https://doi.org/10.2307/2685031
%
% In case of the weighted linear regression, the authors argue that is
% more appropriate not to include the weights in the calculation of R² to
% better relate them to the original variables.
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
% Initial release: 2021-12-30
% Last revision: 2021-12-30

%% Function argument validation
arguments
    X (:,1) double
    Y (:,1) double
    W (:,1) double = ones(size(X))
end

%% Main

% Add a constant term to the data.
n = length(X);
X = [ones(n,1), X];

% Include the weights.
Xstar = sqrt(W) .* X;
Ystar = sqrt(W) .* Y;

% Estimate the coefficients.
c = Xstar \ Ystar;

% Compute the coefficient of determination (R²). For ordinary least squares
% (OLS, equal weights), this follows eq.(1) in Willet and Singer (1988).
% Their modification for weighted least squares (WLS) in eq.(7) is
% identical in this context, as the weights are reflected in the (already
% estimated) OLS or WLS coefficients.
R2 = 1 - (Y - X*c).' * (Y - X*c) / (Y.' * Y - n*mean(Y)^2);

% R² increases with more explanatory variables in the model, i.e., the more
% complex the model, the higher the (apparent) R² value. The adjusted
% coefficient of determination counters this effect by penalizing the
% addition of variables (beyond the constant term). For linear regressions,
% this amounts to (n - number of variables - 1) = (n - 2).
adjR2 = 1 - (1 - R2) * (n - 1) / (n - 2);

end
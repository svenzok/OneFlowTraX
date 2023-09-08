function M = matchpairsSlim(Cost, costUnmatched)
% matchpairsSlim is a trimmed version of MATLAB's >matchpairs< function.
%
% Syntax:
%   M = matchpairsSlim(Cost, costUnmatched)
%
% Input Arguments:
%   (Required)
%   Cost               Cost matrix. Each entry Cost(i,j) specifies the cost
%                      of assigning row i to column j.
%                      (:,:) double
%
%   costUnmatched      Cost of not matching, specified as a scalar.
%                      >matchpairsSlim< compares the value of
%                      2*|costUnmatched| to the entries in |Cost| to
%                      determine whether it is more beneficial for a row or
%                      column to remain unmatched. Use this parameter to
%                      make matches more or less likely in the algorithm.
%                      (1,1) double
%
% Output Arguments:
%   M                  Nx2 assigment list containing the indices of
%                      matching pairs.
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Adapted from MATLAB's >matchpairs< function, trimmed for enhanced speed.
%
% From >matchpairs<: "Solves the linear assignment problem for the rows and
% columns of the |Cost| matrix. Each row is assigned to a column in such a
% way that the global cost is minimized. The scalar |costUnmatched|
% specifies the cost of not assigning a row to a column, or a column to a
% row, at all".
%
% This assignment algorithm solves the same problem as the Hungarian
% algorithm, but has a different complexity. While the Hungarian algorithm
% has complexity O(N^4) and a further developed form has O(N^3), the
% algorithm in >matchpairs< has O(N^3*log(N)). N is the number source
% points.
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
% Initial release: 2021-08-25
% Last revision: 2022-11-16

%% Function argument validation
arguments
    Cost (:,:) double
    costUnmatched (1,1) double
end

%% Main

% Create a larger cost matrix, with dummy rows and columns to account for
% the possibility of not matching.
[m, n] = size(Cost);

paddedCost = Inf(m+n, m+n,'like',Cost);
paddedCost(1:m, 1:n) = Cost;
paddedCost(m+1:end, n+1:end) = Cost.';
for ii=1:m
    paddedCost(ii, n+ii) = 2*costUnmatched;
end
for jj=1:n
    paddedCost(m+jj, jj) = 2*costUnmatched;
end
colToRow = matlab.internal.graph.perfectMatching(paddedCost);

% Real row matched to real column: a matched pair.
matchedPairs = colToRow(1:n) <= m;
M = [colToRow(matchedPairs), find(matchedPairs)];
M = reshape(M, [], 2); % Correct size in empty case.

end
function changedParameter = compareS(oldS, newS)
% compareS finds the changed parameter between two parameter structures.
%
% Syntax:
%   changedParameter = compareS(oldS, newS)
%
% Input Arguments:
%   (Required)
%   oldS               Old (or original) parameter structure (see below).
%                      (1,1) struct
%
%   newS               New (or updated) parameter structure (see below).
%                      (1,1) struct
%
%
% Output Arguments:
%   changedParameter   Changed parameter, with the name of the parameter in the first cell and its
%                      new value in the second cell.
%                      (1,2) cell
%
% Other required m-files: none
% Subfunctions: unpackS, unpackTable
% Additional required MATLAB products: none
%
% Notes:
% This program is intended for the parameter structure S that is created for the OneFlowTraX app, to
% facilitate the update of parameters that the user either changes in the app itself or in the
% accompanying parameter app. Accordingly, the scope of this program is very narrow. First, the
% fields of the structures |oldS| and |newS| are itself structures. Each of these structures can
% then either feature a parameter as a fieldname-value pair or fieldname-table pair (the latter will
% be separated into individual name-value pairs). It is also possible to have another substructure,
% but this substructure must only have fieldname-value pairs (no tables). This allows for some
% flexibility, should the parameter structure, order or names be changed. The output can then be
% used to find the correct GUI element in the OneFlowTraX app (by tagging it with a corresponding name).
% This intermediate step maintains control over the GUI elements within an app. However, this still
% requires that all related GUI elements be tagged with the correct strings. Changes to the
% parameter structure (e.g. their names) must be adjusted manually in the OneFlowTraX app. The detailed
% name in |changedParameter| should help.
%
% Examples for resulting names:
% S.tracking.maxGapClose:
% 'tracking maxGapClose'
%
% S.localization.limits (this is a table), changing the 'upper limit' value for 'photons':
% 'localization limits photons (upper limit)'
%
% S.cluster.Voronoi.trimBoundary (as an example for a substructure):
% 'cluster Voronoi trimBoundary'
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
% Initial release: 2023-06-06
% Last revision: 2022-06-06

%% Function argument validation
arguments
    oldS (1,1) struct
    newS (1,1) struct
end

%% Main

oldC = unpackS(oldS);
newC = unpackS(newS);
changedParameter = newC(cellfun(@(x,y) ~isequal(x,y),oldC(:,2), newC(:,2)),:);

end

%% Subfunctions
function allC = unpackS(S)
% Unpacks the parameter structure into a cell array that only has one name for each parameter (the
% parameter value can, however, have different formats like strings, arrays etc.).
fn = fieldnames(S);
C = struct2cell(S);
allC = cell(numel(C),1);

for c = 1:numel(C)
        subFn = cellstr(string(fn{c}) + " " + string(fieldnames(C{c})));
        subC = struct2cell(C{c});
        for s = 1:numel(subC)
            if isstruct(subC{s})
                sub2Fn = fieldnames(subC{s});
                sub2C = struct2cell(subC{s});

                for s2 = 1:numel(sub2C)     
                        sub2C(s2,2) = sub2C(s2);
                        sub2C(s2,1) = cellstr([subFn{s}, ' ', sub2Fn{s2}]);
                end
                subC{s,1} = sub2C(:,1);
                subC{s,2} = sub2C(:,2);

            elseif istable(subC{s})
                T = unpackTable(subC{s});
                subC{s,1} = cellstr(subFn{s} + " " + string(T(:,1)));
                subC{s,2} = T(:,2);
            else
                subC(s,2) = {subC(s)};
                subC(s,1) = {subFn(s)};
            end
        end
        allC{c} = [vertcat(subC{:,1}), vertcat(subC{:,2})];
end

allC = vertcat(allC{:});

end

function T = unpackTable(T)
% Tables will be unpacked so that the row and variable names form new, unique names for each value.
rowNames = T.Properties.RowNames;
columnNames = T.Properties.VariableNames;
fullNames = string(rowNames) + " (" + string(columnNames) + ")";
fullNames = reshape(fullNames',[],1);
allValues = reshape(table2array(T)',[],1);
T = [cellstr(fullNames), num2cell(allValues)];
end
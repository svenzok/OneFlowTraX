function S = getDefaultS
% getDefaultS sets up the parameter structure S and fills it with some default values (OneFlowTrax app).
%
% Syntax:
%   S = getDefaultS
%
% Output Arguments:
%   S                  Structure holding parameter values for localization, tracking, cluster analysis etc.
%                      (1,1) struct
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
% Initial release: 2023-06-21
% Last revision: 2023-06-21

%% Main

% Sets up the parameter structure S.
S = struct;

% Steps that will be performed in a batch analysis.
T = table('Size', [4 1], 'VariableTypes',{'logical'},'VariableNames',{'doAnalysis'},...
    'RowNames',{'localization','tracking','MSD analysis','cluster analysis'});
T{:,:} = true(4,1);
S.general.analysisSteps = T;
S.general.batchText = "";

% Checks if raw data of the analysis steps (localization lists, track arrays) should be given as an additional output,
% and also which suffixes these files should have.
T = table('Size', [2 2], 'VariableTypes',{'logical','string'},...
    'VariableNames',{'raw data output','file suffix'},'RowNames',{'localization','tracking'});
T(:,:) = {false, "_locs"; false, "_tracks"};
S.general.rawDataOutput = T;

% Input file lists and name of the output file.
S.files.inputFiles = cell(0,1);
S.files.noiseFiles = cell(0,1);
S.files.maskFiles = cell(0,1);
S.files.outputFile = '';

% Localization parameters
S.localization.varianceMap = zeros(128, 'single');
S.localization.offset = 100;
S.localization.conversion = 0.48;
S.localization.pixelSize = 100;
S.localization.filterSize = 1.3;
S.localization.cutoff = 3;
S.localization.PSFROI = 7;
T = table('Size', [3 2], 'VariableTypes',{'double','double'}, ...
    'VariableNames',{'lower limit','upper limit'},'RowNames',{'PSF','photons','precision'});
T(:,:) = {0 Inf;0 Inf;0 Inf};
S.localization.limits = T;

% Tracking parameters.
S.tracking.maxLinkDistance = 150;
S.tracking.maxGapClose = 4;
T = table('Size', [4 2], 'VariableTypes',{'double','double'}, ...
    'VariableNames',{'lower limit','upper limit'},'RowNames',{'N', 'duration','displacement','range'});
T(:,:) = {3 Inf;0 Inf;0 Inf; 0 Inf};
S.tracking.limits = T;
S.tracking.maskPolygons = cell(1,0);

% MSD analysis parameters.
S.MSDanalysis.frameRate = 50;
S.MSDanalysis.minTrackLength = 8;
S.MSDanalysis.adjR2 = 0.4;
S.MSDanalysis.linFitPoints = [2 4];
S.MSDanalysis.splitValue = -1.75;

% Cluster analysis parameters.
S.cluster.algorithm = 'Voronoi';
S.cluster.input = 'track centroids';
S.cluster.minN = 5;
S.cluster.Voronoi.method = 'density';
S.cluster.Voronoi.alpha = 2;
S.cluster.Voronoi.areaThreshold = 20;
S.cluster.Voronoi.trimEdges = true;
S.cluster.Voronoi.trimBoundary = 'AABB';
S.cluster.Voronoi.trimMethod = 'intersect';
S.cluster.DBSCAN.minPts = 4;
S.cluster.DBSCAN.epsilon = 280;
S.cluster.NASTIC.r = 1;
S.cluster.segNASTIC.r = 1;
S.cluster.segNASTIC.nSegmentOverlaps = 2;

end


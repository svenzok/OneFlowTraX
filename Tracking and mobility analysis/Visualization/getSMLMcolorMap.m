function cmap = getSMLMcolorMap(colorMap)
% getSMLMcolorMap creates colormaps, especially for plotting particle tracks as surface objects.
%
% Syntax:
%   cmap = getSMLMcolorMap(colorMap)
%
% Input Arguments:
%   (Required)
%   colorMap           Specifies the track colors (see below for all options). Some choices are
%                      intended to work with certain visualization options:
%                      'cycle X colors' defines a colormap with X colors, which the tracks can then
%                                       be randomly assigned to.
%                      'X/Y' defines two colors X and Y as a dual-color option to paint tracks with
%                            X when they are under, and with Y when they are over a certain
%                            threshold.
%                      'X/Y/Z' defines three colors X,Y and Z as a triple-color option. 
%                      The rest uses continuous color maps or just paints all tracks black or white.
%                      (1,:) char | (1,1) string
%                      {'cycle 4 colors', 'cycle 6 colors', 'cycle 8 colors', 'cycle 10 colors', ...
%                       'orange/blue', 'blue/orange', 'red/teal', 'teal/red', 
%                       'gray/orange/blue', 'gray/blue/orange', 'gray/red/teal', 'gray/teal/red'...
%                       'inferno', 'inferno (reversed)', 'plasma', 'plasma (reversed)', 'hot', ...
%                       'hotInv', 'black', 'white'}
%
% Output Arguments:
%   cmap               Colormap that can be used when plotting particle trajectories as a surface
%                      object (the values in the CData property then determine which color is
%                      chosen from |cmap| for each vertex).
%                      (:,3) double
%
% Other required m-files: getColors.m
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Developed to work in tandem with >prepareTracksPlotData<.
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
% Initial release: 2023-03-21
% Last revision: 2023-03-30

%% Function argument validation
arguments   
    colorMap (1,:) char {mustBeMember(colorMap, ...
        {'cycle 4 colors','cycle 6 colors','cycle 8 colors', 'cycle 10 colors',...
        'orange/blue','blue/orange','red/teal','teal/red',...
        'gray/orange/blue','gray/blue/orange','gray/red/teal','gray/teal/red',...
        'inferno','inferno (reversed)','plasma','plasma (reversed)','hot','hot (reversed)',...
        'black','white'})};
end

%% Main

% Get color-blind friendly colors and colormaps.
Colors = getColors;
colorNames = [Colors.CB.Properties.RowNames; 'gray'];
colorTriplets = [Colors.CB.all;[0.5 0.5 0.5]];

% Special treatment for cycling colors and dual-/triple-color.
if contains(colorMap,'cycle')
    nCycleCols = str2double(extract(colorMap, digitsPattern));
    colorMap = 'cycle';
elseif contains(colorMap,'/')
    multiColorNames = strsplit(colorMap,'/');
    [~,~,multiID] = intersect(multiColorNames, colorNames, 'stable');
    colorMap = 'multi-color';
else
    % empty
end

switch colorMap
    case 'cycle'
        % Use a set of color-blind friendly colors.
        cmap = colorTriplets(1:nCycleCols,:);
    case 'multi-color'
        % Mulitple colors according to the input string.
        cmap = colorTriplets(multiID,:);
    case 'black'
        cmap = zeros(256,3);
    case 'white'
        cmap = ones(256,3);
    otherwise
        % Use one of the stored color maps.
        cmap = Colors.cmap{colorMap,:}{:};
end

end
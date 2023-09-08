function [rectanglePos, labelStr, labelPos] = placeScaleBar(axesLimits, unitExp)
% placeScaleBar gets the position and label for a neatly formatted scale bar.
%
% Syntax:
%   [rectanglePos, labelStr, labelPos] = placeScaleBar(axesLimits, unitExp)
%
% Input Arguments:
%   (Required)
%   axeslimits         Axes limits in the format [xmin, xmax, ymin, ymax].
%                      (1,4) double
%
%   unitExp            Length unit exponent of |axesLimits| if expressed in meters (e.g. -9 for
%                      nanometers). An algorithm will choose an appropriate unit for the scale bar.
%                      (1,1) double
%
% Output Arguments:
%   |rectanglePos|     Vector [x y w h] for the 'Position' property of >rectangle<, with x and y
%                      as coordinates for the upper left corner, and w and h as width and height
%                      (all in data units and with the y-axis reversed, i.e., the axis origin is in
%                      the upper left corner).
%                      (1,4) double
%
%   |labelStr|         String containing the scale bar text.
%                      (1,1) string
%
%   |labelPos|         Vector denoting the x and y coordinates of the scale bar text in data units,
%                      provided that the horizontal alignment of the text will be formatted as
%                      'center', and the vertical alignment as 'baseline'.
%                      (1,2) double
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
%
% The following design principles were chosen:
% The bottom right corner of the rectangle should lie on the imaginary diagonal line of the image,
% at about 7% of its length (measured from the bottom right). The bar height should equal 1.5% of
% the image height, and the bar width 10% to 20% (preferably 20%) of the image width. This interval
% allows some flexibility for a neat formatting of the label text (see below). The baseline of the
% label text should be one bar height above the scale bar. The font size must be set (outside this
% function) using a font size equal to ca. 8% of the image height by first setting 'FontUnits' to
% 'normalized' and then setting 'FontSize' to 0.08.
% Be aware that "images" mentioned here are meant to have real units (not pixel coordinates). The
% scale bar is therefore best placed into an image that has been converted into a surface object to
% reflect units that are not dependent on the pixel size.
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
% Initial release: 2023-03-24
% Last revision: 2023-03-24

%% Function argument validation
arguments
    axesLimits (1,4) double
    unitExp (1,1) double
end

%% Main

% Image width and height.
wIm = axesLimits(2) - axesLimits(1);
hIm = axesLimits(4) - axesLimits(3);

% Neat data unit lengths:
neatLengths = [1 2 5 10 20 50 100 200 500]';

% Collection of bar widths within 10% to 30% of the image width in meters.
barWidthList = (0.1:0.005:0.30)' * wIm * 10^unitExp;

% Divide by exponential factors (engineering notation), covering kilometers to picometers.
unitExpList = 3:-3:-12;
unitStrList = {'km','m','mm','µm','nm','pm'};
barWidthList = barWidthList./10.^unitExpList;

% The array |barWidthList| now contains a list of numbers that could be plotted over the length
% limited scale bars for each unit ('km', 'm' ... as columns) To find the one that is closest to one
% of the neat data unit lengths, divide (3rd dimension) by them to get fractional values.
barWidthList = barWidthList./reshape(neatLengths,1,1,[]);

% Get the linear indices of the three values closest to 1.
[~, IDs] = mink(barWidthList(:) - 1, 3, 'ComparisonMethod', 'abs');

% Find their respective units and neat lengths.
[~, unitID, neatID] = ind2sub(size(barWidthList), IDs);
unitExpFinal = unitExpList(unitID);
barWidthFinal = neatLengths(neatID);

% Calculate their width in the image and choose the one that is closest to 20% of the image width.
tmp = barWidthFinal'.*10.^unitExpFinal / 10^unitExp / wIm;
[~, ID] = min(tmp - 0.2, [], 'ComparisonMethod', 'abs');

% Prepare the label text output.
labelStr = string(barWidthFinal(ID)) + " " + unitStrList{unitID(ID)};

% Determine the width and height of the bar rectangle.
w = barWidthFinal(ID) * 10^unitExpFinal(ID) / 10^unitExp;
h = 0.015 * hIm;

% Calculate the lower left corner coordinates of the rectangle.
x = axesLimits(2) - 0.07*wIm - w;
y = axesLimits(4) - 0.07*hIm - h;

% The horizonally centered text position should be centered above the bar, with one bar height
% between the text's baseline and the bar.
labelX = x + 0.5 * w;
labelY = y - h;

% Combine for the output.
rectanglePos = [x y w h];
labelPos = [labelX labelY];

end
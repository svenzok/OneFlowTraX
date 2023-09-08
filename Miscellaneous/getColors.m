function ColorCollection = getColors
% getColors provides the user with a selection of colors and colormaps, including color blind friendly options, colors
%           according to the light wavelength and a collection of design colors (all as RGB triplets).
%
% Syntax:
%   ColorCollection = getColors
%
% Output Arguments:
%   ColorCollection    Structure with fields
%                      .CB: ten color blind friendly colors as a table with row names:
%                           'red','green','yellow','blue','orange','cyan', 'magenta','teal',
%                           'lavender','maroon'
%                      .lambda: matrix, each row holds the RGB triplet for the wavelength specified
%                               by the row number in nanometers
%                      .UT: design colors, as a table with row names:
%                           'carmine', 'gold', 'charcoal', 'darkBlue', 'blue', 'lightBlue', 
%                           'paleGreen', 'green', 'darkGreen', 'rust', 'purple', 'taupe', 'straw', 
%                           'curry', 'mud'
%                      .cmap: colormaps and their reversed counterparts, as a table with row names:
%                             'inferno', 'inferno (reversed)', 'plasma', 'plasma (reversed)', 
%                             'hot', 'hot (reversed)'
%                       
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% In the following, some examples how colors (as RGB triplets) and colormaps (as 256x3 arrays) can
% be accessed from the structure C:
%
% Example 1: teal color from the CB palette (all the CB format examples also apply to the UT palette)
% C.CB{'teal',:}
%
% Example 2: the first three colors from the CB palette
% C.CB{1:3,:}
%
% Example 3: twenty random CB colors (duplicates allowed)
% C.CB{randi(height(C.CB),20,1),:}
%
% Example 4: three different CB colors (no duplicates allowed)
% C.CB{randperm(height(C.CB),3),:}
%
% Example 5: all colors from the CB palette
% C.CB.all
%
% Example 6: RGB color that matches a light wavelength of 500 nanometers
% C.lambda(500,:)
%
% Example 7: visible colors (light wavelengths from 400-780 nanometers)
% C.lambda(400:780,:)
%
% Example 8: the reversed plasma colormap
% C.cmap{'plasma (reversed)',:}{:}
%
% Example 9: all colormaps, each in one cell of a cell array
% C.cmap.all
%
% Example 10: the names of all colormaps as a cell array of character arrays (this also works for
% the color names in .CB and .UT).
% C.cmap.Properties.RowNames
%
% The color blind friendly colors are taken from:
% Trubetskoy, S. (2017, January 11). List of 20 Simple, Distinct Colors.
% Retrieved from https://sashamaps.net/docs/resources/20-colors/
%
% The wavelength to RGB conversion was adapted from FORTRAN code by:
% Bruton, D. (1996, February 20). Approximate RGB values for Visible Wavelengths.
% Retrieved from http://www.physics.sfasu.edu/astro/color/spectra.html
%
% The colormaps 'inferno' and 'plasma' were adapted from:
% Hunter, J.D. Matplotlib: A 2D Graphics Environment. Comput Sci Eng 9(3), 90-95 (2007).
% https://doi.org/10.1109/MCSE.2007.55.
%
% Tested: MATLAB Version 9.11.0.1769968 (R2021b),
%         Microsoft Windows 10 Pro Version 10.0 (Build 19042)
%
% Author: Sven zur Oven-Krockhaus
%         Institute of Physical and Theoretical Chemistry
%         University of Tuebingen, Tuebingen, Germany
% E-mail: sven.zur-oven-krockhaus@uni-tuebingen.de
%
% GNU placeholder
%
% Initial release: 2021-12-10
% Last revision: 2023-03-30

%% Color blind friendly colors

namesCB = {'red','green','yellow','blue','orange','cyan','magenta', 'teal','lavender','maroon'};

RGB = [230,  25,  75; ...
        60, 180,  75; ...
       255, 225,  25; ...
         0, 130, 200; ...
       245, 130,  48; ...
        70, 240, 240; ...
       240,  50, 230; ...
         0, 128, 128; ...
       220, 190, 255; ...
       128,   0,   0] / 255;

ColorCollection.CB = table(RGB,'VariableNames',{'all'},'RowNames',namesCB);

%% Colors corresponding to light wavelength

% RGB color transitions for visible wavelength sections.
RGB1 = zeros(60,3);
wavelength = (380:439)';
RGB1(:,1) = -(wavelength-440)/(440-380);
RGB1(:,3) = 1;

RGB2 = zeros(50,3);
wavelength = (440:489)';
RGB2(:,2) = (wavelength-440)/(490-440);
RGB2(:,3) = 1;

RGB3 = zeros(20,3);
wavelength = (490:509)';
RGB3(:,2) = 1;
RGB3(:,3) = -(wavelength-510)/(510-490);

RGB4 = zeros(70,3);
wavelength = (510:579)';
RGB4(:,1) = (wavelength-510)/(580-510);
RGB4(:,2) = 1;

RGB5 = zeros(65,3);
wavelength = (580:644)';
RGB5(:,1) = 1;
RGB5(:,2) = -(wavelength-645)/(645-580);

RGB6 = zeros(136,3);
RGB6(:,1) = 1;

% Join all sections.
RGB = vertcat(RGB1, RGB2, RGB3, RGB4, RGB5, RGB6);

% Create intensity factors to mimick human perception.
wavelength = (380:419)';
factor1 = 0.3+0.7*(wavelength-380)/(420-380);

factor2 = ones(281,1);

wavelength = (701:780)';
factor3 = 0.3+0.7*(780-wavelength)/(780-700);

factor = vertcat(factor1, factor2, factor3);

% Adjust the perceived intensity of the colors.
RGB = RGB .* factor;

% Fill into an array so that the row number equals the wavelength.
ColorCollection.lambda = zeros(800,3);
ColorCollection.lambda(380:780,:) = RGB;

%% University of Tuebingen colors

namesUT = {'carmine','gold','charcoal','darkBlue','blue','lightBlue','paleGreen','green', ...
    'darkGreen','rust','purple','taupe','straw','curry','mud'};

% Each row corresponds to one color in the above order.
RGB = [165,  30,  55; ...
       180, 160, 105; ...
        50,  65,  75; ...
        65,  90, 140; ...
         0, 105, 170; ...
        80, 170, 200; ...
       130, 185, 160; ...
       125, 165,  75; ...
        50, 110,  30; ...
       200,  80,  60; ...
       175, 110, 150; ...
       180, 160, 150; ...
       215, 180, 105; ...
       210, 150,   0; ...
       145, 105,  70] / 255;

ColorCollection.UT = table(RGB,'VariableNames',{'all'},'RowNames',namesUT);

%% Colormaps (inferno and plasma colormaps are color blind friendly)
cmInferno = [  0,   0,   4; ...
              40,  11,  84; ...
             101,  21, 110; ...
             159,  42,  99; ...
             212,  72,  66; ...
             245, 125,  21; ...
             250, 193,  39; ...
             252, 255, 164] / 255;
cmInferno = interp1(linspace(1,256,8), cmInferno, 1:256);
cmInfernoRev = flipud(cmInferno);

cmPlasma = [ 13,   8, 135; ...
             84,   2, 163; ...
            139,  10, 165; ...
            185,  50, 137; ...
            219,  92, 104; ...
            244, 136,  73; ...
            254, 188,  43; ...
            240, 249,  33] / 255;
cmPlasma = interp1(linspace(1,256,8), cmPlasma, 1:256);
cmPlasmaRev = flipud(cmPlasma);

cmHot = hot(256);
cmHotRev = flipud(cmHot);

namesCmap = {'inferno', 'inferno (reversed)', 'plasma', 'plasma (reversed)', ...
    'hot', 'hot (reversed)'};

ColorCollection.cmap = table({cmInferno;cmInfernoRev;cmPlasma;cmPlasmaRev;cmHot;cmHotRev}, ...
    'VariableNames',{'all'},'RowNames',namesCmap);

end

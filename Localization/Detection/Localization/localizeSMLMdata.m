function [locList, outlierInfo] = localizeSMLMdata(imageStack, Options)
% localizeSMLMdata performs single-molecule localization analysis for image
%                  data stored as a 3D array (image stack).
%
% Syntax:
%   [locList, outlierInfo] = localizeSMLMdata(imageStack, Options)
%
% Input Arguments:
%   (Required)
%   imageStack         Image stack as a 3D array, one page for each frame
%                      (e.g., of a multipage TIFF file).
%                      (:,:,:) (uint8, uint16, single, double)
%  
% Name-Value Pair Input Arguments (Options):
%   (Optional)
%   varianceMap        Variance map to account for pixel-dependent noise of
%                      sCMOS cameras, must have the same size as each frame
%                      of |imageStack| (default: all zero).
%                      (:,:) double
%
%   offset             Offset ADU value of the camera (default: 100).
%                      (1,1) double
%
%   conversion         Conversion factor (electrons per ADU) of the camera
%                      (default: 0.48).
%                      (1,1) double
%
%   pixelSize          Pixel size in nm, indicating the length of one pixel
%                      in the sample plane (default: 100)
%                      (1,1) double
%
%   filterSize         Spatial band-pass filter sigma value for the initial
%                      peak finding in pixels (default: 1.3).
%                      (1,1) double
%
%   cutoff             Cutoff value in photons to determine which initial
%                      peaks are further processed, i.e., peaks with photon
%                      numbers below this value are ignored (default: 3)
%                      (1,1) double
%   
%   PSFROI             Width of the square area to be cut out in pixels
%                      around each potential PSF (each found peak maximum).
%                      As the center pixel of this area is placed on the
%                      pixel with the peak maximum, the width must be an
%                      odd value (default: 7).
%                      (1,1) double
%
%   useGPU             Logical switch to use the GPU if possible (default: true).
%                      (1,1) logical
%
% Output Arguments:
%   locList            N x 7 array containing the localization data, with
%                      columns:
%                      1 frame (the image frame of the localization)
%                      2 x (x-coordinate in nm)
%                      3 y (y-coordinate in nm)
%                      4 sdv (sigma value of the symmetric Gaussian the PSF was fitted with, in nm)
%                      5 photons (detected photons, estimated by using the sCMOS variance map, offset, conversion
%                                 factor and the MLE algorithm)
%                      6 background (background in photons)
%                      7 localization precision (localization precision according to the Cramer-Rao lower bound, in nm)
%                      (:,7) double
%
%   outlierInfo        Info about ignored localizations due to fitting
%                      errors or irregularities as a 1x5 vector, with
%                      columns:
%                      1 total number of original localizations
%                      2 localizations with any fitting parameter being
%                        either 0, NaN or inf
%                      3 obvious outliers (see code)
%                      4 fitting errors (the algorithm got stuck at a hard
%                        limit, visible as a rounded value in the results)
%                      5 remaining localizations after the errors above are
%                        sorted out (some might exhibit several of these
%                        errors, so the number of remaining localizations
%                        might be higher than the number of original
%                        localizations minus errors)
%                      (1,5) double
%  
% Other required files: GPUmleFit_LM.mexw64, CPUmleFit_LM.mexw64, cudart64_80.dll, 
%                       extractPeakSections.m, find3x3peaksCall.m,
%                       find3x3peaks.mexw64
% Subfunctions: createGaussian2Dfilter
% Additional required MATLAB products: none
%
% Notes:
% Makes extensive use of the localization routine published in:
%
% Ries, J. SMAP: a modular super-resolution microscopy analysis platform
% for SMLM data. Nat Methods 17, 870–872 (2020).
% https://doi.org/10.1038/s41592-020-0938-1.
%
% Their code was rewritten for better readability and simplified to be most
% suitable for a sCMOS camera (default values are based on an ORCA-Flash4.0
% V2 Digital CMOS camera C11440-22CU from Hamamatsu, Japan). This
% adaptation aims at better integration and faster execution, is therefore
% limited to 2D imaging with symmetric PSFs and also disregards other
% advanced features of the original program.
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
% Initial release: 2021-12-31
% Last revision: 2023-08-24

%% Function argument validation
arguments
    imageStack (:,:,:) {mustBeA(imageStack, ...
        ["uint8", "uint16", "single", "double"])}
    Options.varianceMap (:,:) double = zeros(size(imageStack,[1 2]));
    Options.offset (1,1) double = 100
    Options.conversion (1,1) double = 0.48
    Options.pixelSize (1,1) double = 100
    Options.filterSize (1,1) double = 1.3
    Options.cutoff (1,1) double = 3
    Options.PSFROI (1,1) double = 7
    Options.useGPU (1,1) logical = true
end

%% Localization

% Create a difference of Gaussians (DoG) filter, corresponding to the
% chosen filter size.
ROIsize = max(ceil(6*Options.filterSize + 1),3);
DoG = createGaussian2Dfilter(ROIsize, Options.filterSize) - ...
      createGaussian2Dfilter(ROIsize, max(1, 2.5*Options.filterSize));

% Convert the pixel intensites to photons.
imageStackPhotons = (single(imageStack)-Options.offset)*Options.conversion;
Options.varianceMap = Options.varianceMap*Options.conversion^2;

% smallest photon number in each frame
minPhotons = min(imageStackPhotons,[],[1 2]);

% Filter all images with the DoG.
imageStackFiltered = convn(imageStackPhotons - minPhotons, DoG, 'same');

% Find local 2D maxima for all frames (peaks).
peakList = find3x3peaksCall(imageStackFiltered, Options.cutoff, Options.PSFROI);

% Cut out each peak. 
[trimmedPeakList, croppedPeaksStack, croppedVarStack] = ...
    extractPeakSections(peakList, Options.PSFROI, imageStackPhotons, ...
    Options.varianceMap);

% Fit the peak stack on the GPU or CPU.
if Options.useGPU
    [Pcspline, CRLB] = GPUmleFit_LM(croppedPeaksStack,2,50,1.5,croppedVarStack,1,0);
else
    [Pcspline, CRLB] = CPUmleFit_LM(croppedPeaksStack,2,50,1.5,croppedVarStack,1,0);
end

% Put the results into the format: frame, x-coordinates (nm),
% y-coordinates (nm), PSF size (nm), photons, photons background,
% localization precision (CRLB). 
frame = trimmedPeakList(:,3);
dROI = floor(Options.PSFROI/2);
x_nm = (Pcspline(:,2) - dROI + trimmedPeakList(:,1))*Options.pixelSize;
y_nm = (Pcspline(:,1) - dROI + trimmedPeakList(:,2))*Options.pixelSize;
PSFsize = Pcspline(:,5)*Options.pixelSize;
photons = Pcspline(:,3);
photonsBackground = Pcspline(:,4);

% Calculate the localization precision by averaging over the x and y direction.
precisionCRLB = Options.pixelSize*mean(real(sqrt(CRLB(:,[2 1]))),2);
       
locList = [frame, x_nm, y_nm, PSFsize, photons, photonsBackground, precisionCRLB];

%% Determination of outliers

outlierInfo = zeros(1,5);
    
% original number of localizations
outlierInfo(1) = size(locList, 1);

% rows with <=0, NaN or inf
isOutlier1 = any(isnan(locList) | locList<=0 | isinf(locList), 2);
outlierInfo(2) = nnz(isOutlier1);

% Assume that values greater than 10 times the data median can be
% considered outliers (empirical estimation, tested for PSF size, photons,
% localization precision). There is a negligible number of data points
% already at five times the data median, but choose ten times to be sure.
% The photon distribution is heavily tailed, but usually features no
% apparent outliers due to fitting errors. Here, the cut-off can be set at
% 50*median.
med = median(locList(:,[4 5 7]),1,'omitnan');
med = med.*[10 50 10];
isOutlier2 = any(locList(:,[4 5 7])>med, 2);
outlierInfo(3) = nnz(isOutlier2);

% rows with fitting errors
isOutlier3 = any(round(locList(:,2:end))./locList(:,2:end) == 1, 2);
outlierInfo(4) = nnz(isOutlier3);

% remove all outliers from the localization list 
locList(isOutlier1 | isOutlier2 | isOutlier3, :) = [];

% store the number of remaining localizations
outlierInfo(5) = size(locList, 1);

end

%% Subfunctions
function lowpassFilter = createGaussian2Dfilter(ROIsize, filterSize)
% Creates a rotationally symmetric Gaussian lowpass filter of size
% |ROIsize| with standard deviation |filterSize|. This is the same as the
% MATLAB function >fspecial< with the gaussian option.

% create grid of x and y values
sz = ([ROIsize, ROIsize]-1)/2;
[x,y] = meshgrid(-sz(2):sz(2),-sz(1):sz(1));

% function for the filter
lowpassFilter = exp(-(x.*x + y.*y)/(2*filterSize^2));

% set very small values to zero
lowpassFilter(lowpassFilter<eps*max(lowpassFilter(:))) = 0;

% normalize the filter 
sumFilter = sum(lowpassFilter(:));
if sumFilter ~= 0
    lowpassFilter  = lowpassFilter/sumFilter;
end

end
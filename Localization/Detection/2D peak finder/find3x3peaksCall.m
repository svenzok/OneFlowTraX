function peakList = find3x3peaksCall(imageStack, minIntensity, ROIwidth)
% find3x3peaksCall searches for local maxima in a 3x3 neighborhood.
%
% Syntax:
%   peakList = find3x3peaksCall(imageStack, minIntensity, ROIwidth)
%
% Input Arguments:
%   (Required)
%   imageStack         Image stack as a 3D array, one page for each frame
%                      (but also works for a single image).
%                      (:,:,:) single
%
%   minIntensity       Minimal intensity for a pixel to be processed.
%                      (1,1) double
%
%   ROIwidth           Width of the square ROI in pixels that will be
%                      centered on each peak to cut out a respective image
%                      area (later on). This is already integrated in this
%                      function to prevent searching for maxima in the
%                      border region of the image that could be cut out
%                      with the specified ROI width.
%                      (1,1) double
%
% Output Arguments:
%   peakList           List of peak coordinates in the format
%                      [x-coordinate, y-coordinate, frame].
%                      The coordinates are given in pixels.
%                      (:,3) double
%
% Other required files: find3x3peaks.mexw64
% Subfunctions: none
% Additional required MATLAB products: none
%
% Notes:
% Wrapper for find3x3peaks.mexw64, which was coded along the procedure
% described in section 3 (scan-line algorithm for 3x3 NMS) of the
% publication:
% Pham, T.Q. Non-maximum Suppression Using Fewer than Two Comparisons per
% Pixel. In: Blanc-Talon J. et al. (eds). Advanced Concepts for Intelligent
% Vision Systems. ACIVS 2010. Lecture Notes in Computer Science, vol 6474.
% Springer, Berlin, Heidelberg.
% https://doi.org/10.1007/978-3-642-17688-3_41
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
% Initial release: 2022-01-07
% Last revision: 2022-12-19

%% Function argument validation
arguments
    imageStack (:,:,:) single
    minIntensity (1,1) single
    ROIwidth (1,1) double
end

%% Main

% The mex function expects an image stack. To make it work for single
% images, this workaround just adds a second, empty page.
if ndims(imageStack) ~= 3
    imageStack = cat(3,imageStack,imageStack.*0);
end

% The mex function also requires the maximum number of potential peaks 
% for pre-allocation. The following value will be an overestimatation, the
% list will be trimmed later.
maxPeaks = nnz(imageStack >= minIntensity);

% call mex function
peakList = find3x3peaks(imageStack,minIntensity,maxPeaks,ROIwidth);

% trim the list
lastPeak = find(peakList == 0, 1);
if lastPeak ~= 0
    peakList = peakList(1:lastPeak-1,:);
else
    peakList = [];
end

end


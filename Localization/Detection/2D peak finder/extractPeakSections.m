function [trimmedPeakList, croppedPeaksStack, croppedVarStack] = ...
    extractPeakSections(peakList, ROIwidth, imageStack, varianceMap)
% extractPeakSections cuts out square areas of an image stack.
%
% Syntax:
%   [trimmedPeakList, croppedPeaksStack, croppedVarStack] = ...
%    extractPeakSections(peakList, ROIwidth, imageStack, varianceMap)
%
% Input Arguments:
%   (Required)
%   peakList           List of peak coordinates in the format
%                      [x-coordinate, y-coordinate, frame].
%                      The coordinates must be given in pixels.
%                      (:,3) double
%
%   ROIwidth           Width of the square ROI in pixels that will be
%                      centered on each peak to cut out a respective image
%                      area.
%                      (1,1) double
%
%   imageStack         Image stack as a 3D array, one page for each frame.
%                      (:,:,:) single
%
%   varianceMap        Pixel-wise variances (e.g., of the sensor of a
%                      sCMOS camera) with the same x-and y-dimensions as
%                      |imageStack|.
%                      (:,:) single
%
% Output Arguments:
%   trimmedPeakList    Same as the input |peakList|, but without those
%                      peaks that could not be cut out due to their
%                      closeness to the image borders.
%                      (:,3) double
%
%   croppedPeaksStack   Image stack as a 3D array, one page for each peak
%                       that was cut out from |imageStack|, in the same
%                       order as |trimmedPeakList|.
%                       (:,:,:)
%
%   croppedVarStack    Corresponding variance map cutouts for each page in
%                      |croppedPeaksStack|.
%                      (:,:,:)
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products: none
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
% Initial release: 2022-01-04
% Last revision: 2022-03-13

%% Function argument validation
arguments
    peakList (:,3) double
    ROIwidth (1,1) double
    imageStack (:,:,:) single
    varianceMap (:,:) single
end

%% Main

% Remove all peak locations that are too close to the image borders (they
% would not be able to be fully extracted due to the given cutout size.)
isValid = peakList(:,1) > floor(ROIwidth/2) & ...
          peakList(:,1) < size(imageStack,2) - floor(ROIwidth/2) & ...
          peakList(:,2) > floor(ROIwidth/2) & ...
          peakList(:,2) < size(imageStack,1) - floor(ROIwidth/2); %% -1/+1??
trimmedPeakList = peakList(isValid,:);

% Calculate the index ranges for each peak (with respect to the single
% image dimensions).
colID = (trimmedPeakList(:,1) - floor(ROIwidth/2)) * ...
    ones(1,ROIwidth) + (ROIwidth - 1) * linspace(0,1,ROIwidth);
rowID = (trimmedPeakList(:,2) - floor(ROIwidth/2)) * ...
    ones(1,ROIwidth) + (ROIwidth - 1) * linspace(0,1,ROIwidth);

% Replicate the ranges to cover the ROI in x and y, also multiply the
% column ranges with the height of the image to get the offset for the
% linear indices when addressing a column (- 1, because this starts with an
% offset of 0 for the first column).
x = repmat(rowID, 1, ROIwidth); 
y = size(imageStack,1) * (repelem(colID, 1, ROIwidth) - 1);

% Each frame in the image stack creates an additional offset, given by the
% total number of pixels for one frame (-1, the first frame has no offset).
z = (trimmedPeakList(:,3) - 1) * numel(imageStack(:,:,1));

% The linear indices for the image stack can now be calculated by a
% simple addition.
linearIDlist = (x + y + z)';

% The linear indices for the variance map do not need the frame offsets.
linearIDlistVarianceMap = (x + y)';

% Now, each column in |linearIDlist| holds the linear indices for a cutout
% area, corresponding to each respective row in |peakList|. For further
% processing, they can be concatenated.
linearIDlist = linearIDlist(:);

% Pre-allocate the array for the cropped peaks and fill by indexing into
% the original image stack.
croppedPeaksStack = zeros(ROIwidth,ROIwidth,length(trimmedPeakList),...
    'like',imageStack);
croppedPeaksStack(:) = imageStack(linearIDlist);

% Each cropped peak is accompanied by a cutout of the variance map of the
% same image region. Here, the linear indices only refer to the single
% image of |varianceMap|, from which all required regions will be copied.
croppedVarStack = zeros(size(croppedPeaksStack),'like',croppedPeaksStack);
croppedVarStack(:) = varianceMap(linearIDlistVarianceMap(:));

end
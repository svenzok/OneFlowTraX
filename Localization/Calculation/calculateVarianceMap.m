function [varianceMap, offset] = calculateVarianceMap(fileName, Options)
% calculateVarianceMap calculates the variance per pixel of a noise file.
%
% Syntax:
%   varianceMap = calculateVarianceMap(fileName, Name, Value, ...)
%
% Input Arguments:
%   (Required)
%   fileName           Full file name of a .TIF noise file for SMLM. This
%                      is usually a movie (image stack) with the same
%                      number (or more) frames as the corresponding .TIF
%                      data movies (e.g., the same frame acquisition time),
%                      but without any light reaching the camera (e.g., by
%                      screwing on the protective cap).
%                      If this file only consists of one frame, it is 
%                      assumed to be an already calculated varmap and
%                      |offset| will be set to [].
%                      (1,:) char
%
% Name-Value Pair Input Arguments:
%   (Optional)
%   useGPU             Logical switch to use the GPU if possible (default: true).
%                      (1,1) logical
%
% Output Arguments:
%   varianceMap        Pixel-wise intensity variance with the same height
%                      and width as the input file.
%                      (:,:) double
%
%   offset             Offset (average value of all pixels).
%                      (1,1) double
%
% Other required m-files: readTIFFstack
% Subfunctions: none
% Additional required MATLAB products:
% - Parallel Computing Toolbox 7.6 (for default 'GPU' option)
% - MATLAB Parallel Server 7.6 (for default 'GPU' option)
%
% Notes:
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
% Initial release: 2022-11-16
% Last revision: 2023-08-24

%% Function argument validation
arguments
    fileName (1,:) char
    Options.useGPU (1,1) logical = true
end

%% Main
imageStack = readTIFFstack(fileName);

if size(imageStack, 3) == 1
    % Only one frame, the TIFF file already represents the variance map.
    offset = [];
    varianceMap = imageStack;

else
    if Options.useGPU
        imageStackGPU = gpuArray(imageStack);

        offsetG = mean(imageStackGPU, 'all');

        varmapG = sum((single(imageStackGPU) - ...
            sum(imageStackGPU,3) / size(imageStackGPU,3)).^2,3) / ...
            (size(imageStackGPU,3) - 1);

        offset = gather(offsetG);
        varianceMap = gather(varmapG);
    else
        offset = mean(imageStack, 'all');

        varianceMap = sum((single(imageStack) - ...
            sum(imageStack,3) / size(imageStack,3)).^2,3) / ...
            (size(imageStack,3) - 1);

    end
end

end
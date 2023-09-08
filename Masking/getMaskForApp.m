function [maskIm, maskPolygon] = getMaskForApp(filename)
% createMaskForApp tries to fetch a mask for the given file (for the OneFlowTraX app).
%
% Syntax:
%   mask = createMaskForApp(filename)
%
% Input Arguments:
%   (Required)
%   filename           Name of the SMLM file that the mask will be applied to. The mask file itself must be in the same
%                      folder and must have the same name as the SMLM file (with the suffix '_mask'). The file extension
%                      will be ignored, this program checks if the file is readable.
%                      (:,1) char | (1,1) string
%
% Output Arguments:
%   maskIm             Logical mask with the same width and height as as the input file. Empty if no mask is found.
%                      (:,:) logical
%
%   maskPolygon        The same mask as a polyshape object that can be scaled according to the pixel size. This also
%                      enables to use >isinterior< to check which localized point coordinates in an image lie outside
%                      the mask. Empty if no mask is found.
%                      (1,1) polyshape object
%
% Other required m-files: none
% Subfunctions: none
% Additional required MATLAB products:
% - Image Processing Toolbox 11.7
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
% Initial release: 2023-07-14
% Last revision: 2023-07-14

%% Function argument validation
arguments
    filename (1,:) char
end

%% Main

% In this program, polyshapes will be created for lists of vertices with the same entry at the top and bottom of the
% list (signifying a closed polygon). The polyshapes will still be created correctly, but suppress the corresponding
% warning for duplicate vertices.
warning('off','MATLAB:polyshape:repairedBySimplify');

% Create empty outputs (they are replaced when a mask is found).
maskIm = [];
maskPolygon = [];

% Read the first frame of the input file to get the width and height (for later).
inputImage = imread(filename);

% Add '_mask' to the name and collect all corresponding files (ignore the extension for now). This results in a
% structure.
[pathname, filename] = fileparts(filename);
TestFiles = dir(fullfile(pathname,strcat(filename,'_mask.*')));

if ~isempty(TestFiles)
    % Use >imfinfo< to mark all files that can be read by MATLAB as images.
    for iTest = 1:numel(TestFiles)
        try ImageInfo = imfinfo(fullfile(TestFiles(iTest).folder, TestFiles(iTest).name));
            if isequal([ImageInfo.Height, ImageInfo.Height], size(inputImage))
                % The mask must have the same dimensions as the input image.
                TestFiles(iTest).isValidImage = true;
                TestFiles(iTest).isTIF = strcmp(ImageInfo.Format,'tif');
            else
                TestFiles(iTest).isValidImage = false;
            end

        catch
            TestFiles(iTest).isValidImage = false;
        end
    end

    if sum([TestFiles.isValidImage]) > 0
        % Delete all invalid files from the list.
        TestFiles = TestFiles([TestFiles.isValidImage]);

        % Prioritize the TIF file if it exists (the default output format for the mask file), otherwise get the first
        % valid image in the list.
        if any([TestFiles.isTIF])
            iTest = find([TestFiles.isTIF],1);
        else
            iTest = find([TestFiles.isValidImage],1);
        end

        % Convert it into a logical array by thresholding between the darkest and brightest pixels (use >im2gray< to
        % cover RGB images).
        validMaskFilename = fullfile(TestFiles(iTest).folder, TestFiles(iTest).name);
        maskIm = imread(validMaskFilename);
        maskIm = imbinarize(mat2gray(im2gray(double(maskIm))));

        % Store the mask as a polyshape.
        B = bwboundaries(maskIm, 'CoordinateOrder', 'xy')';
        B(2,:) = repmat({[NaN,NaN]},1,numel(B));
        B = cell2mat(B(:));
        maskPolygon = polyshape(B);
    end
end

% Turn the warning back on.
warning('on','MATLAB:polyshape:repairedBySimplify');

end
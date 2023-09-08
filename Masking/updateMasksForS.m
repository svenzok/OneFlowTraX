function S = updateMasksForS(S)
% createMasksForS updates mask related entries in the parameter structure S based on its stored files.
%
% Syntax:
%   S = createMasksS(S)
%
% Input Arguments:
%   (Required)
%   S                  Structure that holds the batch analysis parameters. See |OneFlowTraX.mlapp| for its required
%                      contents.
%                      (1,1) struct
%
% Output Arguments:
%   S                  Updated structure, with S.tracking.maskPolygons now holding (updated) polyshapes in an Nx1 cell
%                      array, with their corresponding filenames in S.files.maskFiles.
%                      (1,1) struct
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
% Last revision: 2023-08-03

%% Function argument validation
arguments
    S (1,1) struct
end

%% Main

% In this program, polyshapes will be created for lists of vertices with the same entry at the top and bottom of the
% list (signifying a closed polygon). The polyshapes will still be created correctly, but suppress the corresponding
% warning for duplicate vertices.
warning('off','MATLAB:polyshape:repairedBySimplify');

% Get the input file list.
inputFiles = S.files.inputFiles;
nFiles = numel(inputFiles);

% Overwrite the mask filenames and polygons in S with empty cells.
S.files.maskFiles = cell(nFiles,1);
S.tracking.maskPolygons = cell(nFiles,1);

% Loop through each file.
for iFile = 1:nFiles

    % Read the first frame of the input file to get the width and height (for later).
    inputImage = imread(inputFiles{iFile});

    % Add '_mask' to the name and collect all corresponding files (ignore the extension for now). This results in a
    % structure.
    [pathname, filename] = fileparts(inputFiles{iFile});
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

            % Prioritize the TIF file if it exists (the default output format for the mask file), otherwise get the
            % first valid image in the list.
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

            % Store the mask as a polyshape and also store its filename.
            B = bwboundaries(maskIm)';
            B(2,:) = repmat({[NaN,NaN]},1,numel(B));
            B = cell2mat(B(:));
            S.tracking.maskPolygons{iFile} = polyshape(B);
            S.files.maskFiles{iFile} = validMaskFilename;
        end
    end
end

% Turn the warning back on.
warning('on','MATLAB:polyshape:repairedBySimplify');

end
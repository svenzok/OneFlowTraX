function [imageStack, Metadata] = readTIFFstack(filename)
% readTiffStack reads in a multipage TIFF file.
%
% Syntax:
%   [imageStack, imageDescription] = readTIFFstack(filename)
%
% Input Arguments:
%   (Required)
%   filename           File name of the multipage TIFF file. If the file is
%                      not in the current folder or in a folder on the
%                      MATLAB path, then provide the full file name.
%                      (1,:) char
%
% Output Arguments:
%   imageStack         Image stack as a 3D array, one page for each frame
%                      of the multipage TIFF file. The output class depends
%                      on the input data. 
%                      (:,:,:) (uint8, uint16, single, double)
%
%   Metadata           The metadata as read from the first image of the
%                      image stack, additionally, the number of frames will
%                      be stored in the field .Frames
%                      struct
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
% Initial release: 2021-12-29
% Last revision: 2023-09-04

%% Function argument validation
arguments
    filename (1,:) char
end

%% Main

% Disable tifflib warnings to ignore custom TIFF tags.
warning('off', 'imageio:tiffmexutils:libtiffWarning');

% Use an internal MATLAB function to quickly get the image info for all
% frames.
Metadata = matlab.io.internal.imagesci.imtifinfo(filename);

% stack dimensions
nPixelsX = Metadata.Width;
nPixelsY = Metadata.Height;
nFrames = numel(Metadata);
rowsPerStrip = Metadata.RowsPerStrip;
nStrips = nFrames*nPixelsY/rowsPerStrip;

% data precision
BitsPerSample = Metadata.BitsPerSample;

% byte order
byteOrder = Metadata.ByteOrder;
switch byteOrder
    case 'big-endian'
        byteOrder = 'b';
    case 'little-endian'
        byteOrder = 'l';
end

% sample format
if ~isfield(Metadata, 'SampleFormat')
    sampleFormat = 'Unsigned integer';
else
    sampleFormat = Metadata.SampleFormat;
end

switch sampleFormat
    case 'Unsigned integer'
        dataType = strcat('uint',num2str(BitsPerSample));

    case 'IEEE floating point'
        if BitsPerSample == 32
            dataType = 'single';
        else
            error('incorrect TIFF file bit depth descriptors');
        end

    otherwise
        error('incorrect TIFF file bit depth descriptors');
end

% Get the offset for each strip. Images can consist of several strips, we load all strips first and reshape the image
% stack later.
stripOffsets = [Metadata.StripOffsets];

% Pre-allocate according to bit depth.
imageStrips = zeros(rowsPerStrip, nPixelsX, nStrips, dataType);

% Construct precision string for fread (e.g., *uint16).
dataPrecision = strcat('*',dataType);

% Open the file.
fileID = fopen(filename, 'r');

for iStrip = 1:nStrips
    % Jump to the next strip offset to find the data.
    fseek(fileID, stripOffsets(iStrip), 'bof');
    
    % Store each strip.
     imageStrips(:, :, iStrip) = fread(fileID, [rowsPerStrip, nPixelsX], dataPrecision, byteOrder);
end

fclose(fileID);

% Reshape the strips to the image stack.
imageStack = permute(reshape(imageStrips, nPixelsY, nPixelsX, nFrames), [2 1 3]);

% Re-enable tifflib warnings
warning('on', 'imageio:tiffmexutils:libtiffWarning');

% Only store the metadata of the first frame and add the number of frames.
Metadata = Metadata(1);
Metadata.Frames = nFrames;

end
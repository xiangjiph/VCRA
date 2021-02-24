function [outputArgs] = mytiffwrite(outputFileName,imgdata)
%MYTIFFWRITE Summary of this function goes here
%
% [OUTPUTARGS] = MYTIFFWRITE(INPUTARGS) Explain usage here
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/09/04 15:55:07 $	$Revision: 0.1 $
% Copyright: HHMI 2015
if 1
    imwrite(imgdata(:, :, 1), outputFileName, 'WriteMode', 'overwrite',  'Compression','none');
    for K=2:length(imgdata(1, 1, :))
        imwrite(imgdata(:, :, K), outputFileName, 'WriteMode', 'append',  'Compression','none');
    end
else
    dims = size(imgdata);
    obj = Tiff(outputFileName,'w8');
    
    tagstruct.ImageLength = size(imgdata,1)
    tagstruct.ImageWidth = size(imgdata,2)
    tagstruct.Photometric = Tiff.Photometric.RGB
    tagstruct.BitsPerSample = 16
    tagstruct.SamplesPerPixel = 2
    %     tagstruct.RowsPerStrip = 16
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky
    tagstruct.Software = 'MATLAB'
    %     tagstruct.SubIFD = 2  % required to create subdirectories
    obj.setTag(tagstruct)
    
    
end

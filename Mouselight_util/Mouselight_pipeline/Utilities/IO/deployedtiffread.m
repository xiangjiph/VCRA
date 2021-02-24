function [Iout] = deployedtiffread(fileName,slices)
%DEPLOYEDTIFFREAD Summary of this function goes here
% 
% [OUTPUTARGS] = DEPLOYEDTIFFREAD(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/21 12:26:16 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(fileName, 'tif');
if nargin<2
    slices = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(slices);
Iout  = zeros(hIm, wIm, numIm,'uint16');

for i=1:numIm
    Iout(:,:,i) = imread(fileName,'Index',slices(i),'Info',info);
end
warning on
end
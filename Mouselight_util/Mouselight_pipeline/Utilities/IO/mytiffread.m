function [Iout] = mytiffread(FileTif,frames)
%MYTIFFREAD Summary of this function goes here
%
% [OUTPUTARGS] = MYTIFFREAD(INPUTARGS) Explain usage here
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2015/08/11 15:37:18 $	$Revision: 0.1 $
% Copyright: HHMI 2015
warning off
info = imfinfo(FileTif, 'tif');
if nargin<2
    frames = 1:length(info);
end
wIm=info(1).Width;
hIm=info(1).Height;
numIm = numel(frames);
Iout  = zeros(hIm, wIm, numIm,'uint16');

FileID = tifflib('open',FileTif,'r');
for i=1:numIm
    tifflib('setDirectory',FileID,frames(i)-1);
    Iout(:,:,i) = tifflib('readEncodedStrip',FileID,0);
end
warning on
tifflib('close',FileID);
end

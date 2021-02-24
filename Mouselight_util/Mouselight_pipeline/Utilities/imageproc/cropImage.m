function [Isub] = cropImage(I,bbox)
%CROPIMAGE Crops a patch from image cell array I based on bbox 
%
% [OUTPUTARGS] = CROPIMAGE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/09/01 14:32:58 $	$Revision: 0.1 $
% Copyright: HHMI 2015

if iscell(I)
	for j=1:length(I)
	    Isub{j} = I{j}(bbox(1):bbox(2),bbox(3):bbox(4),bbox(5):bbox(6));
	end
else
	Isub = I(bbox(1):bbox(2),bbox(3):bbox(4),bbox(5):bbox(6));
end
end

function [rt1] = normImage(rt1)
%NORMIMAGE Summary of this function goes here
% 
% [OUTPUTARGS] = NORMIMAGE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/09/08 14:46:07 $	$Revision: 0.1 $
% Copyright: HHMI 2015

rt1 = (rt1-min(rt1(:)))/(max(rt1(:))-min(rt1(:)));


end

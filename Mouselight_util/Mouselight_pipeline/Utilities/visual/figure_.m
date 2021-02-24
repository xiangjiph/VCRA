function [hfig] = figure_(idx)
%FIGURE_ Summary of this function goes here
% 
% [OUTPUTARGS] = FIGURE_(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/10/16 17:51:23 $	$Revision: 0.1 $
% Copyright: HHMI 2015
if nargin<1
    hfig = figure;
else
    hfig = figure(idx);
end
subplot('position', [0 0 1 1]);

end

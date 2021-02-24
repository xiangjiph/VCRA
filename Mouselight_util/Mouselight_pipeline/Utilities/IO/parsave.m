function [outputArgs] = parsave(varargin)
%PARSAVE saves inside a parfor loop
% 
% [OUTPUTARGS] = PARSAVE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/09/17 15:29:30 $	$Revision: 0.1 $
% Copyright: HHMI 2015

savefile = varargin{1}; % first input argument

for i=2:nargin
    savevar.(inputname(i)) = varargin{i}; % other input arguments
end
save(savefile,'-struct','savevar')
end

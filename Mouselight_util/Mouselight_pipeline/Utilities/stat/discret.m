function [X] = discret(X,disc)
%DISCRET discretize a integer valued signal
% 
% [OUTPUTARGS] = DISCRET(INPUTARGS) discretize a signal with a given
% interval
% 
% Inputs: 
%   X [Nx1] : 1D signal
%   disc (s): discretization level
% 
% Outputs: 
%   X [Nx1] : 1D signal
% Examples: 
%   X = uint8(rand(100,1)*255)
%   Xest = discret(X,4) % discretize signals with 4 as the intercal width
%   % display discretezation error per index
%   X-Xest
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/09/03 16:53:31 $	$Revision: 0.1 $
% Copyright: HHMI 2015

maxbin = max(X(:));
minbin = min(X(:));
minbin = floor(minbin/disc)*disc;
edges = [minbin:disc:maxbin-1,maxbin];
X = round(edges(discretizemex(X(:),edges,false)));

end

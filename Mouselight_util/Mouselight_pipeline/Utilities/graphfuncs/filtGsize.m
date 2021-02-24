function [G,subs] = filtGsize(G,subs,argin)
%FILTGSIZE Summary of this function goes here
% 
% [OUTPUTARGS] = FILTGSIZE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/05 16:35:55 $	$Revision: 0.1 $
% Copyright: HHMI 2017

sizeThr = argin{1};
CompsC = conncomp(G,'OutputForm','cell');
% S = length(CompsC);
Y = cellfun(@length,CompsC);
deletethese = (Y<sizeThr);
cpl = [CompsC{deletethese}];
G = rmnode(G,cpl);
subs(cpl,:)=[];
end

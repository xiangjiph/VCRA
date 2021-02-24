function [out] = downSamplePath(xyz,sp)
%DOWNSAMPLEPATH Summary of this function goes here
%
% [OUTPUTARGS] = DOWNSAMPLEPATH(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/10/22 11:07:25 $	$Revision: 0.1 $
% Copyright: HHMI 2015
if nargin<2
    sp = 3;
end
if size(xyz,1)<sp % sampling size is bigger then the number of nodes on the path
    out = [1 size(xyz,1)]; % return start and end points
    return
end
p2 = sum((xyz(1:end-1,:)-xyz(2:end,:)).^2,2);
dists = sqrt(p2);
cdist = cumsum(dists);

if 1
    set = 1:length(p2)+1;
    out = set(1);
    dist = 0;
    for jj=1:length(dists)
        if dist+dists(jj)>sp
            out = [out set(1+jj)];
            dist = 0;
        else
            dist = dist + dists(jj);
        end
    end
    if length(out)>2 & dists(end)<3 % overwrite last one. This is for cosmatics
        out(end) = set(end);
    else
        out = [out set(end)];
    end
else
    out = linspace(0,cdist(end),round(cdist(end)/sp)+1);
end



end

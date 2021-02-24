function [intree] = smoothtree(intree,opt)
%SMOOTHTREE Summary of this function goes here
%
% [OUTPUTARGS] = SMOOTHTREE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/03/29 10:17:20 $	$Revision: 0.1 $
% Copyright: HHMI 2016
[L,list] = getBranches(intree.dA);


if isfield(opt,'sizethreshold')
    filtwidth = opt.sizethreshold;
else
    myhist = histc(list,unique(list));
    filtwidth = max(3,round(median(myhist)/10));
end

XYZ = [intree.X intree.Y intree.Z];
R = intree.R;
D = intree.D;

for ii=1:length(L)
    set_ii = L(ii).set(2:end-1);
    if isempty(set_ii)
        continue
    end
    xyz = XYZ(set_ii,:);
    for jj=3 % only on z
        X = xyz(:,jj);
        XYZ(set_ii,jj) = medfilt1(medfilt1(X,filtwidth),filtwidth);
    end
    intree.Z(set_ii) = XYZ(set_ii,jj);
    % smooth radius
    if 1
        R(set_ii) = medfilt1(medfilt1(R(set_ii),3),3);
    end
end
end

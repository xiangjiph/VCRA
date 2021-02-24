function [inupdate,deleteThese] = prunTree(in_,lengthThr,res)
%PRUNTREE Pruns a tree with a given length threshold
% 
% [OUTPUTARGS] = PRUNTREE(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/10/28 10:42:45 $	$Revision: 0.1 $
% Copyright: HHMI 2015

A = in_.dA;
XYZ = [in_.X*res(1) in_.Y*res(2) in_.Z*res(3)];% in (um)
[L,list] = getBranches(A);
numBranches = length(L);
% find leaf branches
termnodes = find(sum(A)==0);
leafbranches = zeros(1,numBranches);
deleteThese = [];
for ii=1:numBranches
    Liiset = L(ii).set;
    if ~isempty(Liiset) & any(termnodes==Liiset(1)) %& length(Liiset)<sizeThr
        lenBranch = sum(sqrt(sum(diff(XYZ([Liiset L(ii).parentnode],:)).^2,2)));
        if lenBranch<lengthThr
            leafbranches(ii) = 1;
            deleteThese = [deleteThese L(ii).set];
        end
    end
end

inupdate = in_;
inupdate.dA(deleteThese,:) = [];
inupdate.dA(:,deleteThese) = [];
inupdate.X(deleteThese) = [];
inupdate.Y(deleteThese) = [];
inupdate.Z(deleteThese) = [];
inupdate.R(deleteThese) = [];
inupdate.D(deleteThese) = [];


end

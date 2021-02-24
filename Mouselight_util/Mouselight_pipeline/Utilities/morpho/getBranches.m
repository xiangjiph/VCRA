function [branch,list] = getBranches(G,rootnode)
%GETBRANCHES Given directed graph G, returns list of branch segments
% 
% [OUTPUTARGS] = GETBRANCHES(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/10/21 17:00:54 $	$Revision: 0.1 $
% Copyright: HHMI 2015

% force lower triangle
% G_ = G.*tril(G,-1);

% critical nodes
if nargin<2
    rootnode = find(sum(G,2)==0);
end
branchnodes = find(sum(G)>1);
% root is branch if it splits into 3
if sum(G(:,rootnode))<3
    branchnodes = setdiff(branchnodes,rootnode);
end
termnodes = find(sum(G)==0);
% from all branch and termnodes, get segmentsuntill another critical point
% is found
criticalset = [rootnode, branchnodes, termnodes];

% get pred
[DISC,PRED,CLOSE] = graphtraverse(G',rootnode,'DIRECTED',true);
%[DIST2,PATH2,PRED2] = graphshortestpath(G',rootnode,'DIRECTED',true);

%%
criticfalPRED = PRED;
criticfalPRED(criticalset) = 0;
iter = 0;
clear branch
list = zeros(1,size(G,1));
for idx = criticalset
    %%
    set = idx;
    next = PRED(idx);
    while next
        set = [set next];
        next = criticfalPRED(next);
    end
    iter = iter+1;
    branch(iter).set = set(1:end-1); % exclude the parent critical point
    branch(iter).parentnode = set(end);
    if isempty(branch(iter).set)
        list(idx) = iter; % root node
    else
        list(branch(iter).set) = iter;
    end
end
% parent branch ids
for ii = 1:length(branch)
    branch(ii).parentbranch = list(branch(ii).parentnode);
end
%%

end

























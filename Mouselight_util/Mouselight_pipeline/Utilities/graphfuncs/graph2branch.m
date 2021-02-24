function [branches,nodeBrid] = graph2branch(G,subs)
%GRAPH2BRANCH Given graph, finds all branches
% 
% [OUTPUTARGS] = GRAPH2BRANCH(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/03/29 13:22:54 $	$Revision: 0.1 $
% Copyright: HHMI 2017
% G=graph(Ain);
%%
numnodes = size(subs,1);
% first find components
CompsC = conncomp(G,'OutputForm','cell');
Comps = conncomp(G);
A = G.adjacency;
S = length(CompsC);
Y = cellfun(@length,CompsC);
leafnodes = find(sum(G.adjacency,2)==1);
regnodes = find(sum(G.adjacency,2)==2);
juncnodes = find(sum(G.adjacency,2)>2);
cricset=union(leafnodes,juncnodes);
%%
clear Gcomps BR;
mytic = tic;
try;parfor_progress(0);catch;end
tic
parfor_progress(S)
parfor idxC = 1:S
    parfor_progress;
    A_ = A(CompsC{idxC},CompsC{idxC});
    subs_ = subs(CompsC{idxC},:);
    BR{idxC} = graphfuncs.findbranch(A_,subs_);
end
parfor_progress(0);
sprintf('branches found in %d sec for %d comps',round(toc(mytic)),S)

%%
branches = [];
nodeBrid = cell(1,numnodes);
iter = 1;
for idxC = 1:S
    for ii=1:length(BR{idxC})
        B = BR{idxC}(ii);
        for fn = fieldnames(B)'
            branches(iter).(fn{1}) = B.(fn{1});
        end
        branches(iter).idxC = idxC;
        branches(iter).inds = CompsC{idxC}(branches(iter).set);
        branches(iter).tips = CompsC{idxC}(branches(iter).set([1 end]));
        branches(iter).tipsubs =  branches(iter).subs([1 end],:);
        for jj=1:length(branches(iter).inds)
            nodeBrid{branches(iter).inds(jj)}(end+1) = iter;
        end
        iter = iter + 1;
    end
end

% add some statistics
for ii=1:length(branches)
    branches(ii).cent = mean(branches(ii).tipsubs);%median(branches(ii).subs,1);
end




end

function [eout] = buildgraph(Asub,source)
algorithmkeys = {'spb','dij','bel','spa'};
algorithm = 2;
debug_level = 0;
directed = 0;
if nargin<2
    source = find(sum(Asub)==0,1);
end
if strcmp(version('-release'),'2015b')
    [dist,pred] = graphalgs(algorithmkeys{algorithm},debug_level,directed,...
        max(Asub,Asub'),source);
elseif strcmp(version('-release'),'2017a')
    AA = graph(max(Asub,Asub'));
    TR = shortestpathtree(AA,source);
    %%
    %pred = zeros(1,TR.numnodes);
    %for ii=1:TR.numnodes
    %    pred(TR.successors(ii))=ii;
    %end
    %%
    pre = cell(1,TR.numnodes);
    parfor ii=1:TR.numnodes
        pre{ii} = TR.successors(ii);
    end
    pred = zeros(1,TR.numnodes);
    for ii=1:TR.numnodes
        pred(pre{ii}) = ii;
    end
else
    [dist,path,pred] = graphshortestpath(max(Asub,Asub'),source,'directed','false');
end

eout = [([1:size(Asub,1)]') (pred(:))];
eout(eout(:,2)==0,:) = [];



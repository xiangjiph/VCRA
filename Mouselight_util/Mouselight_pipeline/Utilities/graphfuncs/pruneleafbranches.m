function [Gin,subsin,pruned] = pruneleafbranches(Gin,subsin,sizeThr,deletesmallloops)
if nargin<4
    deletesmallloops=0;
end
leafnodes = find(sum(Gin.adjacency,2)==1);
numleafs = length(leafnodes(:));

regnodes = find(sum(Gin.adjacency,2)==2);
juncnodes = find(sum(Gin.adjacency,2)>2);
A = Gin.adjacency;
% delete junction node edges
A(juncnodes,:) = 0;
A(:,juncnodes) = 0;
%%
CompsC = conncomp(graph(A),'OutputForm','cell');
S = length(CompsC);
clustsize = cellfun(@length,CompsC);
Comps = nan(1,size(A,1));
for idx=1:S
   Comps(CompsC{idx}) = idx;
end
leafbranchclusterID = Comps(leafnodes);
leafbranchclustsize = clustsize(leafbranchclusterID);

%%
% leafadjacency: if there is a nearby leaf, dont delete it
CompsGin = conncomp(Gin);
YGin = histcounts(CompsGin,1:max(CompsGin)+1);
IDX = rangesearch(subsin(leafnodes,:),subsin(leafnodes,:),sizeThr);
len = cellfun(@length,IDX);

%%
lenfilter = len;
for ii=1:numleafs
    qclust = CompsGin(leafnodes(ii));
    % 1:length(CompsC)
    if len(ii)>1
        idxq = IDX{ii}(1:end); 
        leafclustid = CompsGin(leafnodes(idxq));
        [ia,ib,ic] = unique(leafclustid);
        %[adjleafs] = setdiff((CompsGin(leafnodes(idxq))),qclust);
        if length(ia)<2 % all nearby leafs belong to same object, prune
            lenfilter(ii) = 1;
        elseif max(clustsize(leafbranchclusterID(idxq(ic~=find(ia==qclust))))) > sizeThr % check max size of branch of adjecent cluster
            % keep
        else
            lenfilter(ii) = 1;
        end
    end
end

deletetheseleafs = find(leafbranchclustsize<sizeThr & lenfilter(:)'==1);
deletethesebranches = leafbranchclusterID(deletetheseleafs);
%%
cpl = [CompsC{deletethesebranches}];
pruned = length(cpl);
Gin = rmnode(Gin,cpl);
subsin(cpl,:)=[];
if deletesmallloops % check if end points are shared
end












function [L] = extractBranches(edge,numNodes)
%EXTRACTBRANCHES Summary of this function goes here
% 
% [OUTPUTARGS] = EXTRACTBRANCHES(INPUTARGS) Explain usage here
% 
% Examples: 
% 
% Provide sample usage code here
% 
% See also: List related files here

% $Author: base $	$Date: 2015/08/17 12:48:27 $	$Revision: 0.1 $
% Copyright: HHMI 2015

if size(edge,1)==size(edge,2) & size(edge,1)>2 % square adjecency matrix
    [link(:,1),link(:,2)] = find(edge);
else
    link = edge;
end
numNodes = max(link(:));
% get the branch nodes
eS = sparse(link(:,1),link(:,2),1,numNodes,numNodes);
% find branch points
nodeType = sum(eS);
nodeType(nodeType>2) = 2;
eS(nodeType==2,:)=0;
eS(:,nodeType==2)=0;
eS = eS+speye(numNodes,numNodes);
[CC,ss] = graphconncomp(eS+eS');
clear L
% node can be
% anchor -> single set
% continuous-> single
for idx = 1:CC
    set = find(ss==idx);
    % parent node
    anchor_node = link(find(min(set)==link(:,1)),2);
    term_nodes = link(find(max(set)==link(:,2)),1)';
    if idx==1 & isempty(term_nodes)
        % no branch
        L(idx).set = 1:size(link,1)+1;
        L(idx).anchorNode = min(L(idx).set);
        L(idx).termNode = max(L(idx).set);
        return
    end
    if isempty(anchor_node)
        continue
    end
    if length(set) >1
        % default linking
        %child nodes (not needed to create set)
        L(idx).set = [anchor_node set term_nodes];
    else
        if nodeType(set)==1 % continuous
            L(idx).set = [anchor_node set term_nodes];
        elseif nodeType(anchor_node)==1 % skip
            continue
        else
            L(idx).set = [anchor_node set];
        end
    end
    L(idx).anchorNode = min(L(idx).set);
    L(idx).termNode = max(L(idx).set);
end
% delete any empty branches
L(cellfun(@isempty,{L.set}))=[];
%%
pairs = [L.anchorNode;L.termNode];
% built branch tree
CC=[];
for ix = 1:length(pairs)
    anc = pairs(:,ix);
    for iy = ix+1:length(pairs)
        tar = pairs(:,iy);
        % check if they are connected
        test = ~(anc-tar([2 1]));
        if any(test)
            if test(1)
                CC(end+1,:)= [ix iy];
            else
                CC(end+1,:)= [iy ix];
            end
        end
    end
end
% brE = sparse(CC(:,1),CC(:,2),1,length(L),length(L));
for idx = 1:size(CC,1)
    L(CC(idx)).parent = CC(idx,2);
end









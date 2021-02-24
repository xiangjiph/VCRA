function [G,subs,totdel] = deleteLoops(G,subs,argin)
%DELETELOOPS Summary of this function goes here
%
% [OUTPUTARGS] = DELETELOOPS(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2017/04/05 16:42:19 $	$Revision: 0.1 $
% Copyright: HHMI 2017
%%
% first find components
CompsC = conncomp(G,'OutputForm','cell');
A = G.adjacency;
S = length(CompsC);
Y = cellfun(@length,CompsC);
clear Gcomps;
deletethese = cell(1,S);
deleteedges = cell(1,S);
%%
mytic=tic;
%%
try;parfor_progress(0);catch;end
tic
parfor_progress(S)
parfor idxC = 1:S
    parfor_progress;
    indsC = CompsC{idxC};
    subs_ = subs(indsC,:);
    if length(indsC)<500 % x100 faster for large A & small |indsC|
        [aa,bb]=ndgrid(indsC,indsC);
        A_ = reshape(A(sub2ind(size(A),aa(:),bb(:))),length(indsC),length(indsC));
    else
        A_ = A(indsC,indsC);
    end
    if nnz(A_)/2>size(A_,1)-1
        % there is a loop        % check for small loops
        [inds,edges] = graphfuncs.findloops(A_,subs_);
        if any(inds); delinds = indsC(inds); deletethese{idxC} = delinds; end
        if any(edges) deledges = indsC(edges); deleteedges{idxC} = deledges; end
    end
end
parfor_progress(0);
deleteedges = cat(1,deleteedges{:});
delthese = [deletethese{:}];
toc
% %%
% [A_,subs_] = deal([]);
% parfor idxC = 1:S
%     if ~rem(idxC,round(S/100))
%         idxC
%     end
%     indsC = CompsC{idxC};
%     if length(indsC)<500 % x100 faster
% %         tic
%         [aa,bb]=ndgrid(indsC,indsC);
%         A_{idxC} = reshape(A(sub2ind(size(A),aa(:),bb(:))),length(indsC),length(indsC));
% %         toc1 = toc;
%     else
% %         tic
%         A_{idxC} = A(indsC,indsC);
% %         toc2 = toc;
%     end
%     subs_{idxC} = subs(indsC,:);
% end
% %%
% parfor idxC = 1:S
%     indsC = CompsC{idxC};
%     if nnz(A_)/2>size(A_,1)-1
%         % there is a loop        % check for small loops
%         [inds,edges] = graphfuncs.findloops(A_{idxC},subs_{idxC});
%         if any(inds); delinds = indsC(inds); deletethese{idxC} = delinds; end
%         if any(edges) deledges = indsC(edges); deleteedges{idxC} = deledges; end
%     end
% end
% deleteedges = cat(1,deleteedges{:});
% delthese = [deletethese{:}];

%%
if ~isempty(deleteedges)
    [deleteedges_] = unique([deleteedges;deleteedges(:,[2 1])],'rows');
    G = rmedgeBig(G,deleteedges_(:,1),deleteedges_(:,2));
end
if ~isempty(delthese)
    G = rmnode(G,delthese);
    subs(delthese,:)=[];
end
totdel = length(delthese);
sprintf('it took %d secs to delete %d loops in %d [num edges]',round(toc(mytic)),length(deleteedges),numedges(G))
%%
end

function G = rmedgeBig(G, s, t)
%RMEDGE Remove edges from a large graph. This is a fix for matlab's
%built in rmedge function that is memory efficient for large scale graphs
nn = numnodes(G);
A = G.adjacency;
Amask = sparse([s;t],[t;s],1,nn,nn)>0;
G = graph(A-A.*Amask);
end

function [intree] = im2graph(IM,debug)
%IM2GRAPH Converts an image to an affinity graph
%
% [OUTPUTARGS] = IM2GRAPH(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2015/10/20 10:11:33 $	$Revision: 0.1 $
% Copyright: HHMI 2015
intree = [];
if nargin<2
    debug = 0;
end
IM = IM>0;
% pad zeros
IM = padarray(IM,[1 1 1],0,'both');
dims = size(IM);
nout = length(dims);

%Compute linear indices
k = [1 cumprod(dims(1:end-1))];
x = [-1:1];
out = zeros(nout^nout,nout);
siz = nout*ones(1,nout);
for i=1:nout
    s = ones(1,nout);
    s(i) = numel(x);
    x = reshape(x,s);
    s = siz;
    s(i) = 1;
    dum = repmat(x,s);
    out(:,i) = dum(:);
end
idxneig = out*k';
idxneig((numel(idxneig)+1)/2)=[];

%%
%Create affinity matrix
idx = find(IM);
if numel(IM)/2<length(idx)
    warning('dense image, might take sime time to build graph. Consider to make input IM sparse')
end

1
clear E
it = 1;
IM = skel;
for i=idx(:)'
    if debug
        if ~rem(it,max(1,round(numel(idx)/10)))
            sprintf('Building graph: %%%d',round(100*it/numel(idx)))
        end
    end
    id = i+idxneig;
    neigs = find(IM(id));
    
    E{it} = [[i*ones(length(neigs),1) id(neigs)] ]';%-(sum(k)+1)
    it = it+1;
end
edges = [E{:}]';
%%
[keepthese,ia,ic] = unique(edges(:));
[subs(:,2),subs(:,1),subs(:,3)] = ind2sub(dims,keepthese);
edges_ = reshape(ic,[],2);
if isempty(edges_)
    return
end
%%
% permute edges
A = sparse(edges_(:,1),edges_(:,2),1,max(edges_(:)),max(edges_(:)));
A = max(A',A);
A(find(speye(size(A)))) = 0;

[S,C] = graphconncomp(A);
Y = histcounts(C,1:S+1);
A_ = A;
A_(find(triu(A_,0)))=0;

%%
ly = 1:length(Y);%find(Y>1);
clear Eout
for ii=1:length(ly)
    %%
    subidx = find(C==ly(ii));
    Asub = A_(subidx,subidx);
    leafs = find(sum(Asub,2)==0);%find(sum(max(Asub,Asub'))==1,1);
    [dist,path,pred] = graphshortestpath(Asub,leafs(1),'directed','false');
%     [sorteddist,idxdist]=sort(dist);
%     Asub_ = Asub(idxdist,idxdist);
%     [dist,path,pred] = graphshortestpath(Asub_+Asub_',1);
    
%     % sort pred to make it lower triu
%     currind=0;
%     while currind<max(pred)+1
%         inds=find(pred==currind);
%         eo{currind+1} = inds
%         
%         
%     end
%     for jj=1:max()
%     [sortedpred,sortedindx] = sort(pred);
%     eout = [([1:size(Asub,1)]') ([1:size(Asub,1)]')-1];
%     eout(eout(:,2)==0,:) = [];
%     % permute subs and indicies
%     subidx = subidx(sortedindx)
    % create an affinity map
    eout = [([1:size(Asub,1)]') (pred(:))];
    eout(eout(:,2)==0,:) = [];
%     Eout{ii} =subidx(sort(eout,2,'descend'))';
    %if any(any(subidx(eout)'==14820));ii,end
    Eout{ii} =subidx(eout)';
end
% make sure that graph is lower triangular
subs_ = subs;
edgesout = [Eout{:}]';
%%
% while true
%     idx = find(edgesout(:,1)<edgesout(:,2),1);
%     if isempty(idx)
%         break
%     end
%     [aa] = edgesout(idx,:);
%     % flip axis
%     edgesout(idx,:) = aa([2 1]);
%     tmp=subs_(aa(1),:);
%     subs_(aa(1),:) = subs_(aa(2),:);
%     subs_(aa(2),:) = tmp;
% end


%%
% Io = Io(2:end-1,2:end-1,2:end-1);
%
Aout = sparse(edgesout(:,1),edgesout(:,2),1,max(edgesout(:)),max(edgesout(:)));
% [S,C] = graphconncomp(Aout,'directed','false');

intree.dA = Aout;
intree.R = ones(size(A_,1),1);
intree.D = ones(size(A_,1),1);
intree.X = subs_(:,1)-1; % subtract padsize
intree.Y = subs_(:,2)-1; % subtract padsize
intree.Z = subs_(:,3)-1; % subtract padsize

























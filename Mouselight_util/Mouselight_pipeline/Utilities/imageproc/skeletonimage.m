function [skel,A,subs_,edges_,weights_] = skeletonimage(Io,opt)
%SKELETONIMAGE Summary of this function goes here
%   Detailed explanation goes here

[skel,A,subs_,edges_] = deal([]);

if isfield(opt,'thr')
    probThr = opt.thr;
else
    % non-zero min 
    non_zero_min = min(Io(Io>0));
    probThr = getThresh(Io(Io>non_zero_min));
end

if isfield(opt,'fullh')
    fullh = opt.fullh;
else
    fullh = min(size(Io),5);
end

if isfield(opt,'anisotropy')
    anisotropy = opt.anisotropy;
else
    anisotropy = [3 3 1];
end

if isfield(opt,'sizethreshold')
    sizethreshold = opt.sizethreshold;
else
    sizethreshold = 0;
end
% smooth image
Io = smooth3(Io,'gaussian',anisotropy);
Io = Io>probThr;

%%
if ~any(Io(:));warning('Empty volume');return;end

%% cleanup image
s  = regionprops(Io, 'centroid','PixelIdxList','Area');
% memory efficient
for ii=1:length(s)
    if s(ii).Area<sizethreshold
        Io(s(ii).PixelIdxList)=0;
    end
end

%%
if ~any(Io(:));warning('Empty volume');return;end

%%
% binarize it before skeletionization
Io = Io>0;
% run skeletonization
% if size(Io) is big limit memory by using less number of nodes
skel = block3d({Io},[200 200 200],fullh,1,@Skeleton3D,[]);
% estimate radius
skel = padarray(skel,ones(1,3),0,'both');
Io = padarray(Io,ones(1,3),0,'both');
bIo=bwdist(~Io);
radskel = double(bIo.*single(skel));

%%
% get the edge pairs
dims = size(skel);
skelinds = find(skel);
if isempty(skelinds);warning('Empty volume');return;end
%%
nout = length(dims);
%Compute linear indices
k = [1 cumprod(dims(1:end-1))];
x = [-1:1];
per = zeros(nout^nout,nout);
siz = nout*ones(1,nout);
for i=1:nout
    s = ones(1,nout);
    s(i) = numel(x);
    x = reshape(x,s);
    s = siz;
    s(i) = 1;
    dum = repmat(x,s);
    per(:,i) = dum(:);
end
idxneig = per*k';
idxneig((numel(idxneig)+1)/2)=[];

%% get edge pairs
E = [];
it = 1;
for idx = skelinds(:)'
    inds = idx + idxneig;
    hits = inds(skel(inds));
    rad = radskel(idx);
    % crop back to original size
    E{it} = [[idx*ones(length(hits),1) hits(:) rad*ones(length(hits),1)]]';
    it = it+1;
end
edges = [E{:}]'; clear E
%% map onto original graph
for ii=1:2
    [xx,yy,zz] = ind2sub(dims,edges(:,ii)); % subs on appended crop
    subs = [xx(:),yy(:),zz(:)]-1; % to compansate crop;
    edges(:,ii) = sub2ind(dims-1,subs(:,1),subs(:,2),subs(:,3));
end
% [xx,yy,zz] = ind2sub(dims,edges(:,1)); % subs on appended crop
% subs = [xx(:),yy(:),zz(:)]-1; % to compansate crop;

clear subs_
[keepthese,ia,ic] = unique(edges(:,[1 2]));
[subs_(:,1),subs_(:,2),subs_(:,3)] = ind2sub(dims-1,keepthese);
edges_ = reshape(ic,[],2);
weights_ = edges(ia,3:end);
A = sparse(edges_(:,1),edges_(:,2),1,max(edges_(:,2)),max(edges_(:,2)));
skel = skel(2:end-1,2:end-1,2:end-1);
end



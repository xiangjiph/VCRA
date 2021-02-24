function [pixinit, nummatches, varargout] = initTiles(featmap,directions,scopeloc,numthr)
% estimates initial shift values 
if nargin<4
    numthr = 50;
end
numtiles = length(featmap);
pixinit = nan(numtiles,3);
nummatches = zeros(numtiles,1);
%%
for ii=1:numtiles
    feat = featmap(ii).(genvarname(directions));
    if isempty(feat)
    elseif isempty(feat.paireddescriptor.Y)
    else
        difvec = feat.paireddescriptor.X-feat.paireddescriptor.Y;
        pixinit(ii,:) = median(difvec,1);
        nummatches(ii) = size(difvec,1);
    end
end
if nargout == 3
    varargout{1} = pixinit;
end
%%
% build kdtree
% gridix = scopeloc.gridix(:,1:3);
inliers = find(all(isfinite(pixinit),2) & nummatches>numthr);
anchors = scopeloc.gridix(inliers,1:3);
queries = scopeloc.gridix(:,1:3);
IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]
pixinit = pixinit(inliers(IDX),:);

function bbox = createOverlapBox(dims,blockSize,fullh)
if nargin<1
    dims = [2000 2000 2000];
    blockSize = [200 200 200];
    fullh = 15;
end
%%
if length(fullh)<3
    fullh = ones(1,3)*fullh;
end
h = (fullh-1)/2;
numDims = length(dims);
[BB{1:numDims}]=deal([]);
for i=1:numDims
    st = 1:(blockSize(i)-fullh(i)):dims(i)-h(i)-1;
    en = min(dims(i),st+(blockSize(i)-1));
    BB{i} = [st(:) en(:)];
end
numb = cellfun(@(x) size(x,1),BB);
if length(numb)==2
    numb = [numb 1];
end
[aa,bb,cc] = ndgrid(1:numb(1),1:numb(2),1:numb(3));
perms = [aa(:),bb(:),cc(:)];
numblocks = size(perms,1);
bbox = zeros(numblocks,2*numDims);
for i=1:numblocks
    bbox(i,:) = [BB{1}(perms(i,1),:) BB{2}(perms(i,2),:) BB{3}(perms(i,3),:)];
end

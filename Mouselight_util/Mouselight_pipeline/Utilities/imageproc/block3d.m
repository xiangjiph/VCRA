function varargout = block3d(I,blockSize,fullh,parallel,blockfun, blockfunarg)
%BLOCK3D Summary of this function goes here
%
% [OUTPUTARGS] = BLOCK3D(INPUTARGS) runs pixelwise operations on
% overlapping 3d blocks.
%
% Inputs:
%   {I}: image cell array that will be processed
%   blockSize [H,W,D] : size of each block
%   fullh (s) : filter size in each dimension
%   parallel {"0",1} : flag to parallelize process
%   blockfun (@) : filter function handle
%   blockfunarg {} : filter arguments that will be passed

%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: blockproc

% $Author: base $	$Date: 2015/09/01 14:31:51 $	$Revision: 0.1 $
% Copyright: HHMI 2015

h = (fullh-1)/2;
if ~iscell(I)
    I = {I};
end
dimsori = size(I{1});
numDims = ndims(I{1});
% outSize = [dimsori nbins];
numChannel = length(I);
for i=1:numChannel
    % I{i} = padarray(I{i},h*ones(1,numDims),0,'both');
    I{i} = padarray(I{i},h*ones(1,numDims),'symmetric','both');
end
dims = size(I{i});
%%
% split data into overlapping blocks
[BB{1:numDims}]=deal([]);
for i=1:numDims
    st = 1:(blockSize(i)-2*h):dims(i)-h-1;
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
%%
[OUT_{1:numblocks}]=deal([]);
if parallel
    parfor i=1:numblocks
        % crop image
        Isub = cropImage(I,bbox(i,:));
        % run the function for this block
        if numChannel>1
            OUT_{i} = feval(blockfun,Isub,blockfunarg);
        else
            OUT_{i} = feval(blockfun,Isub{1},blockfunarg);
        end
    end
    if isdeployed
        poolobj = gcp('nocreate');
        delete(poolobj)
    end
else
    for i=1:numblocks
        if ~rem(i,round(numblocks/10))&0
            disp(['running: ',num2str(i), ' out of ', num2str(numblocks)])
        end
        % crop image
        Isub = cropImage(I,bbox(i,:));
        % run the function for this block
        if numChannel>1
            OUT_{i} = feval(blockfun,Isub,blockfunarg);
        else
            OUT_{i} = feval(blockfun,Isub{1},blockfunarg);
        end
        if ~rem(i,round(numblocks/10))&0
            disp(['finished',num2str(i)])
        end
    end
end
%%
numOut = nargout;
output = cell(1,numOut);
if iscell(OUT_{1}) % multiple output
    for iout = 1:numOut
        output{iout} = zeros(dimsori,class(OUT_{iout}));
    end
else
    numOut = 1;
    if islogical(I{1})
        [output] = deal(false(dimsori));
    else
        [output] = deal(zeros(dimsori,class(OUT_{1})));
    end
end

% stitch back to original size
for i=1:numblocks
    if numOut>1
        for j=1:numOut
            output{j}(bbox(i,1):bbox(i,2)-2*h,...
                bbox(i,3):bbox(i,4)-2*h,...
                bbox(i,5):bbox(i,6)-2*h) = OUT_{i}{j}(h+1:end-h,h+1:end-h,h+1:end-h);
        end
    else
        output(bbox(i,1):bbox(i,2)-2*h,...
            bbox(i,3):bbox(i,4)-2*h,...
            bbox(i,5):bbox(i,6)-2*h) = OUT_{i}(h+1:end-h,h+1:end-h,h+1:end-h);
    end
end
if ~iscell(output)
    varargout{1} = output;
end
end

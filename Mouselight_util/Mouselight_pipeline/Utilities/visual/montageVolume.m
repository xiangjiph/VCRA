function [outputArgs] = montageVolume(inVolume,labels,params)
%MONTAGEVOLUME Creates a montaged image for browsing. This is extended
%version of Piotr's montage2 function
%
% [OUTPUTARGS] = MONTAGEVOLUME(INPUTARGS)
%
% Examples:
%
% Provide sample usage code here
%
% See also: MONTAGE, MONTAGE2

% $Author: base $	$Date: 2015/08/27 14:05:27 $	$Revision: 0.1 $
% Copyright: HHMI 2015
if nargin<2
    labels = [];
end

% make sure input is in single precision
inVolume = single(inVolume);

% setup output
dims = size(inVolume);
if length(dims)==3
    % montage gray scale image
    [ImOut numColumns]= tileIm(inVolume);
    
elseif length(dims)==4
    % montage RGB image
    projection = @max;
    
    % normalize images
    if params.normalize
        for i=1:dims(4)
            in = single(inVolume(:,:,:,i));
            inVolume(:,:,:,i) = (in-min(in(:)))/(max(in(:))-min(in(:)));
        end
    end
    
    Im1=squeeze(projection(inVolume,[],1));
    Im2=squeeze(projection(inVolume,[],2));
    Im3=squeeze(projection(inVolume,[],3));
    
    [tiledIm1 numColumns] = tileIm(Im1);
    tiledIm2 = tileIm(Im2);
    tiledIm3 = tileIm(Im3);
    ImOut = cat(3,tiledIm1,tiledIm2,tiledIm3);
    
else
    
end
%%
pos=get(0,'MonitorPositions');
try
    h = figure_;
catch
    h = figure;
    subplot('position', [0 0 1 1])
end
imshow(ImOut,[])
set(h,'Position',[round(pos(2,1)/2)-10 1 abs(round(pos(2,3))/2) abs(round(pos(2,4)))])
%%
% overlay lines
% draw lines separating frames
if( params.showLines )
    montageWd = numColumns * dims(1) + .5;
    montageHt = numColumns * dims(1) + .5;
    for i=1:numColumns-1
        ht = i*dims(1) +.5; 
        line([.5,montageWd],[ht,ht],'Color','w');
        wd = i*dims(1) +.5; 
        line([wd,wd],[.5,montageHt],'Color','w');
    end
    
    % cross out unused frames
    [nns,mms] = ind2sub( [numColumns,numColumns], dims(4)+1 );
    for i=mms-1:numColumns-1
        for j=nns-1:numColumns-1,
            rStr = i*dims(1); rs = [rStr,rStr+dims(1)];
            cStr = j*dims(2); cs = [cStr,cStr+dims(2)];
            line( cs, rs,'Color','r' );  line( fliplr(cs), rs,'Color','r' );
        end
        nns = 1; % reset the column counter to start from beginning
    end
end
%%
% link the callback function
datacursormode on
if ~isempty(labels)
    labelData.labels = labels;
    labelData.numColumns = numColumns;
    labelData.dims = dims;
    %labelData.frames = params.frames;
    %labelData.scales = params.scale;
    %labelData.refloc = 
    set(gca,'UserData',labelData);
end
dcm_obj = datacursormode(gcf);
set(dcm_obj,'Enable','on',...
    'UpdateFcn',@montageCallback)





end

%%
function [Iout numColumns] = tileIm(imIn)
dims = size(imIn);
numColumns = ceil(sqrt(dims(3)));

% creat checkerboard
dimOut = [dims(1)*numColumns dims(2)*numColumns];
Iout = zeros(dimOut);
for i=1:dims(3)
    [x,y] = ind2sub([numColumns,numColumns],i);
    xsp = (x-1)*dims(2)+1:x*dims(2);
    ysp = (y-1)*dims(1)+1:y*dims(1);
    Iout(ysp,xsp) = imIn(:,:,i);
end
end









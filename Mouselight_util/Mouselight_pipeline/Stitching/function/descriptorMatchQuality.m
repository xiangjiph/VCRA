function [outputArgs] = descriptorMatchQuality(regpts,scopeparams,scopeloc,videofile)
%DESCRIPTORMATCHQUALITY Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCHQUALITY(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/11/15 11:13:34 $	$Revision: 0.1 $
% Copyright: HHMI 2016

%%
if nargin<4
    videofile = 'myvideo.avi';
end
% populate afftile
numTiles = size(regpts,2)
afftile = zeros(3, 4, numTiles);
afftransform = scopeparams.affineglFC;
mflip = eye(3);%mflip(3,3)=1;
afftransform = afftransform*mflip;
for idxt = 1:numTiles
    %%
    % form an affine matrix
    tform = [afftransform scopeloc.loc(idxt,:)'*1e6];
    afftile(:,:,idxt) = tform;
end
%%
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
checkthese = [1 4 5 7]; % 0 - below
D = zeros(1,size(neighbors,1));
for ii=1:size(neighbors,1)
    if any(isnan(neighbors(ii,[1 7])))
        D(ii) = nan;
        continue
    end
    D(ii)=diff(scopeloc.loc(neighbors(ii,[1 7]),3));
end

%%
latticeZRange = unique(scopeloc.gridix(:,3));
clear xyz_t xyz_tp1 XYZ_t XYZ_tp1
for t = latticeZRange(1:end-1)'
    %%
    disp(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
    ix = (scopeloc.gridix(:,3)'==t);
    %%
    cnt = 0;
    Npts = 0;
    clear layer layerp1 poslayer_t poslayer_tp1
    [layer layerp1 poslayer_t poslayer_tp1]=deal([]);
    for id_ix = find(ix)
        %%
        if ~isfield(regpts{id_ix},'X') | isempty(regpts{id_ix}.X)
            continue
        end
        %%
        cnt = cnt+1;
        neigs = regpts{id_ix}.neigs;
        layer{cnt} = regpts{id_ix}.X;
        layerp1{cnt} = regpts{id_ix}.Y;
        % apply tforms
        Layer_t = [regpts{id_ix}.X ones(size(regpts{id_ix}.X,1),1)]*afftile(:,:,neigs(1))';
        Layer_tp1 = [regpts{id_ix}.Y ones(size(regpts{id_ix}.Y,1),1)]*afftile(:,:,neigs(4))';
        Npts = Npts+size(Layer_t,1);
        poslayer_t{cnt} = Layer_t;
        poslayer_tp1{cnt} = Layer_tp1;
    end
    %%
    if isempty(layer)
        continue
    end
    xyz_t{t} = cat(1,layer{:});
    xyz_tp1{t} = cat(1,layerp1{:});
    XYZ_t{t} = cat(1,poslayer_t{:});
    XYZ_tp1{t} = cat(1,poslayer_tp1{:});
end
%%
trmax=cellfun(@max,XYZ_t(latticeZRange(1):end),'UniformOutput',false);
trmin=cellfun(@min,XYZ_t(latticeZRange(1):end),'UniformOutput',false);
if length(scopeparams)<2
    imsize_um = scopeparams.imsize_um;
else
    imsize_um = scopeparams(1).imsize_um;
end

Rmax = max(cat(1,trmax{:}))+imsize_um*1e3.*[1 1 0];
Rmin = min(cat(1,trmin{:}))-imsize_um*1e3.*[1 1 0];
% plot desctiptors
myfig = 100
figure(myfig), cla, clf
hold on
loops = latticeZRange(end)-latticeZRange(1);
F(loops) = struct('cdata',[],'colormap',[]);
iter=1;
for t = latticeZRange(1:end-1)'
    %%
%     idxtest = sliceinds(375)
    ix = (scopeloc.gridix(:,3)'==t);
    sliceinds = find(ix);
    if any(sliceinds)
    end
    if t>size(XYZ_t,2)
        F(iter) = getframe;
        iter=iter+1;    
        continue
    end
    %%
    figure(myfig), cla, clf, hold on
    view(0,90)
    set(gca,'Color',[1 1 1]*.8)
    if isempty( XYZ_t{t})
        continue
    end
        
    X = XYZ_t{t}(:,1);
    Y = XYZ_t{t}(:,2);
    Z = XYZ_t{t}(:,3);
    U = XYZ_tp1{t}(:,1)-XYZ_t{t}(:,1);
    V = XYZ_tp1{t}(:,1)-XYZ_t{t}(:,2);
    W = XYZ_tp1{t}(:,1)-XYZ_t{t}(:,3);
    % plot boxes
    disp(['    Layer ' num2str(t) ' of ' num2str(max(scopeloc.gridix(:,3)))]);
%     quiver3(X,Y,Z,U,V,W)
    myplot3(XYZ_t{t},'.'),
    %myplot3(XYZ_tp1{t},'+'),
%     legend('x','y')
    x = scopeloc.loc(ix,1)*1e6;
    y = scopeloc.loc(ix,2)*1e6;
    w = imsize_um(1)*ones(sum(ix),1)*1e3;
    h = imsize_um(2)*ones(sum(ix),1)*1e3;
    for ii=1:sum(ix)
        rectangle('Position', [x(ii) y(ii) w(ii) h(ii)])
        mystr = sprintf('%05d\n(%d)\n%d:%d',sliceinds(ii),ii,scopeloc.gridix(sliceinds(ii),1),scopeloc.gridix(sliceinds(ii),2));
        text(x(ii)+w(ii)/2,y(ii)+h(ii)/2,mystr,'Color','r','HorizontalAlignment','center')
    end
    set(gca,'Ydir','reverse')
    xlim([Rmin(1) Rmax(1)])
    ylim([Rmin(2) Rmax(2)])
    title([num2str(t),' - ', num2str(max(D(ix)))])
    text(Rmin(1)+7e5,Rmax(2)-5e5,num2str(t),...
        'FontSize',60,'Color','k','HorizontalAlignment','center')    
    drawnow
    F(iter) = getframe;
    iter=iter+1;
    
end
%%
v = VideoWriter(videofile,'Motion JPEG AVI');
% v.CompressionRatio = 3;
v.FrameRate=5;
open(v)
writeVideo(v,F)
close(v)
%%
end










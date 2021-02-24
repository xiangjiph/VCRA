function [X,sdisp] = createDataMatrix(params,curvemodel,paireddescriptor,scopeloc,idxinlayer)
%%
% fix field curvature, then build a system of linear equations
Ntiles_inlayer = length(idxinlayer);
dims = params.imagesize;
order = params.order; % power to weight shift
imsize_um = params.imsize_um;
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];
flagfc = 1;
X = [];
sdisp=[];
for ic = 1:Ntiles_inlayer
    %%
    % on x Ac*Xc-An*Xn=D
    ic_c = idxinlayer(ic);
    paireddescriptor_ic = paireddescriptor{ic_c};
    ic_n = paireddescriptor_ic.neigs(2);
    icx_inlayer = find(ic_n==idxinlayer);
    Xc = paireddescriptor{ic_c}.onx.X;
    Xn = paireddescriptor{ic_c}.onx.Y;
    sampleSizeforTile = size(Xc,1);
    %%
    if isnan(ic_n) | ~sampleSizeforTile
        continue
    end
    %%
    % fix FC
    subXc = Xc+1;
    subXn = Xn+1;
    [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,ic_c),order,xy,dims,subXc);
    fcXc = locs-1;
    [locs,xshift2D,yshift2D] = util.fcshift(curvemodel(:,:,ic_n),order,xy,dims,subXn);
    fcXn = locs-1;
    %%
    if ~sampleSizeforTile;continue;end
    %%
    sdisp{ic} = [1000*(scopeloc.loc(ic_c,:)-scopeloc.loc(ic_n,:))]'*ones(1,sampleSizeforTile);
    Xsamp = zeros(3,Ntiles_inlayer); % append [0;0;0] for z+1 match
    Xtile = [];
    for ii=1:sampleSizeforTile
        if flagfc
            Xsamp(:,ic) = fcXc(ii,:);
            Xsamp(:,icx_inlayer) = -fcXn(ii,:);
        else
            Xsamp(:,ic) = Xc(ii,:);
            Xsamp(:,icx_inlayer) = -Xn(ii,:);
        end
        Xtile{ii} = Xsamp(:);
    end
    X{ic} = cat(2,Xtile{:});
end
xX = cat(2,X{:});
xsdisp = cat(2,sdisp{:});
%%
X = [];
sdisp=[];
for ic = 1:Ntiles_inlayer
    %%
    % on x Ac*Xc-An*Xn=D
    ic_c = idxinlayer(ic);
    paireddescriptor_ic = paireddescriptor{ic_c};
    ic_n = paireddescriptor_ic.neigs(3);
    if isnan(ic_n)|~paireddescriptor{ic_c}.count(2);continue;end
    icx_inlayer = find(ic_n==idxinlayer);
    Xc = paireddescriptor{ic_c}.ony.X;
    Xn = paireddescriptor{ic_c}.ony.Y;
    sampleSizeforTile = size(Xc,1);
    if ~sampleSizeforTile;continue;end
    sdisp{ic} = [1000*(scopeloc.loc(ic_c,:)-scopeloc.loc(ic_n,:))]'*ones(1,sampleSizeforTile);
    Xsamp = zeros(3,Ntiles_inlayer);
    Xtile = [];
    for ii=1:sampleSizeforTile
        Xsamp(:,ic) = Xc(ii,:);
        Xsamp(:,icx_inlayer) = -Xn(ii,:);
        Xtile{ii} = Xsamp(:);
    end
    X{ic} = cat(2,Xtile{:});
end
yX = cat(2,X{:});
ysdisp = cat(2,sdisp{:});
%%
if 1
    % %% since we dont have bead, just 0/250 pad for z
    numZ = floor((size(xX,2)+size(yX,2))/2/Ntiles_inlayer);
    zX = zeros(3*Ntiles_inlayer,numZ*Ntiles_inlayer);
    for it=1:Ntiles_inlayer
        zX((it-1)*3+[1:3],(it-1)*numZ+[1:numZ]) = [0;0;dims(3)-1]*ones(1,numZ);
    end
    zsdisp = zeros(3,size(zX,2)); zsdisp(end,:) = -(dims(3)-1);
    %%
    X = [xX yX zX];
    sdisp = [xsdisp ysdisp zsdisp];
else
    X = [xX yX];
    sdisp = [xsdisp ysdisp];
end

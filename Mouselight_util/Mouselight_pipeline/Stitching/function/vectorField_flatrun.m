function [vecfield] = vectorField_flatrun(params,scopeloc,scopeparams,affcase)

%VECTORFIELD Summary of this function goes here
%
% [OUTPUTARGS] = VECTORFIELD(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/10/04 12:22:05 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%%
disp('Calculating vector fields');
if nargin<4
    affcase = 2;
end
noopt = 0;
htop = params.htop;
Nlayer = params.Nlayer;
Npts = (params.Ndivs+1).^2;
dims = params.imagesize;
%%
[tileneighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% tileneighbors = tileneighbors(:,checkthese);
%% % initialize points
numTiles = length(scopeloc.filepath);
[corrctrlpnttmp,xlim_cntrl,ylim_cntrl] = util.getCtrlPts(params.imagesize(1),params.imagesize(2),params);
% zlimdefaults = [5 25 dims(3)-26 dims(3)-6];
zlimdefaults = [2 20 dims(3)-21 dims(3)-3];
subcorrctrlpnttmp = corrctrlpnttmp+1;
idxctrl = sub2ind(dims([2 1]),subcorrctrlpnttmp(:,2),subcorrctrlpnttmp(:,1));
mflip = -eye(3);%mflip(3,3)=1;
scopeparams.xshift = zeros(params.imagesize(1:2));
scopeparams.yshift = zeros(params.imagesize(1:2));

switch affcase
    case 0
        afftransform = diag([333 333 1000]);
        mflip=-eye(3);
        mflip(3,3)=1;
    case 1 % affine, no FC
        afftransform = scopeparams.affine;
        scopeparams.xshift(:)=0;
        scopeparams.yshift(:)=0;
    case 2 % affineFC
        afftransform = scopeparams.affineglFC;
    case 3 % affinexy
        afftransform = scopeparams.affinexy;
        scopeparams.xshift(:)=0;
        scopeparams.yshift(:)=0;
    case 4 % affinexyFC
        afftransform = scopeparams.affinexyFC;
        afftransform(1,2)=afftransform(1,2)/2;
        afftransform(2,1)=afftransform(2,1)/2;
    case 5 % affineFC no y
        afftransform = scopeparams.affineFC;
        scopeparams.yshift(:)=0;
    case 6 % affinexyFC, no y
        afftransform = scopeparams.affinexyFC;
        scopeparams.xshift(:)=0;
    case 7 % affine w/o FC, then FC
        afftransform = scopeparams.affine;
        mflip(3,3)=1;
    case 8 % affine w/o FC, then FC
        afftransform= 1e3*repmat([[379.62286 673.389 240]./([1024 1920 241]-1)],3,1)*eye(3).*eye(3)
        mflip=-eye(3);
        mflip(3,3)=1;
    otherwise
end
%%
% flip x/y dimensions to match imaging
% mflip = eye(3);
afftransform = afftransform*mflip;
xshift = scopeparams.xshift(idxctrl);
yshift = scopeparams.yshift(idxctrl);
corrctrlpnttmp(:,1) = corrctrlpnttmp(:,1) + xshift; % correctedcontrolpointtemplate, based on 0 index
corrctrlpnttmp(:,2) = corrctrlpnttmp(:,2) + yshift;
%%
control = zeros(Npts*Nlayer, 3, numTiles);
topcntrl1 = zeros(Npts, 3, numTiles);
topcntrl2 = zeros(Npts, 3, numTiles);
botcntrl1 = zeros(Npts, 3, numTiles);
botcntrl2 = zeros(Npts, 3, numTiles);
zlim_cntrl = zeros(Nlayer, numTiles);
afftile = zeros(3, 4, numTiles);
pixstats = nan(numTiles,8);

filep = strsplit(scopeloc.filepath{1},'/');
vecfield.root = fullfile('/',filep{1:end-4});

for idxt = 1:numTiles
    %%
    % copy path
    filepath = fileparts(scopeloc.filepath{idxt});
    vecfield.path{idxt} = filepath(length(vecfield.root)+1:end);
    % form an affine matrix
    tform = [afftransform scopeloc.loc(idxt,:)'*1e6];
    afftile(:,:,idxt) = tform;
    % top layer at z=15: as we have some noisy data at the top
    zlim_1 = zlimdefaults(1);
    corrctrlpnttmp(:,3)= zlim_1;
    contr_t_top1 = [corrctrlpnttmp ones(Npts,1)]*tform';
    
    % second layer at z=25
    zlim_2 = zlimdefaults(2);
    corrctrlpnttmp(:,3)= zlim_2;
    contr_t_top2 = [corrctrlpnttmp ones(Npts,1)]*tform';
    
    % third layer at z=251-26
    zlim_3 = zlimdefaults(3);
    corrctrlpnttmp(:,3)= zlim_3;
    contr_t_bot1 = [corrctrlpnttmp ones(Npts,1)]*tform';
    
    % forth layer at z=251-16
    zlim_4 = zlimdefaults(4);
    corrctrlpnttmp(:,3)= zlim_4;
    contr_t_bot2 = [corrctrlpnttmp ones(Npts,1)]*tform';
    
    topcntrl1(:,:,idxt) = contr_t_top1;
    topcntrl2(:,:,idxt) = contr_t_top2;
    botcntrl1(:,:,idxt) = contr_t_bot1;
    botcntrl2(:,:,idxt) = contr_t_bot2;
    
    control(1:2*Npts,:,idxt) = [contr_t_top1;contr_t_top2];
    control(2*Npts+1:end,:,idxt) = [contr_t_bot1;contr_t_bot2];
    zlim_cntrl(:,idxt) = [zlim_1 zlim_2 zlim_3 zlim_4];
end

bbox = zeros(8,3,numTiles);
origin = zeros(numTiles,3);
sz = zeros(numTiles,3);
for i = 1:numTiles
    [bbox(:,:,i), origin(i,:), sz(i,:)] = util.bboxFromCorners(control(:,:,i));
end
%% update affines
disp('Update affines based on control points')
tform = zeros(5,5,numTiles);
numX = size(control(:,:,1),1);
for ii=1:numTiles
    %%
    X = control(:,:,ii)'/1000;
    x = xlim_cntrl;
    y = ylim_cntrl;
    z = zlim_cntrl(:,ii)';
    [xx,yy,zz] = ndgrid(x,y,z);
    x = [xx(:),yy(:),zz(:)]';
    X = [X;ones(1,numX)];
    x = [x;ones(1,numX)];
    %A = (X*x')/(x*x');
    A = X/x;
    Aest = eye(5);
    Aest(1:3,1:3) = A(1:3,1:3)*1000;
    Aest(1:3,5) = A(1:3,4)*1000;
    Aest(5,1:3) = A(4,1:3)*1000;
    tform(:,:,ii) = Aest;
end
%%
vecfield.control = control;
vecfield.xlim_cntrl = xlim_cntrl;
vecfield.ylim_cntrl = ylim_cntrl;
vecfield.zlim_cntrl = zlim_cntrl;
vecfield.afftile = afftile;
vecfield.tform = tform;
vecfield.corrctrlpnttmp = corrctrlpnttmp;
vecfield.ctrlpnttmp = subcorrctrlpnttmp-1;
vecfield.bbox = bbox;
vecfield.origin = origin;
vecfield.sz = sz;
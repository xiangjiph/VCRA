function [paireddescriptor,R,curvemodel] = xymatch(descriptors,neigs,scopeloc,params)
%ESTIMATESCOPEPARAMETERS Summary of this function goes here
%
% [OUTPUTARGS] = ESTIMATESCOPEPARAMETERS(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/09/12 10:38:28 $	$Revision: 0.1 $
% Copyright: HHMI 2016
%%
addpath(genpath('./thirdparty'))
debug = 0;
res = 0; 
viz=0;
fignum = 101;

projectionThr = 20; % distance between target and projected point has to be less than this number
dims = params.imagesize;
imsize_um = params.imsize_um;
% slid = [[75 960];[0 dims(2)];[0 dims(3)]];
% expensionshift = [0 0 20]; % HEURISTICS:: tissue expends, so overlap is bigger between tiles

%%
model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
% opt.method='nonrigid_lowrank';
opt.method='nonrigid';
opt.beta=6;            % the width of Gaussian kernel (smoothness), higher numbers make transformation more stiff
opt.lambda=16;          % regularization weight
opt.viz=0;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
opt.fgt=0;              % do not use FGT (default)
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
%     opt.max_it=100;         % max number of iterations
%     opt.tol=1e-10;          % tolerance
matchparams.model = model;
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;
matchparams.debug = debug;
matchparams.viz = viz;
matchparams.fignum = fignum;
matchparams.opt.beta=2;
% matchparams.opt.method = 'nonrigid';
matchparams.init(1,:)=[765 1e-5 867];
matchparams.init(2,:)=[537 -1e-5 1445];

%%
% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% indicies are 1 based,e.g. x = 1:dims(1), not 0:dims(1)-1
% xyz_umperpix = zeros(size(neigs,1),3);
curvemodel = zeros(3,3,size(neigs,1));
R = zeros(3,3,size(neigs,1));
paireddescriptor = cell(size(neigs,1),1);
% initialize
for ix = 1:size(neigs,1)
    paireddescriptor{ix}.onx.valid = 0;
    paireddescriptor{ix}.onx.X = [];
    paireddescriptor{ix}.onx.Y = [];
    paireddescriptor{ix}.ony.valid = 0;
    paireddescriptor{ix}.ony.X = [];
    paireddescriptor{ix}.ony.Y = []; 
    paireddescriptor{ix}.neigs = neigs(ix,:);
    paireddescriptor{ix}.count = [0 0];
end

%
paireddesctemp=[];
paireddesctemp{1}.valid = 0;
paireddesctemp{1}.X = [];
paireddesctemp{1}.Y = [];
paireddesctemp{2}.valid = 0;
paireddesctemp{2}.X = [];
paireddesctemp{2}.Y = [];
%
Ntiles = size(neigs,1);
%%
try parfor_progress(0);catch;end
parfor_progress(Ntiles)

parfor ineig = 1:Ntiles%5163%[5162 5163 5164]%Ntiles%find(neigs(:,1)==5463)%1:size(neigs,1)%5162
    %% load descriptor pairs X (center) - Y (adjacent tile)
    idxcent = neigs(ineig,1);
    descent = descriptors{idxcent};
    if isempty(descent);continue;end
    descent = double(descent(:,1:3));
    if size(descent,1)<3;continue;end

    descent = util.correctTiles(descent,dims); % flip dimensions
    mout = zeros(3,3);
    paireddescriptor_ = paireddesctemp;
    R_ = zeros(3);
    %%
    for iadj = 1:size(neigs,2)-2 %1:x-overlap, 2:y-overlap, 3:z-overlap
        %%
        % idaj : 1=right(+x), 2=bottom(+y), 3=below(+z)
        idxadj =  neigs(ineig,iadj+1);
        
        if isnan(idxadj);continue;end
        descadj = descriptors{idxadj};
        if isempty(descadj);continue;end

        descadj = double(descadj(:,1:3)); % descadj has x-y-z-w1-w2 format
        if size(descadj,1)<3;continue;end
        
        descadj = util.correctTiles(descadj,dims); % flip dimensions
        stgshift = 1000*(scopeloc.loc(idxadj,:)-scopeloc.loc(idxcent,:));
        pixshift = round(stgshift.*(dims-1)./(imsize_um));
        descadj = descadj + ones(size(descadj,1),1)*pixshift; % shift with initial guess based on stage coordinate

        %%
        nbound = [0 0];
        nbound(1) = max(pixshift(iadj),min(descadj(:,iadj)));
        nbound(2) = min(dims(iadj),max(descent(:,iadj)))+0;
        X = descent(descent(:,iadj)>nbound(1)&descent(:,iadj)<nbound(2),:);
        Y = descadj(descadj(:,iadj)>nbound(1)&descadj(:,iadj)<nbound(2),:);
        if size(X,1)<3 | size(Y,1)<3;continue;end
        
        % get descpair
        [X_,Y_] = match.descriptorMatch(X,Y,matchparams);
        if size(X_,1)<3 | size(Y_,1)<3;continue;end
        Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP

        % get field curvature model
        [X_,Y_,out,valid] = match.fcestimate(X_,Y_,iadj,matchparams);
        % flip back dimensions
        X_ = util.correctTiles(X_,dims);
        Y_ = util.correctTiles(Y_,dims);
        
        % store pairs
        mout(iadj,:) = out;
        paireddescriptor_{iadj}.valid = valid;
        paireddescriptor_{iadj}.X = X_;
        paireddescriptor_{iadj}.Y = Y_;
        %R(:,iadj,ineig) = round(median(X_-Y_));
        R_(:,iadj) = round(median(X_-Y_));
    end
    R(:,:,ineig) = R_;
    curvemodel(:,:,ineig) = mout;

    paireddescriptor{ineig}.onx.valid = paireddescriptor_{1}.valid;
    paireddescriptor{ineig}.onx.X = paireddescriptor_{1}.X;
    paireddescriptor{ineig}.onx.Y = paireddescriptor_{1}.Y;
    
    paireddescriptor{ineig}.ony.valid = paireddescriptor_{2}.valid;
    paireddescriptor{ineig}.ony.X = paireddescriptor_{2}.X;
    paireddescriptor{ineig}.ony.Y = paireddescriptor_{2}.Y;
    paireddescriptor{ineig}.count = [size(paireddescriptor_{1}.X,1) size(paireddescriptor_{2}.X,1)];
    parfor_progress;
end
parfor_progress(0);
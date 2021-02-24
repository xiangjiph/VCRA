function [X_,Y_,out] = descriptorMatch(X,Y,pixshift,iadj,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/09/23 14:09:29 $	$Revision: 0.1 $
% Copyright: HHMI 2016

opt = params.opt;
model = params.model;
optimopts = params.optimopts;
projectionThr = params.projectionThr;
debug = params.viz;
fignum = params.fignum;
%%
% initial match based on point drift
[Transform, C]=cpd_register(X,Y,opt);
%% check if match is found
pD = pdist2(X,Transform.Y);
[aa1,bb1]=min(pD,[],1);
[aa2,bb2]=min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)'<projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
if isempty(X_) | isempty(Y_)
    out = zeros(1,3);
    return
end
Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
%%
% displacement field between follows a field curve on x&y due to
% optics and deformation curve due to tissue and cut force on z
dispvec = X_-Y_;
x = dispvec(:,iadj);
% for non focus axis reject outliers based on vector norm. This should be
% (roughly) constant for non curvature directions
vcomp = dispvec(:,setdiff(1:3,iadj));
medvcomp = median(vcomp);
normvcomp = vcomp-ones(size(vcomp,1),1)*medvcomp;
normvcomp = sqrt(sum(normvcomp.*normvcomp,2));

validinds = 1:length(normvcomp);
X_ = X_(validinds,:);
Y_ = Y_(validinds,:);
x = x(validinds,:);

if iadj==1 % x-neighbor
    % x : x-displacement
    % y : y-location
    y = X_(:,2);
    bw = [2 220];
elseif iadj==2 % y-neighbor
    % x : y-displacement
    % y : x-location
    y = X_(:,1);
    bw=[3 50];
else % z-neighbor
    % x : z-displacement
    % y : y-location (not too much on x as cut is on y direction)
    y = X_(:,2);
    bw=[2 220];
end
% build a probabilistic model of displacement vectors
N = 101;
gridx = linspace(min(x),max(x),N);
gridy = linspace(min(y),max(y),N);
[density,bw] = util.ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));

if debug
    %%
    figure(7)
    cla
    imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
    axis tight
    hold on,
    plot(x,y,'m.')
end
%%
[xmin,ix] = min(pdist2(x,gridx'),[],2);
[ymin,iy] = min(pdist2(y,gridy'),[],2);
idx = sub2ind([N,N],iy,ix);
prob_inliers = density(idx)>max(density(idx))*.25;
x_inline = x(prob_inliers,:);
y_inline = y(prob_inliers,:);
%
% fit curve model
[~,im] = max(density,[],iadj);
sgn = 2*((max(im)==im(end) | max(im)==im(1))-.5);
if isfield(params,'init')
    pinit = params.init(iadj,:);
else
    pinit = [median(y) sgn*1e-5 median(x)];
end
% pinit = [800 sgn*1e-5 median(x)];
% turn off warning
warning off
out = nlinfit(y_inline, x_inline, model, pinit,optimopts);
%%
if abs(out(2))<sqrt(eps) % mostlikely a line, p(1) is not reliable
    out(2) = 0; % to prevent scaling error in fc images
end
% if sign(pinit(2))~=out(2) | abs(pinit(1)-out(1))>40
%     % try with reverse sign
%     out = nlinfit(y_inline, x_inline, model, pinit.*[1 -1 1],optimopts);
% 
%     BW = imregionalmax(density);
%     [aa,bb,~] = find(BW);
%     dens = 0*aa;
%     for ii = 1:length(aa(:))
%         dens(ii) = density(aa(ii),bb(ii));
%     end
%     [gridy(aa);gridx(bb)]
% turn on warning
warning on

% outlier rejection based on parametric model
xest = feval(model,out,y);
outliers = abs(x-xest)>2;
X_ = X_(~outliers,:);
Y_ = Y_(~outliers,:);
if debug
%%
    tts{1} = [1 4];
    tts{2} = [2 3];
    tts{3} = [7];
    xgridest = feval(model,out,gridy);
    figure(fignum),
    subplot(3,3,tts{iadj}),cla
    imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
    axis tight
    hold on,
    plot(x,y,'m.')
    plot(x_inline,y_inline,'mo')
    plot(x(~outliers),y(~outliers),'gd')
    plot(xgridest,gridy,'r-')
    
    subplot(3,3,[5 6 8 9]),
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'b+')
    plot3(Y(:,1),Y(:,2),Y(:,3),'r.')
    plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'ro')
    pause(.5)
    axis equal
    drawnow
    
end
end

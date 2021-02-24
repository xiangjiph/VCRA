function [rate,X_,Y_,tY_] = vessel_descriptorMatchforz(X,Y,pixshift,iadj,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
%
% Inputs:
%   X, Y: two 2D real, double marices, specifying the position of the point
%   in the point cloud.
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
out = [];
opt = params.opt;
model = params.model;
projectionThr = params.projectionThr;
debug = params.viz;
%% Initial match based on point drift
[Transform, C] = cpd_register(X,Y,opt);
%% check if match is found
% Compute the pairwise euclidean distance between two input array
pD = pdist2(X,Transform.Y);
[aa1,bb1] = min(pD,[],1);
[aa2,bb2] = min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)' < projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
tY_= Transform.Y(keeptheseY(disttrim),:);
% Rate is the ratio of the number of matched pair of distance less than
% projectionThr over the total number of matched pairs
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
% if rate < .5 % dont need to continue
%     [X_,Y_,out] = deal(0);
%     return
% end
%%
% Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
Y_ = bsxfun(@minus, Y_, pixshift);
% Y_ = Y_ - ones(size(Y_,1),1)*pixshift;% move it back to original location after CDP

% %%
% % displacement field between follows a field curve on x&y due to
% % optics and deformation curve due to tissue and cut force on z
% dispvec = X_-Y_;
% x = dispvec(:,iadj);
% if iadj==1 % x-neighbor
%     % x : x-displacement
%     % y : y-location
%     y = X_(:,2);
%     bw = [2 220];
% elseif iadj==2 % y-neighbor
%     % x : y-displacement
%     % y : x-location
%     y = X_(:,1);
%     bw=[3 100];
% else % z-neighbor
%     % x : z-displacement
%     % y : y-location (not too much on x as cut is on y direction)
%     y = X_(:,2);
%     bw=[2 220];
% end
% % build a probabilistic model of displacement vectors
% N = 101;
% gridx = linspace(min(x),max(x),N);
% gridy = linspace(min(y),max(y),N);
% [density,bw] = ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));
% [xmin,ix] = min(pdist2(x,gridx'),[],2);
% [ymin,iy] = min(pdist2(y,gridy'),[],2);
% idx = sub2ind([N,N],iy,ix);
% prob_inliers = density(idx)>max(density(idx))*.25;
% x_inline = x(prob_inliers,:);
% y_inline = y(prob_inliers,:);
% %
% % fit curve model
% [~,im] = max(density,[],2);
% sgn = -2*((max(im)==im(end) | max(im)==im(1))-.5);
% pinit = [median(y) sgn*1e-5 median(x)];
% warning off
% out = nlinfit(y_inline, x_inline, model, pinit,optimopts);
% warning on
% % outlier rejection based on parametric model
% xest = feval(model,out,y);
% outliers = abs(x-xest)>2;
% X_ = X_(~outliers,:);
% Y_ = Y_(~outliers,:);
% tY_ = tY_(~outliers,:);
% if debug
%     xgridest = feval(model,out,gridy);
%     figure(100),
%     subplot(2,2,iadj),cla
%     imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
%     axis tight
%     hold on,
%     plot(x,y,'m.')
%     plot(x_inline,y_inline,'mo')
%     plot(x(~outliers),y(~outliers),'gd')
%     plot(xgridest,gridy,'r-')
%     
%     subplot(2,2,4),
%     cla
%     hold on
%     plot3(X(:,1),X(:,2),X(:,3),'b+')
%     plot3(Y(:,1),Y(:,2),Y(:,3),'m.')
%     plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'ro')
%     pause(.5)
%     drawnow
%     
% end
end
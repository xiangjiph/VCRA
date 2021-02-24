function [X_,Y_,rate_,pixshift,nonuniformity] = searchpair(descent,descadjori,pixshiftinit,iadj,dims,matchparams)
%SEACHPAIR Summary of this function goes here
%
% [OUTPUTARGS] = SEACHPAIR(INPUTARGS) Explain usage here
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

% $Author: base $	$Date: 2016/11/03 16:02:56 $	$Revision: 0.1 $
% Copyright: HHMI 2016

pixshift = pixshiftinit;
% pixshift = zeros(1,size(descent,2));
% pixshift(1:length(pixshiftinit)) = pixshiftinit;
%search
flag = 0;
iter = 0;
R = zeros(1,10);
% nonuniformity = zeros(1,10);
clear nonuniformity

[X_,Y_,neigs_,rate_] = deal([]);
while ~flag & iter<50% run a search
    %%
    iter = iter + 1;
    descadj = descadjori(:,1:3) + ones(size(descadjori,1),1)*pixshift ;
    
    nbound = [0 0];
    nbound(1) = max(pixshift(iadj),min(descadj(:,iadj)));
    nbound(2) = min(dims(iadj),max(descent(:,iadj)))+3;
    X = descent(descent(:,iadj)>nbound(1)&descent(:,iadj)<nbound(2),:);
    Y = descadj(descadj(:,iadj)>nbound(1)&descadj(:,iadj)<nbound(2),:);
    X = X(:,1:3);
    Y = Y(:,1:3);

    %%
    if size(X,1)<3 | size(Y,1)<3% not enough sample to match
        [X_,Y_,rate_,pixshift,nonuniformity] = deal([]);
        flag = 1;
    else
        %%
        % check uniformity of data
        nbins = [2 2];
        edges = [];
        for ii=1:2%length(dims)%[1 2 3],
            minx = 0;
            maxx = dims(ii);
            binwidth = (maxx - minx) / nbins(ii);
            edges{ii} = minx + binwidth*(0:nbins(ii));
        end
        [accArr] = hist3([X(:,1:2);Y(:,1:2)],'Edges',edges);
        accArr = accArr(1:2,1:2);
        if ~all(sum(accArr>mean(accArr(:))) & sum(accArr>mean(accArr(:)),2)')
            % non uniform over quad-representation
            nonuniformity(iter) = 1;
        else
            nonuniformity(iter) = 0;
        end
        
        try
            %%
            [rate,X_,Y_,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,matchparams);
            if size(X_,1)<3
                rate = 0; % overparametrized system
            end
            R(iter) = rate;
            if iter>1 & R(iter)-R(iter-1)<0
                flag = 1;
                X_ = X_t_1;
                Y_ = Y_t_1;
                rate = R(iter-1);
            else
                X_t_1 = X_;
                Y_t_1 = Y_;
                if rate<.95 & iadj ==3% no match
                    pixshift = pixshift + [0 0 5]; % expand more
                    flag = 0;
                    error('increase shift')
                else % match found
                    flag = 1;
                end
            end
            % store pairs
            rate_ = rate;
        catch
            X_ = [];
            Y_ = [];
            disp('error')
        end
    end
end
% [iter R(end)]

%%

end
function [rate,X_,Y_,tY_] = descriptorMatchforz_relaxed(X,Y,pixshift,iadj,params)
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
tY_ = [];
opt = params.opt;
model = params.model;
optimopts = params.optimopts;
projectionThr = params.projectionThr;
debug = params.viz;
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
tY_= Transform.Y(keeptheseY(disttrim),:);
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
if rate < .5 % dont need to continue
    [X_,Y_,out] = deal(0);
    return
end
%%
% Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
Y_ = Y_- ones(size(Y_,1),1)*pixshift;% move it back to original location after CDP

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
if debug
    xgridest = feval(model,out,gridy);
    figure(100),
    subplot(2,2,iadj),cla
    imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
    axis tight
    hold on,
    plot(x,y,'m.')
    plot(x_inline,y_inline,'mo')
    plot(x(~outliers),y(~outliers),'gd')
    plot(xgridest,gridy,'r-')
    
    subplot(2,2,4),
    cla
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'b+')
    plot3(Y(:,1),Y(:,2),Y(:,3),'m.')
    plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'ro')
    pause(.5)
    drawnow
    
end
end

function [rate,X_,Y_,tY_] = descriptorMatchforz(X,Y,pixshift,iadj,params)
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
tY_ = [];
out = [];
opt = params.opt;
model = params.model;
optimopts = params.optimopts;
projectionThr = params.projectionThr;
debug = params.viz;
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
tY_= Transform.Y(keeptheseY(disttrim),:);
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
if rate < .5 % dont need to continue
    [X_,Y_,out] = deal(0);
    return
end
%%
% Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
Y_ = Y_- ones(size(Y_,1),1)*pixshift;% move it back to original location after CDP

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
if debug
    xgridest = feval(model,out,gridy);
    figure(100),
    subplot(2,2,iadj),cla
    imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
    axis tight
    hold on,
    plot(x,y,'m.')
    plot(x_inline,y_inline,'mo')
    plot(x(~outliers),y(~outliers),'gd')
    plot(xgridest,gridy,'r-')
    
    subplot(2,2,4),
    cla
    hold on
    plot3(X(:,1),X(:,2),X(:,3),'b+')
    plot3(Y(:,1),Y(:,2),Y(:,3),'m.')
    plot3(Transform.Y(:,1),Transform.Y(:,2),Transform.Y(:,3),'ro')
    pause(.5)
    drawnow
    
end
end


function [X_,Y_,out,valid] =  fcestimate(X_,Y_,iadj,params)
%FCESTIMATE Summary of this function goes here
%   Detailed explanation goes here
valid = 0;
viz = 1| params.viz;
% model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
model = params.model;
optimopts = params.optimopts;
if isfield(params,'dims')
    dims = params.dims;
else
    dims = [1024 1536 251];
end
dispvec = X_-Y_;
y = dispvec(:,iadj);

% for non focus axis reject outliers based on vector norm. This should be
% (roughly) constant for non curvature directions
if 1
    validinds = 1:length(y);
elseif length(y)>20 
    vcomp = dispvec(:,setdiff(1:3,iadj));
    medvcomp = median(vcomp);
    normvcomp = vcomp-ones(size(vcomp,1),1)*medvcomp;
    normvcomp = sqrt(sum(normvcomp.*normvcomp,2));
    validinds = normvcomp<util.get1DThresh(normvcomp,20,.95);
end

X_ = X_(validinds,:);
Y_ = Y_(validinds,:);
y = y(validinds,:);

x = X_(:,setdiff([1 2],iadj));
dimcent = dims(setdiff([1:2],iadj))/2; % center of image along curvature axis
% polynomial coeeficients (p3-p2(y-p1)^2):
% p(1) : imaging center ~ dims/2 +/- %10
% p(2) : curvature: -/+[1e-5 1e-6] for x & y, use initialization if
% avaliable, curvature might flip sign based on objective or reduce to 0 as
% medium changes (due to temperature or adding liquid solution, etc.)
% p(3): avarage displacement: between [[1-%overlap]*dims dims],
% initialization might not be useful, as this reduces to mean descriptor
% displacement
if isfield(params,'init')
    pinit = params.init(iadj,:);
    pinit(3) = median(y); 
else
    pinit = [dimcent -dimcent^-2 median(y)]; 
end
% set upper and lower boundaries
% imaging center is around image center
lb1 = dimcent-dimcent*0.1;
ub1 = dimcent+dimcent*0.1;
% curvature should rely on initialization as magnitude and sign might change
% percent ratios do not make sense here as this number get squared, so
% provide a large range
lb2 = -1e-4; % (dims/2).^2*lb2 ~ [25 60] pixels in x & y
ub2 = 1e-4; 
% mean displacement is initialized based on descriptors
lb3 = pinit(3)-pinit(3)*0.1;
ub3 = pinit(3)+pinit(3)*0.1;

ub = [ub1 ub2 ub3];
lb = [lb1 lb2 lb3];
%%
if 1
    fun = @(p) sum((y-feval(model,p,x)).^2);
    sqerr = @(p) sum((y-feval(model,p,x)).^2);
    
    options = optimoptions('fmincon','Display', 'off');
    options.ConstraintTolerance = 1e-9;
    options.OptimalityTolerance = 1e-9;
    out = fmincon(fun,pinit,[],[],[],[],lb,ub,[],options);
else
    % turn off warning
    warning off
    %TODO: very inefficient, find a robust way to pick inflection sign
    [out1,r1] = nlinfit(x, y, model, pinit,optimopts);
    [out2,r2] = nlinfit(x, y, model, pinit.*[1 -1 1],optimopts);
    if norm(r1)<norm(r2);out = out1;else out = out2;end
    
    if abs(out(2))<(dims(iadj)/2)^-2 % mostlikely a line, p(1) is not reliable
        out(2) = 0; % to prevent scaling error in fc images
    end
    warning on
end
% [ out.*[1 1e6 1] sqerr(out)]
% modd = params.init(:,:);
% modd(iadj,:) = out;
% match.vizCurvature(modd)

%%
% outlier rejection based on parametric model
yest = feval(model,out,x);
x_range = min(x):max(x);
y_range = feval(model,out,x_range);
outliers = abs(y-yest)>2;
%%
% if percentage of otliers is large, dont do correction!!
if sum(outliers)/length(outliers) < .25
    X_ = X_(~outliers,:);
    Y_ = Y_(~outliers,:);
    valid=1;
else
    out = out*0;
    valid = 0;
end

%%
if viz
    %%
    if iadj==1
        % range
        ran = [min(dispvec);max(dispvec)];xl = ran(:,iadj);xl=xl(:)';
        yl = [0 dims(setdiff([1 2],iadj))];
        

        figure(304),
        subplot(2,2,[1 3])
        cla
        xlim(xl), ylim(yl)
        plot(y,x,'+')
        hold on
        plot(y(outliers),x(outliers),'ro')
        plot(y_range,x_range,'g-')
        plot(feval(model,out,x),x,'go')
        daspect([1 50 1])
        
    elseif iadj==2
        figure(304),
        subplot(2,1,iadj)
        cla
        plot(x,y,'+')
        hold on
        plot(x(outliers),y(outliers),'ro')
        plot(x_range,y_range,'g-')
        plot(x,yest,'go')
        %     daspect([1 100 1])
    end
end
end


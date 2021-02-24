function fc_params =  fun_fcestimate_vessel(tile_1_xyz,tile_2_xyz,iadj,params)
%FCESTIMATE Summary of this function goes here
%   Detailed explanation goes here
viz = params.viz;
% model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
model = params.model;
% optimopts = params.optimopts;
if isfield(params,'dims')
    dims = params.dims;
else
    dims = [1024 1536 251];
end
dispvec = tile_1_xyz - tile_2_xyz;
disp_in_adj = dispvec(:,iadj);
num_points = size(disp_in_adj,1);

% for non focus axis reject outliers based on vector norm. This should be
% (roughly) constant for non curvature directions
if 1
    validinds = 1:length(disp_in_adj);
elseif length(y)>20 
    vcomp = dispvec(:,setdiff(1:3,iadj));
    medvcomp = median(vcomp);
    normvcomp = vcomp-ones(size(vcomp,1),1)*medvcomp;
    normvcomp = sqrt(sum(normvcomp.*normvcomp,2));
    validinds = normvcomp<util.get1DThresh(normvcomp,20,.95);
end

tile_1_xyz = tile_1_xyz(validinds,:);
tile_2_xyz = tile_2_xyz(validinds,:);
disp_in_adj = disp_in_adj(validinds,:);
x_sub_idx = setdiff([1 2],iadj);
x = tile_1_xyz(:,x_sub_idx);
dimcent = dims(setdiff([1:2],iadj))/2; % center of image along curvature axis
% polynomial coeeficients (p3-p2(y-p1)^2):
% p(1) : imaging center ~ dims/2 +/- %10
% p(2) : curvature: -/+[1e-5 1e-6] for x & y, use initialization if
% avaliable, curvature might flip sign based on objective or reduce to 0 as
% medium changes (due to temperature or adding liquid solution, etc.)
% p(3): average displacement: between [[1-%overlap]*dims dims],
% initialization might not be useful, as this reduces to mean descriptor
% displacement
if isfield(params,'init')
    pinit = params.init(iadj,:);
    pinit(3) = median(disp_in_adj); 
else
    pinit = [dimcent -dimcent^-2 median(disp_in_adj)]; 
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
%% Fitting the curve by nonlinear minimization 
if 1
    fun = @(p) sum((disp_in_adj-feval(model,p,x)).^2);
    options = optimoptions('fmincon','Display', 'off');
    options.ConstraintTolerance = 1e-9;
    options.OptimalityTolerance = 1e-9;
    out= fmincon(fun,pinit,[],[],[],[],lb,ub,[],options);
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
%% outlier rejection based on parametric model
yest = feval(model,out,x);
x_range = min(x):max(x);
y_range = feval(model,out,x_range);
max_disp_error = 2.0;

fc_params.model = model;
fc_params.fitting_params = out;
fc_params.x_sub_idx = x_sub_idx;
fc_params.y_sub_idx=  iadj;
fc_params.X_in = tile_1_xyz;
fc_params.Y_in = tile_2_xyz;
% if percentage of otliers is large, dont do correction!!
outliers_Q = abs(disp_in_adj - yest) > max_disp_error; 
fc_params.inliers_ratio = nnz(~outliers_Q) / num_points;
if fc_params.inliers_ratio > 0.75
    fc_params.valid = true;
else
    max_disp_error = max_disp_error + 1;
    outliers_Q = abs(disp_in_adj - yest) > max_disp_error; 
    fc_params.inliers_ratio = nnz(~outliers_Q) / num_points;
    if fc_params.inliers_ratio > 0.9
        fc_params.valid = true;
    else
        fc_params.valid = false;
    end
end
if fc_params.valid
    fc_params.X_in = fc_params.X_in(~outliers_Q,:);
    fc_params.Y_in = fc_params.Y_in(~outliers_Q,:);
end
fc_params.max_disp_error = max_disp_error;
%% Visualization
if viz
    figure;
    if iadj==1
        ran = [min(dispvec);max(dispvec)];xl = ran(:,iadj);xl=xl(:)';
        yl = [0 dims(setdiff([1 2],iadj))];
        subplot(2,2,[1 3])
        cla
        xlim(xl), ylim(yl)
        plot(disp_in_adj,x,'+')
        hold on
        plot(disp_in_adj(outliers_Q),x(outliers_Q),'ro')
        plot(y_range,x_range,'g-')
        plot(feval(model,out,x),x,'go')
        daspect([1 50 1])   
    elseif iadj==2
        subplot(2,1,iadj)
        cla
        plot(x,disp_in_adj,'+')
        hold on
        plot(x(outliers_Q),disp_in_adj(outliers_Q),'ro')
        plot(x_range,y_range,'g-')
        plot(x,yest,'go')
    end
end
end


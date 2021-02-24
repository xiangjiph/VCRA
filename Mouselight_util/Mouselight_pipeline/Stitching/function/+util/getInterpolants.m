function [Fxt, Fyt, Fzt, Fxtp1, Fytp1, Fztp1,XYZ_t_ori,XYZ_tp1_ori,outliers] = ...
    getInterpolants(ix,regpts,afftile,params,curvemodel)
% Note: 
% 
%
%
%% Initialization 
viz = params.viz;
dims = params.imagesize;
cnt = 0;
Npts = 0;
interpolation_method = 'linear';
extrapolation_method = 'nearest';

if isfield(params,'order')
    order = params.order;
else
    order = 1;
end
if isfield(params,'applyFC') % If true, apply field curvature correction 
    applyFC = params.applyFC;
else
    applyFC = 0;
end
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

[poslayer_t,poslayer_tp1,Fxt, Fyt, Fzt, Fxtp1, Fytp1, Fztp1,XYZ_t_ori,XYZ_tp1_ori,outliers] = deal([]);
median_descriptor_shift = nan(3, nnz(ix));
%% Correct the displacement of the descriptor due to field curvature and homography
ind_list = find(ix);
for iter_pair = 1 : numel(ind_list)
    id_ix = ind_list(iter_pair);
    if isempty(regpts{id_ix}.X)
        % If no descriptor, skip.
        continue
    end
    neigs = regpts{id_ix}.neigs;
    % regpts is a sturcture with matched pairs in Z direction. X and Y
    % stand for tile 1 and tile 2 respectively. 
    layer = regpts{id_ix}.X;
    layerp1 = regpts{id_ix}.Y;
    % Median displacement between matched descriptor
    median_descriptor_shift(:, iter_pair) = round(median(layer - layerp1, 1))';
    % apply FC - correct for the field curvature
    if applyFC
        [layer] = util.fcshift(curvemodel(:,:,neigs(1)),order,xy,dims,layer + 1);
        layer = layer - 1;
        [layerp1] = util.fcshift(curvemodel(:,:,neigs(4)),order,xy,dims,layerp1 + 1);
        layerp1 = layerp1 - 1;
    end
    
    % apply affine transformation, using the matrix computed in
    % fun_estimateaffine_vessel - correct for the homography
    layer_t = [layer ones(size(layer,1),1)] * afftile(:,:,neigs(1))';
    layer_tp1 = [layerp1 ones(size(layerp1,1),1)] * afftile(:,:,neigs(4))';
    
    cnt = cnt + 1;
    Npts = Npts + size(layer_t,1);
    poslayer_t{cnt} = layer_t;
    poslayer_tp1{cnt} = layer_tp1;
end
if isempty(poslayer_t),return,end
%% Further reject outlier in the descriptor list 
% 1. For each descriptor, search for the K-nearest neighbor descriptors in
% the tile above, record their list indices: IDX
% 2. Compute the displacement for the matched descriptor pairs: diffXY
% 3. Compute the 2D Euclidean distance between descritpor and its KNN
% descriptors: dists
XYZ_t = cat(1, poslayer_t{:});
XYZ_tp1 = cat(1, poslayer_tp1{:});
% Store the matched descriptor positions before outlier rejection
XYZ_t_ori = XYZ_t;
XYZ_tp1_ori = XYZ_tp1;
if 1
    is_inliers_Q = fun_feature_match_smooth_displacement_field(XYZ_t, XYZ_tp1);
else
    diffZ = XYZ_t(:,3)-XYZ_tp1(:,3);diffZ=diffZ(IDX);
    % compare first to rest
    inliers = abs(diffZ(:,1)) <= 2e3 | abs(diffZ(:,1)) <= 3*abs(median(diffZ(:,2:end),2));
end
%% Visualization 
if viz
    figure(35), cla
    hold on
    plot3(XYZ_t(:,1),XYZ_t(:,2),XYZ_t(:,3),'r.') % layer t
    plot3(XYZ_tp1(:,1),XYZ_tp1(:,2),XYZ_tp1(:,3),'k.') % layer tp1
    %text(XYZ_t(outliers,1),XYZ_t(outliers,2),num2str(st(outliers,2:4)))
    myplot3(XYZ_t(outliers,:),'md') % layer t
    myplot3(XYZ_tp1(outliers,:),'gd') % layer t
    %plot3(XYZ_t(IDX(idx,1),1),XYZ_t(IDX(idx,1),2),XYZ_t(IDX(idx,1),3),'r*') % layer t
    %plot3(XYZ_t(IDX(idx,2:end),1),XYZ_t(IDX(idx,2:end),2),XYZ_t(IDX(idx,2:end),3),'bo') % layer t
    set(gca,'Ydir','reverse')
    drawnow    
end
%% Remove outliers
XYZ_t = XYZ_t(is_inliers_Q,:);
XYZ_tp1 = XYZ_tp1(is_inliers_Q,:);
%% Visualization
if viz
    figure(33), cla
    hold on
    plot3(XYZ_t(:,1),XYZ_t(:,2),XYZ_t(:,3),'k.') % layer t
    plot3(XYZ_tp1(:,1),XYZ_tp1(:,2),XYZ_tp1(:,3),'r.') % layer tp1
    myplot3([layer ones(size(layer,1),1)]*afftile(:,:,neigs(1))','g.') % layer t
    myplot3([layerp1 ones(size(layer,1),1)]*afftile(:,:,neigs(4))','m.') % layer tp1
    set(gca,'Zdir','reverse');
    legend('layer t','layer t+1')
    set(gca,'Ydir','reverse')
    daspect([1,1,1]);
    %         plot3(XYZ_t(inliers,1),XYZ_t(inliers,2),XYZ_t(inliers,3),'go') % layer t
    %         plot3(XYZ_tp1(inliers,1),XYZ_tp1(inliers,2),XYZ_tp1(inliers,3),'mo') % layer tp1
    %         plot3(XYZ_t(16248,1),XYZ_t(16248,2),XYZ_t(16248,3),'g*') % layer t
    %         plot3(XYZ_tp1(16248,1),XYZ_tp1(16248,2),XYZ_tp1(16248,3),'m*') % layer tp1
end
%% Compute deformation field
% displacement vectors of the reliable descriptor paris. 
% Split the deformation into two parts, half for tile 1 and another half for
% tiel 2 ( -1/2 vecdif and +1/2 vecdif below). The resulting is a vector
% field. 
% Notice that this deformation field is the nonlinear local deformation
% field, after field curvature correction and affine transformation. The
% effect accommodated in field curvature correction is unclear to me, but
% the following deformation field certainly encode the nonlinear
% defformation of the tissue. 
% Q: why not compute the overall decormation field directly? 
vecdif = XYZ_t - XYZ_tp1;
rt = params.expensionratio/(1+params.expensionratio);
% layer t
Fxt = scatteredInterpolant(XYZ_t, -vecdif(:,1) * (1-rt), interpolation_method, extrapolation_method);
Fyt = scatteredInterpolant(XYZ_t, -vecdif(:,2) * (1-rt), interpolation_method, extrapolation_method);
Fzt = scatteredInterpolant(XYZ_t, -vecdif(:,3) * (1-rt), interpolation_method, extrapolation_method);
% layer t+1
Fxtp1 = scatteredInterpolant(XYZ_tp1, vecdif(:,1) * rt, interpolation_method, extrapolation_method);
Fytp1 = scatteredInterpolant(XYZ_tp1, vecdif(:,2) * rt, interpolation_method, extrapolation_method);
Fztp1 = scatteredInterpolant(XYZ_tp1, vecdif(:,3) * rt, interpolation_method, extrapolation_method);
end
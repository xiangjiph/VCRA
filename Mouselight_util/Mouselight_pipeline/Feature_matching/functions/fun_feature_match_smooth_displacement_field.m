function is_inliers_Q = fun_feature_match_smooth_displacement_field(XYZ_t, XYZ_tp1)


%% build kd-tree based on xy location for outlier rejection
num_nearest_neighbor = min(20,round(sqrt(size(XYZ_t,1))));
IDX = knnsearch(XYZ_t(:,1:2), XYZ_t(:,1:2), 'K', num_nearest_neighbor);
% Find the K-nearest neighbor in xy direciton, whose displacement vectors
% between two tiles are expected to be consistent. 
% Select inliers in XYZ_t and XYZ_tp1 ( t plus 1)
% diffXYZ = XYZ_t - XYZ_tp1;
% unit_diffXYZ = normr(diffXYZ - median(diffXYZ));
diffXY = XYZ_t(:,1:2) - XYZ_tp1(:,1:2);
diffXY_norm = vecnorm(diffXY, 2, 2);
unit_diffXY = diffXY ./ diffXY_norm;
% interpolate vector from nearest K samples
bins_cos = [linspace(-1,1,21) - 0.05 1.05];

bins_xy_disp = [linspace(0.1,1,10) - 0.05 1.05];
st = zeros(size(diffXY,1),4);
max_dist_2_neighbor = 1e6;% 1e6 nm = 1e3 um = 1 mm
for idx = 1 : size(diffXY, 1)
    % dists is the euclidean distance between point idx and the K nearest
    % neighbors among all the matched voxels in tile 1 in x-y plane
    l2_dist_to_knn = [0; sqrt(sum((ones(num_nearest_neighbor - 1, 1) * XYZ_t(idx,1:2) - XYZ_t(IDX(idx, 2:end), 1:2)) .^ 2 ,2))]; % can be used as weighting
    
    % inner product of the normalized relative displacement on the xy plane
    % and the KNN displacement. i.e. Cos(v_idx, v_KNN)
    inner_prod = [1 unit_diffXY(idx,:) * unit_diffXY(IDX(idx, 2 : end),:)'];
    
    % Consistance in local displacement field on XY plane
    % cos theta
    [~, idx_max_theta] = max(histc(inner_prod(l2_dist_to_knn > 0 & l2_dist_to_knn < max_dist_2_neighbor), bins_cos));
    % idxmaxtheta is the mode cos value
    % st(idx, 2) is the mode of the projection angle between vector idx and
    % its KNN points. This value should be closed to 1 since the local
    % displacement should be consistent.
    st(idx,2) = bins_cos(idx_max_theta) + 0.05;
    
    
    % difference between the XY displacement vector of point idx and the
    % displacement vectores of its KNN points.
    dV = ones(num_nearest_neighbor, 1) * diffXY(IDX(idx,1),:) - diffXY(IDX(idx, :),:);
    
    dV_norm = sqrt(sum(dV .^ 2, 2)); % L2 norm
    dV_norm = exp(-dV_norm/diffXY_norm(idx)); % normalized by the L2 norm of the displacement vector
    assert(dV_norm(1) == 1);
    
    aa = histc(dV_norm(l2_dist_to_knn > 0 & l2_dist_to_knn < max_dist_2_neighbor), bins_xy_disp); % Max-likely
    idx_maxdist = find(aa == max(aa), 1, 'last');
    st(idx,3) = bins_xy_disp(idx_maxdist) + 0.05; % mode normalized deviation in the displacement vector
    
    st(idx,4) = (bins_xy_disp - 0.05) * aa(:) / sum(aa); % average of the normalized deviation, weighted by the count of KNN displacement vector deviation
end
outliers = (st(:,2) < 0.8) | (st(:,3) < 0.5 & st(:,4) < 0.5);
is_inliers_Q = ~outliers;
%% Visualization
% figure;
% quiver3((XYZ_t(:, 1) - min(XYZ_t(:, 1))) ./1000, (XYZ_t(:, 2) - min(XYZ_t(:, 2))) ./1000,...
%     (XYZ_t(:, 3) - min(XYZ_t(:, 3))) ./1000, normalized_diffXYZ(:, 1), ...
%     normalized_diffXYZ(:, 2), normalized_diffXYZ(:, 3))
% daspect([1,1,1]);
end

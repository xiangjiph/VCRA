function lf = fun_graph_get_link_w_1_ep_features(link_cc, vessel_image, vessel_mask_dt)
% fun_graph_get_link_w_1_ep_features computes the features derived from the
% position and orientation of endpoints and their nearest neighbor
% Input: 
%   ep_sub: N-by-3 numerical array, where N is the number of endpoint and 3
%   is the 3D subscript of the endpoint in the array
%   ep_vec: 3-by-N numerical array, the normalize direction vector of the
%   endpoint, which is computed as the first normalized principle component
%   of the link that the endpoint attached to, which is pointing along the link and
%   ends at the endpoint. 
% Output: 
%   lf: table with fields
%       num_voxel: number of voxels in each link connected component
%       length: Euclidean length of the link
%       ep2ep_dist: Euclidean distance between the two endpoints of the
%       link
%       dt_max: maximum distance transform(DT) value of the link voxels
%       dt_min: minimum DT value of the link voxels 
%       dt_mean: mean value of the DT of the link voxels
%       dt_std: standard deviation of the DT of the link voxel
%       dt_diff_ep2ep: absolute value of the difference between the 
%       DT value of the two link endpoints 
%       int_mean: mean value of the intensity of the link voxels
%       int_std: standard deviation of the intensity of the link voxels
%       int_median: median value of the intensity of the link voxels
%       int_diff_ep2ep: absolute value of the difference between the
%       intensity of the two link endpoints
%       straightness: endpoint to endpoint Euclidean distance divided by
%       the lenght of the link
%       dt_e2e_2_length: ratio between maximum DT of the link voxel and the
%       length of the link 
%       dt_mmxx_2_length: ratio between the maximum and minimum DT of the
%       link voxels and the link length
%       dt_std_n: standard deviation of the link DT, normalized by the mean
%       link DT
%       dt_e2e_2_ep_dist: ratio between the link endpoints DT difference
%       and the link endpoint-endpoint distance
%       int_std_n: standard deviation of the link intensity divided by the
%       mean link intensity
%       nearest_ep_dist: distance to the closest endpoint 
%       nearest_ep_idx: index of the closest endpoint in the ep_sub list
%       ep_ep_disp_vec_n: normalized displacement vector between two
%       endpoints, pointing from the one to its nearest neighbor
%       inner_product_epv_epv: inner product of two direction vectors of
%       the endpoint
%       inner_product_epv1_dispv: inner product of the first endpoint
%       direction vector and the endpoint displacement vector
%       inner_product_epv2_dispv: inner product of the second endpoint
%       direction vector and the endpoint displacement vector (flip
%       direction)

% Parameters
local_bbox_expansion = 10;
pca_max_num_voxel = 10;
image_size = size(vessel_image);
% Initialization
lf = struct;
num_l = numel(link_cc);
lf.num_voxel = cellfun(@length, link_cc);
lf.length = nan(num_l, 1);
lf.ep2ep_dist = ones(num_l, 1);
lf.dt_max = nan(num_l, 1);
lf.dt_min = nan(num_l, 1);
lf.dt_mean = nan(num_l, 1);
lf.dt_std = nan(num_l, 1);
lf.dt_diff_ep2ep = nan(num_l, 1);
lf.int_mean = nan(num_l, 1);
lf.int_median = nan(num_l, 1);
lf.int_std = nan(num_l, 1);
lf.int_diff_ep2ep = nan(num_l, 1);
lf.int_min = nan(num_l, 1);
% lf.node_int = nan(num_l, 1);
lf.ep_direction_vec = nan(3, num_l);
ep_sub = nan(3, num_l);
% Information about the background statistics
lf.bg_mean = nan(num_l, 1);
lf.bg_std = nan(num_l, 1);
lf.mask_mean = nan(num_l, 1);
lf.mask_std = nan(num_l, 1);
vessel_mask = vessel_mask_dt > 0;
%% Compute single link features
for iter_link = 1 : num_l
    tmp_ind = link_cc{iter_link};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    tmp_ep_sub = tmp_sub(end,:);
    ep_sub(:, iter_link) = tmp_ep_sub';
    lf.length(iter_link) = fun_graph_sub_to_length(tmp_sub, 1);
    if size(tmp_sub,1) > 1
        lf.ep2ep_dist(iter_link) = sqrt(sum((tmp_sub(1, :) - tmp_sub(end,:)).^2));
    end
    if size(tmp_sub,1) > 2
        tmp_1_sub = tmp_sub(max(1, (end-min(pca_max_num_voxel, size(tmp_sub,1)) + 1)):end,:);
        [ep_1_vec] = pca(tmp_1_sub, 'Algorithm', 'svd', 'Centered',true, 'NumComponents', 1);
        % make sure the vector is point along the direction of the vessel: 
        ep_voxel_vec = tmp_1_sub(end, :) - tmp_1_sub(1, :);
        tmp_inner_produce = ep_voxel_vec * ep_1_vec;
        if tmp_inner_produce < 0
            lf.ep_direction_vec(:, iter_link) = - ep_1_vec;
        else
            lf.ep_direction_vec(:, iter_link) = ep_1_vec;
        end
    end 
    tmp_dt = vessel_mask_dt(tmp_ind);
    tmp_int = single(vessel_image(tmp_ind));
%     tmp_corr = corrcoef(tmp_dt, tmp_int);
%     if numel(tmp_dt) > 1
%         lf.dt_int_correlation(iter_link) = tmp_corr(1,2);
%     end
    tmp_dt = tmp_dt(tmp_dt>0);
    if ~isempty(tmp_dt)
        lf.dt_max(iter_link) = max(tmp_dt);
        lf.dt_min(iter_link) = min(tmp_dt);
        lf.dt_mean(iter_link) = mean(tmp_dt);
        lf.dt_std(iter_link) = std(tmp_dt);
        lf.dt_diff_ep2ep(iter_link) = abs(tmp_dt(end) - tmp_dt(1));
    end
    lf.int_mean(iter_link) = mean(tmp_int);
    lf.int_median(iter_link) = median(tmp_int);
    lf.int_std(iter_link) = std(tmp_int);
    lf.int_diff_ep2ep(iter_link) = abs(tmp_int(end) - tmp_int(1));
    lf.int_min(iter_link) = min(tmp_int);
    
    % Background information
    tmp_bg_str = fun_graph_get_link_feature_bg(tmp_sub, vessel_image, vessel_mask, local_bbox_expansion);
    lf.bg_mean(iter_link) = tmp_bg_str.bg_mean;
    lf.bg_std(iter_link) = tmp_bg_str.bg_std;
    lf.mask_mean(iter_link) = tmp_bg_str.mask_mean;
    lf.mask_std(iter_link) = tmp_bg_str.mask_std;
end
lf.straightness = lf.ep2ep_dist ./ lf.length;
lf.dt_e2e_2_length = lf.dt_diff_ep2ep ./ lf.length;
lf.dt_max_2_length = lf.dt_max ./ lf.length;
lf.dt_mmxx_2_length = (lf.dt_max - lf.dt_min) ./ lf.length;
lf.dt_std_n = lf.dt_std ./ lf.dt_mean;
lf.dt_e2e_2_ep_dist = lf.dt_diff_ep2ep ./ lf.ep2ep_dist;
lf.int_std_n = lf.int_std ./ lf.int_mean;

lf.int_above_bg = lf.int_mean - lf.bg_mean;
lf.int_SNR = lf.int_above_bg ./ lf.bg_std;
%% Compute endpoint features
lf.ep_direction_vec = lf.ep_direction_vec';
lf.ep_dir_vec_z = abs(lf.ep_direction_vec(:, 3));
ep_sub = ep_sub';
lf.ep_sub = ep_sub;
dist_ep2ep = squareform(pdist(ep_sub));
num_ep = size(ep_sub, 1);
if num_ep > 1
    dist_ep2ep = dist_ep2ep + diag(ones(num_ep,1) * inf);
    % Find the nearest endpoint for each endpoint, not necessary mutually
    % matched - otherwise, there will be missing features
    [lf.nearest_ep_dist, lf.nearest_ep_idx] = min(dist_ep2ep, [], 2);

    lf.nearest_link_num_voxel = lf.num_voxel(lf.nearest_ep_idx);
    % Inner produce between the two endpoint vector: 
    tmp_matched_ep_vec = lf.ep_direction_vec(lf.nearest_ep_idx, :);
    tmp_matched_ep_sub = ep_sub(lf.nearest_ep_idx, :);
    tmp_ep_ep_disp_vec = tmp_matched_ep_sub - ep_sub;
    tmp_ep_ep_disp_vec = tmp_ep_ep_disp_vec ./ vecnorm(tmp_ep_ep_disp_vec, 2, 2);
    lf.ep_ep_disp_vec_n = tmp_ep_ep_disp_vec;

    lf.inner_product_epv_epv = sum(lf.ep_direction_vec .* tmp_matched_ep_vec, 2);
    lf.inner_product_epv1_dispv = sum(tmp_ep_ep_disp_vec .* lf.ep_direction_vec, 2);
    lf.inner_product_epv2_dispv = sum(- tmp_ep_ep_disp_vec .* tmp_matched_ep_vec, 2);
else
    [lf.nearest_ep_dist, lf.nearest_ep_idx, lf.nearest_link_num_voxel, ...
        lf.inner_product_epv_epv, lf.inner_product_epv1_dispv, ...
        lf.inner_product_epv2_dispv] = deal(nan);
    lf.ep_ep_disp_vec_n = nan(1, 3);    
end
if num_l == 1
   lf = struct2table(lf, 'AsArray', true);  
else
   lf = struct2table(lf); 
end
end
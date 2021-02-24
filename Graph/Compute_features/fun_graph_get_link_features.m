function lf = fun_graph_get_link_features(link_cc, vessel_image, vessel_mask_dt, link_feature)
% fun_graph_get_link_w_1_ep_features computes the features derived from the
% position and orientation of endpoints and their nearest neighbor
% Input: 
%   ep_sub: N-by-3 numerical array, where N is the number of endpoint and 3
%   is the 3D subscript of the endpoint in the array
%   ep_vec: 3-by-N numerical array, the normalize direction vector of the
%   endpoint, which is computed as the first normalized principle component
%   of the link that the endpoint attached to, which is pointing along the link and
%   ends at the endpoint. 
%   link_feature: cell array consisit of the combination of {'geometry', 'dt',
%   int', bg'}
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
% Note: 
% 1. 2/3 of the time is spent on PCA, the other about 1/3 of the time is
% for cropping the arrays for estimating the background level. 
if nargin < 4
    link_feature = {'geometry', 'dt', 'int', 'bg'};
end
    

% Parameters
local_bbox_expansion = 10;
% pca_max_num_voxel = 10;
if ~isempty(vessel_image)
    image_size = size(vessel_image);
    if ~isempty(vessel_mask_dt)
        vessel_mask = vessel_mask_dt > 0;
    end
elseif ~isempty(vessel_mask_dt)
    image_size = size(vessel_mask_dt);
    vessel_mask = vessel_mask_dt > 0;
else
    error('Both vessel_image and vessel_mask_dt are empty. For only geometric properties, use fun_graph_get_link_geometric_features');
end
% Initialization
lf = struct;
num_l = numel(link_cc);
lf.num_voxel = cellfun(@length, link_cc);
lf.length = nan(num_l, 1);
lf.ep2ep_dist = nan(num_l, 1);
lf.dt_ep1 = nan(num_l, 1);
lf.dt_ep2 = nan(num_l, 1);
lf.dt_max = nan(num_l, 1);
lf.dt_min = nan(num_l, 1);
lf.dt_mean = nan(num_l, 1);
lf.dt_median = nan(num_l, 1);
lf.dt_std = nan(num_l, 1);
lf.dt_diff_ep2ep = nan(num_l, 1);
lf.int_mean = nan(num_l, 1);
lf.int_median = nan(num_l, 1);
lf.int_min = nan(num_l, 1);
lf.int_max = nan(num_l, 1);
lf.int_middle_point = nan(num_l, 1);
lf.int_min_idx_n = nan(num_l, 1);
lf.int_min_from_mid = nan(num_l, 1);
lf.int_max_at_end = nan(num_l, 1);
lf.int_ep1 = nan(num_l, 1);
lf.int_ep2 = nan(num_l, 1);

lf.int_std = nan(num_l, 1);
lf.int_diff_ep2ep = nan(num_l, 1);
% lf.node_int = nan(num_l, 1);
lf.link_com = nan(3, num_l);
% ep_sub = nan(3, num_l);
% Information about the background statistics
lf.bg_mean = nan(num_l, 1);
lf.bg_std = nan(num_l, 1);
lf.mask_mean = nan(num_l, 1);
lf.mask_std = nan(num_l, 1);

lf.ep1_direction_vec = nan(3, num_l);
lf.ep2_direction_vec = nan(3, num_l);
lf.ep1_to_ep2_direction_vec = nan(3, num_l);
% lf.cc_direction_vec = nan(3, num_l);
% lf.dt_int_correlation = nan(lf.num_l, 1);
pca_max_num_voxel = 10;
%% Compute single link features
for iter_link = 1 : num_l
    tmp_ind = link_cc{iter_link};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    tmp_num_voxel = lf.num_voxel(iter_link);
    tmp_middle_point = floor((tmp_num_voxel - 1)/2) + 1;
    if ismember('geometry', link_feature)
        lf.link_com(:, iter_link) = mean(tmp_sub, 1)';
        %     tmp_ep_sub = tmp_sub(end,:);
        %     ep_sub(:, iter_link) = tmp_ep_sub';
        lf.length(iter_link) = fun_graph_sub_to_length(tmp_sub, 1);        
        if tmp_num_voxel > 1
            tmp_ep1_to_ep2_voxel_vec = tmp_sub(end, :) - tmp_sub(1, :);
            tmp_ep1_to_ep2_vec_norm = sqrt(sum((tmp_ep1_to_ep2_voxel_vec ).^2));
            lf.ep2ep_dist(iter_link) = tmp_ep1_to_ep2_vec_norm;            
            tmp_ep1_to_ep2_voxel_vec = tmp_ep1_to_ep2_voxel_vec ./ tmp_ep1_to_ep2_vec_norm;
            lf.ep1_to_ep2_direction_vec(:, iter_link) = tmp_ep1_to_ep2_voxel_vec;
        end
        if size(tmp_sub,1) > 2
            % unit orientation vector of two ends and between the two end
            % make sure the vector is point along the direction of the vessel:
            tmp_1_sub = tmp_sub(1 : min(pca_max_num_voxel, size(tmp_sub,1)),:);
            [ep_1_vec] = pca(tmp_1_sub, 'Algorithm', 'svd', 'Centered',true, 'NumComponents', 1);
            tmp_inner_produce = tmp_ep1_to_ep2_voxel_vec * ep_1_vec;
            if tmp_inner_produce > 0
                lf.ep1_direction_vec(:, iter_link) = - ep_1_vec;
            else
                lf.ep1_direction_vec(:, iter_link) = ep_1_vec;
            end
            
            tmp_2_sub = tmp_sub(max(1, (end-min(pca_max_num_voxel, size(tmp_sub,1)) + 1)):end,:);
            [ep_2_vec] = pca(tmp_2_sub, 'Algorithm', 'svd', 'Centered',true, 'NumComponents', 1);
            tmp_inner_produce = tmp_ep1_to_ep2_voxel_vec * ep_2_vec;
            if tmp_inner_produce < 0
                lf.ep2_direction_vec(:, iter_link) = - ep_2_vec;
            else
                lf.ep2_direction_vec(:, iter_link) = ep_2_vec;
            end
        end
    end
    if ismember('dt', link_feature)
        tmp_dt = vessel_mask_dt(tmp_ind);
        lf.dt_ep1(iter_link) = tmp_dt(1);
        lf.dt_ep2(iter_link) = tmp_dt(end);
        lf.dt_diff_ep2ep(iter_link) = abs(tmp_dt(end) - tmp_dt(1));
        % Some link might be out of the original mask 
        tmp_dt = tmp_dt(tmp_dt > 0);
        if ~isempty(tmp_dt)
            lf.dt_max(iter_link) = max(tmp_dt);
            lf.dt_min(iter_link) = min(tmp_dt);
            lf.dt_mean(iter_link) = mean(tmp_dt);
            lf.dt_std(iter_link) = std(tmp_dt);
            lf.dt_median(iter_link) = median(tmp_dt);
        end
    end
    
    if ismember('int', link_feature)
        tmp_int = single(vessel_image(tmp_ind));
        lf.int_mean(iter_link) = mean(tmp_int);
        lf.int_median(iter_link) = median(tmp_int);
        lf.int_std(iter_link) = std(tmp_int);
        lf.int_ep1(iter_link) = tmp_int(1);
        lf.int_ep2(iter_link) = tmp_int(end);
        lf.int_diff_ep2ep(iter_link) = abs(tmp_int(end) - tmp_int(1));
        lf.int_max(iter_link) = max(tmp_int);
        [lf.int_min(iter_link), tmp_idx]= min(tmp_int);
        lf.int_min_idx_n(iter_link) = tmp_idx / tmp_num_voxel;
        lf.int_min_from_mid(iter_link) = abs(0.5 - lf.int_min_idx_n(iter_link))/0.5;
        lf.int_middle_point(iter_link) = tmp_int(tmp_middle_point);
        lf.int_max_at_end(iter_link) = max(tmp_int(1), tmp_int(end));
    end
    % Background information
    if ismember('bg', link_feature)
        tmp_bg_str = fun_graph_get_link_feature_bg(tmp_sub, vessel_image, vessel_mask, local_bbox_expansion);
        lf.bg_mean(iter_link) = tmp_bg_str.bg_mean;
        lf.bg_std(iter_link) = tmp_bg_str.bg_std;
        lf.mask_mean(iter_link) = tmp_bg_str.mask_mean;
        lf.mask_std(iter_link) = tmp_bg_str.mask_std;
    end
end
lf.straightness = lf.ep2ep_dist ./ lf.length;
lf.dt_e2e_2_length = lf.dt_diff_ep2ep ./ lf.length;
lf.dt_max_2_length = lf.dt_max ./ lf.length;
lf.dt_mmxx_2_length = (lf.dt_max - lf.dt_min) ./ lf.length;
lf.dt_std_n = lf.dt_std ./ lf.dt_mean;
lf.dt_e2e_2_ep_dist = lf.dt_diff_ep2ep ./ lf.ep2ep_dist;
lf.dt_ep1_plus_ep2 = lf.dt_ep1 + lf.dt_ep2;
lf.dt_ep_sum_2_length = lf.dt_ep1_plus_ep2 ./ lf.length;
lf.dt_ep_sum_2_ep_dist = lf.dt_ep1_plus_ep2 ./ lf.ep2ep_dist;

lf.int_std_n = lf.int_std ./ lf.int_mean;
lf.bg_std_n = lf.bg_std ./ lf.bg_mean;
lf.mask_std_n = lf.mask_std ./ lf.mask_mean;
lf.int_mid_2_max = lf.int_middle_point ./ lf.int_max;
lf.int_min_2_max = lf.int_min ./ lf.int_max;
lf.int_mid_2_max_at_end = lf.int_middle_point ./ lf.int_max_at_end;
lf.int_ep_avg = (lf.int_ep1 + lf.int_ep2) ./ 2;
lf.int_diff_epavg_mid = lf.int_ep_avg - lf.int_middle_point;
lf.int_diff_epavg_mid_n = lf.int_diff_epavg_mid ./ lf.int_ep_avg;

lf.mid_point_SNR = (lf.int_middle_point - lf.bg_mean) ./ lf.bg_std;

lf.mask_SNR = (lf.mask_mean - lf.bg_mean) ./ lf.bg_std;
lf.skl_SNR = (lf.int_mean - lf.bg_mean) ./ lf.bg_std;
lf.link_com = lf.link_com';

lf.ep1_direction_vec = lf.ep1_direction_vec';
lf.ep2_direction_vec = lf.ep2_direction_vec';
lf.ep1_to_ep2_direction_vec = lf.ep1_to_ep2_direction_vec';

if numel(lf.num_voxel) == 1
    lf = struct2table(lf, 'AsArray', true);
else
    lf = struct2table(lf);
end
end
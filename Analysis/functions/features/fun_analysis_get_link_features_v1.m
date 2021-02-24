function lf = fun_analysis_get_link_features_v1(link_cc, vessel_image, vessel_mask_dt, link_feature)
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
if nargin < 4
    link_feature = {'geometry', 'dt', 'int', 'bg'};
end
if ismember('geometry', link_feature)
    compute_geometryQ = true;
else
    compute_geometryQ = false;
end
if ismember('dt', link_feature)
    compute_dtQ = true;
else
    compute_dtQ = false;
end
if ismember('int', link_feature)
    compute_intQ = true;
else
    compute_intQ = false;
end
if ismember('bg', link_feature)
    compute_bgQ = true;
else
    compute_bgQ = false;
end
% Parameters
local_bbox_expansion = 10;
% pca_max_num_voxel = 10;
if ~isempty(vessel_image)
    image_size = size(vessel_image);
    vessel_image = single(vessel_image);
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
if compute_geometryQ
    pca_max_num_voxel = 10;
    lf.length = zeros(num_l, 1);
    lf.link_com = zeros(3, num_l);
    lf.ep2ep_dist = ones(num_l, 1);
    lf.ep1_to_ep2_direction_vec = zeros(3, num_l);
    lf.ep1_direction_vec = zeros(3, num_l);
    lf.ep2_direction_vec = zeros(3, num_l);
    lf.cc_sub_pca1_vec = zeros(3, num_l);
    lf.cc_sub_pca2_vec = zeros(3, num_l);
    lf.cc_sub_pca3_vec = zeros(3, num_l);
    lf.cc_sub_cov_eig_val = zeros(3, num_l);
end
if compute_dtQ
    lf.dt_ep1 = ones(num_l, 1);
    lf.dt_ep2 = ones(num_l, 1);
    lf.dt_max = zeros(num_l, 1);
    lf.dt_min = zeros(num_l, 1);
    lf.dt_mean = zeros(num_l, 1);
    lf.dt_median = zeros(num_l, 1);
    lf.dt_std = zeros(num_l, 1);
    lf.dt_diff_ep2ep = zeros(num_l, 1);
end
if compute_intQ
    lf.int_mean = zeros(num_l, 1);
    lf.int_median = zeros(num_l, 1);
    lf.int_min = zeros(num_l, 1);
    lf.int_max = zeros(num_l, 1);
    lf.int_middle_point = zeros(num_l, 1);
    lf.int_min_idx_n = zeros(num_l, 1);
    lf.int_max_at_end = zeros(num_l, 1);
    lf.int_std = zeros(num_l, 1);
    lf.int_diff_ep2ep = zeros(num_l, 1);
end
% Information about the background statistics
if compute_bgQ
    lf.bg_mean = zeros(num_l, 1);
    lf.bg_std = zeros(num_l, 1);
    lf.mask_mean = zeros(num_l, 1);
    lf.mask_std = zeros(num_l, 1);
end

%% Compute single link features
for iter_link = 1 : num_l
    tmp_ind = link_cc{iter_link};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    if compute_geometryQ
        tmp_sub_mean = mean(tmp_sub, 1)';
        lf.link_com(:, iter_link) = tmp_sub_mean;
        %     tmp_ep_sub = tmp_sub(end,:);
        %     ep_sub(:, iter_link) = tmp_ep_sub';
        lf.length(iter_link) = fun_graph_ind_to_length(tmp_sub, 1);
        tmp_num_voxel = numel(tmp_ind);
        
        if tmp_num_voxel > 1
            tmp_ep1_to_ep2_voxel_vec = tmp_sub(end, :) - tmp_sub(1, :);
            tmp_ep1_to_ep2_vec_norm = sqrt(sum((tmp_ep1_to_ep2_voxel_vec ).^2));
            lf.ep2ep_dist(iter_link) = tmp_ep1_to_ep2_vec_norm;            
            tmp_ep1_to_ep2_voxel_vec = tmp_ep1_to_ep2_voxel_vec ./ tmp_ep1_to_ep2_vec_norm;
            lf.ep1_to_ep2_direction_vec(:, iter_link) = tmp_ep1_to_ep2_voxel_vec;
        end
        tmp_middle_point = floor((tmp_num_voxel - 1)/2) + 1;
        if tmp_num_voxel > 2
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
            
            % PCA by SVD
            tmp_sub_centered = tmp_sub - tmp_sub_mean';
            cov_mat = (tmp_sub_centered' * tmp_sub_centered) ./ (tmp_num_voxel - 1);
            [svd_vec, svd_val, ~] = svd(cov_mat);
            lf.cc_sub_pca1_vec(:, iter_link) = svd_vec(:,1);
            lf.cc_sub_pca2_vec(:, iter_link) = svd_vec(:,2);
            lf.cc_sub_pca3_vec(:, iter_link) = svd_vec(:,3);
            lf.cc_sub_cov_eig_val(:, iter_link) = svd_val([1, 5, 9]);
%             [cc_pca_vecs, ~, cc_pac_var] = pca(tmp_sub, 'Algorithm', 'svd', 'Centered', true, 'NumComponents', 3);
%             [cc_vec] = pac(tmp_sub, 'Algorithm', 'svd', 'Centered', true, 'NumCompoents', 1);
        end
    end
    if compute_dtQ
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
    
    if compute_intQ
        tmp_int = single(vessel_image(tmp_ind));  
        lf.int_mean(iter_link) = mean(tmp_int);
        lf.int_median(iter_link) = median(tmp_int);
        lf.int_std(iter_link) = std(tmp_int);
        lf.int_diff_ep2ep(iter_link) = abs(tmp_int(end) - tmp_int(1));
        lf.int_max(iter_link) = max(tmp_int);
        [lf.int_min(iter_link), tmp_idx]= min(tmp_int);
        lf.int_min_idx_n(iter_link) = tmp_idx / tmp_num_voxel;
        lf.int_middle_point(iter_link) = tmp_int(tmp_middle_point);
        lf.int_max_at_end(iter_link) = max(tmp_int(1), tmp_int(end));
    end
    % Background information
    if compute_bgQ
        tmp_min = max([1,1,1],  tmp_sub - local_bbox_expansion);
        tmp_max = min(image_size, tmp_sub + local_bbox_expansion);
        tmp_mask = vessel_mask(tmp_min(1):tmp_max(1), tmp_min(2):tmp_max(2), tmp_min(3):tmp_max(3));
        tmp_mask = tmp_mask(:);
        tmp_bg_mask = ~tmp_mask;
        tmp_image = vessel_image(tmp_min(1):tmp_max(1), tmp_min(2):tmp_max(2), tmp_min(3):tmp_max(3));
        tmp_image = tmp_image(:);
        tmp_fg_voxel = nnz(tmp_mask);
        tmp_bg_voxel = nnz(tmp_bg_mask);

        lf.bg_mean(iter_link) = sum(tmp_image(:) .* tmp_bg_mask(:)) /tmp_bg_voxel;
        lf.bg_std(iter_link) = sqrt(sum(tmp_image(:).^2 .* tmp_bg_mask(:)) /tmp_bg_voxel - lf.bg_mean(iter_link)^2);
        lf.mask_mean(iter_link) = sum(tmp_image .* tmp_mask) /tmp_fg_voxel;
        lf.mask_std(iter_link) = sqrt(sum(tmp_image.^2 .* tmp_mask) /tmp_fg_voxel - lf.bg_mean(iter_link)^2);
    end
end

if compute_geometryQ
    lf.straightness = lf.ep2ep_dist ./ lf.length;
    lf.link_com = lf.link_com';
    lf.ep1_direction_vec = lf.ep1_direction_vec';
    lf.ep2_direction_vec = lf.ep2_direction_vec';
    lf.ep1_to_ep2_direction_vec = lf.ep1_to_ep2_direction_vec';
    lf.cc_sub_pca1_vec = lf.cc_sub_pca1_vec';
    lf.cc_sub_pca2_vec = lf.cc_sub_pca2_vec';
    lf.cc_sub_pca3_vec = lf.cc_sub_pca3_vec';
    lf.cc_sub_cov_eig_val = lf.cc_sub_cov_eig_val';
    % Convert the orientation end to end vector to spherical coordinate
    % Fix the z-component of the vector to be nonnegative
    tmp_vec = lf.ep1_to_ep2_direction_vec;
    negative_z_Q = tmp_vec(:,3) < 0;
    tmp_vec(negative_z_Q, :) = - tmp_vec(negative_z_Q, :);
    [lf.ep2ep_angle_azimuth_deg, lf.ep2ep_angle_elevation_deg, ~] = cart2sph(tmp_vec(:,2), tmp_vec(:,1), tmp_vec(:,3));
    % Convert to degree
    lf.ep2ep_angle_azimuth_deg = lf.ep2ep_angle_azimuth_deg .* 180 ./ pi;
    lf.ep2ep_angle_elevation_deg = lf.ep2ep_angle_elevation_deg .* 180 ./ pi;
end
if compute_dtQ 
    lf.dt_std_n = lf.dt_std ./ lf.dt_mean;
    lf.dt_e2e_2_ep_dist = lf.dt_diff_ep2ep ./ lf.ep2ep_dist;
    lf.dt_ep1_plus_ep2 = lf.dt_ep1 + lf.dt_ep2;
    if compute_geometryQ
        lf.dt_e2e_2_length = lf.dt_diff_ep2ep ./ lf.length;
        lf.dt_max_2_length = lf.dt_max ./ lf.length;
        lf.dt_mmxx_2_length = (lf.dt_max - lf.dt_min) ./ lf.length;
        lf.dt_ep_sum_2_length = lf.dt_ep1_plus_ep2 ./ lf.length;
    end
end
if compute_intQ
    lf.int_std_n = lf.int_std ./ lf.int_mean;
    lf.int_mid_2_max = lf.int_middle_point ./ lf.int_max;
    lf.int_min_2_max = lf.int_min ./ lf.int_max;
    lf.int_mid_2_max_at_end = lf.int_middle_point ./ lf.int_max_at_end;
end
if compute_bgQ
    lf.mask_std_n = lf.mask_std ./ lf.mask_mean;
    lf.bg_std_n = lf.bg_std ./ lf.bg_mean;
    lf.mid_point_SNR = (lf.int_middle_point - lf.bg_mean) ./ lf.bg_std;
    lf.mask_SNR = (lf.mask_mean - lf.bg_mean) ./ lf.bg_std;
    lf.skl_SNR = (lf.int_mean - lf.bg_mean) ./ lf.bg_std;
end
if num_l == 1
    lf = struct2table(lf, 'AsArray', true);
else
    lf = struct2table(lf);
end
end
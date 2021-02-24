function lf = fun_analysis_get_link_features_by_label(vessel_graph, link_label, link_feature)
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
if nargin < 2
    link_label = [];
    link_feature = {'geometry', 'dt'};
elseif nargin < 3
    link_feature = {'geometry', 'dt'};
end
if isempty(link_label)
    link_label = 1 : vessel_graph.link.num_cc;
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
image_size = vessel_graph.num.mask_size;
pca_max_num_voxel = 10;
% Initialization
link_cc = vessel_graph.link.cc_ind(link_label);
vessel_mask_dt = vessel_graph.radius;
lf = struct;
num_l = numel(link_cc);
if compute_geometryQ
    pca_max_num_voxel = 10;
    if isfield(vessel_graph, 'info')
        compute_global_ind_Q = true;
        [lf.mid_global_ind, lf.ep1_global_ind, lf.ep2_global_ind] = deal(nan(num_l, 1));
         lf.mid_global_sub = nan(num_l, 3);
    else
        compute_global_ind_Q = false;
    end
    [lf.length, lf.ep2ep_dist] = deal(nan(num_l, 1));
    
    [lf.link_com, lf.ep1_sub, lf.ep2_sub, lf.ep1_to_ep2_direction_vec, ...
        lf.ep1_direction_vec, lf.ep2_direction_vec, lf.cc_sub_pca1_vec, ...
        lf.cc_sub_pca2_vec, lf.cc_sub_pca3_vec, lf.cc_sub_cov_eig_val, ...
        lf.mid_sub] = deal(nan(3, num_l));
end
if compute_dtQ
    lf.dt_ep1 = nan(num_l, 1);
    lf.dt_ep2 = nan(num_l, 1);
    lf.dt_max = nan(num_l, 1);
    lf.dt_min = nan(num_l, 1);
    lf.dt_mean = nan(num_l, 1);
    lf.dt_median = nan(num_l, 1);
    lf.dt_std = nan(num_l, 1);
end
%% Compute single link features
for iter_link = 1 : num_l
    tmp_ind = link_cc{iter_link};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    if compute_geometryQ
        tmp_sub_mean = mean(tmp_sub, 1)';
        lf.link_com(:, iter_link) = tmp_sub_mean;
        lf.length(iter_link) = fun_graph_sub_to_length(tmp_sub, 1);
        tmp_num_voxel = numel(tmp_ind);        
        lf.ep1_sub(:, iter_link) = tmp_sub(1, :);
        lf.ep2_sub(:, iter_link) = tmp_sub(end, :);
        if tmp_num_voxel > 1            
            tmp_ep1_to_ep2_voxel_vec = tmp_sub(end, :) - tmp_sub(1, :);
            tmp_ep1_to_ep2_vec_norm = sqrt(sum((tmp_ep1_to_ep2_voxel_vec ).^2));
            lf.ep2ep_dist(iter_link) = tmp_ep1_to_ep2_vec_norm + 1;   % Add 1 to make it the same for straight line.          
            tmp_ep1_to_ep2_voxel_vec = tmp_ep1_to_ep2_voxel_vec ./ tmp_ep1_to_ep2_vec_norm;
            lf.ep1_to_ep2_direction_vec(:, iter_link) = tmp_ep1_to_ep2_voxel_vec;
        end
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
            
            % Record middle point local and global coordinate
            tmp_mid_idx = ceil(tmp_num_voxel / 2);
            lf.mid_sub(:, iter_link) = tmp_sub(tmp_mid_idx, :);  
        end
    end
    if compute_dtQ
        tmp_dt = full(vessel_mask_dt(tmp_ind));
        lf.dt_ep1(iter_link) = tmp_dt(1);
        lf.dt_ep2(iter_link) = tmp_dt(end);
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
end

if compute_geometryQ
    lf.straightness = lf.ep2ep_dist ./ lf.length;
    lf.link_com = lf.link_com';
    
    lf.ep1_sub = lf.ep1_sub';

    lf.ep2_sub = lf.ep2_sub';
    
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
    % Compute global location of the endpoints and the middle point of the link
    lf.mid_sub = lf.mid_sub';
    if compute_global_ind_Q
        tmp_sub_global = lf.ep1_sub + vessel_graph.info.bbox_mmxx(1:3) - 1;
        lf.ep1_global_ind = sub2ind(vessel_graph.info.dataset_size, tmp_sub_global(:, 1), ...
            tmp_sub_global(:, 2), tmp_sub_global(:, 3));
        tmp_sub_global = lf.ep2_sub + vessel_graph.info.bbox_mmxx(1:3) - 1;
        lf.ep2_global_ind = sub2ind(vessel_graph.info.dataset_size, tmp_sub_global(:, 1), ...
            tmp_sub_global(:, 2), tmp_sub_global(:, 3)); 
        lf.mid_global_sub = lf.mid_sub + vessel_graph.info.bbox_mmxx(1:3) - 1;
        lf.mid_global_ind = sub2ind(vessel_graph.info.dataset_size, ...
            lf.mid_global_sub(:, 1), lf.mid_global_sub(:, 2), lf.mid_global_sub(:, 3));
    end
end
if compute_dtQ 
    if compute_geometryQ
        % Surface area
        lf.surface_area = 2 * pi * lf.length .* lf.dt_median;
        % Volume
        lf.volume = pi * lf.length .* (lf.dt_median .^ 2);
    end
end
if num_l == 1
    lf = struct2table(lf, 'AsArray', true);
else
    lf = struct2table(lf);
end
end
function lf = fun_graph_get_link_geometric_features(link_cc, image_size)
% fun_graph_get_link_w_1_ep_features computes the features derived from the
% position and orientation of endpoints and their nearest neighbor
% Input: 
%   link_cc: cell array, each cell contains a numerical vector, which is
%   the indices of the voxels in a link 
%   image_size: size of the image
% Output: 
%   lf: table with fields
%
% 

% Parameters
pca_max_num_voxel = 10;
% Initialization
lf = struct;
num_l = numel(link_cc);
lf.num_voxel = cellfun(@length, link_cc);
lf.length = zeros(num_l, 1);
lf.ep2ep_dist = ones(num_l, 1);
lf.link_com = zeros(3, num_l);
% ep_sub = zeros(3, num_l);
lf.ep1_direction_vec = zeros(3, num_l);
lf.ep2_direction_vec = zeros(3, num_l);
lf.ep1_to_ep2_direction_vec = zeros(3, num_l);
% lf.cc_direction_vec = zeros(3, num_l);
% lf.dt_int_correlation = nan(lf.num_l, 1);
%% Compute single link features
for iter_link = 1 : num_l
    tmp_ind = link_cc{iter_link};
    tmp_sub = fun_ind2sub(image_size, tmp_ind);
    lf.link_com(:, iter_link) = mean(tmp_sub, 1)';
%     tmp_ep_sub = tmp_sub(end,:);
%     ep_sub(:, iter_link) = tmp_ep_sub';
    lf.length(iter_link) = fun_graph_ind_to_length(tmp_sub, 1);
    tmp_num_voxel = lf.num_voxel(iter_link);
    
    if tmp_num_voxel > 1
        lf.ep2ep_dist(iter_link) = sqrt(sum((tmp_sub(1, :) - tmp_sub(end,:)).^2));
    end
    if size(tmp_sub,1) > 2
        % unit orientation vector of two ends and between the two end
        % make sure the vector is point along the direction of the vessel: 
        tmp_ep1_to_ep2_voxel_vec = tmp_sub(end, :) - tmp_sub(1, :);
        tmp_ep1_to_ep2_voxel_vec = tmp_ep1_to_ep2_voxel_vec ./ vecnorm(tmp_ep1_to_ep2_voxel_vec);
        lf.ep1_to_ep2_direction_vec(:, iter_link) = tmp_ep1_to_ep2_voxel_vec;
        
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
lf.straightness = lf.ep2ep_dist ./ lf.length;
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
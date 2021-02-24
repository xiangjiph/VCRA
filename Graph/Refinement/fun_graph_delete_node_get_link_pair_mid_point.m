function mid_point_str = fun_graph_delete_node_get_link_pair_mid_point(cc_sub_1, cc_sub_2, node_sub, node_dt, num_voxel_fit)
% fun_graph_delete_node_get_link_pair_mid_point fit a 3D line to the
% skeleton link segments subscript and find the closest point on this line
% to the node. 
% Input: 
%   cc_ind_1: N-by-3 numerical array, index of the first link segment. The first
%   element is in the neighbor of the node voxel 
%   cc_ind_2: same, for the second link segment
%   node_sub: N-by-3 numerical array, subscript of the node. 
%   node_dt: distance transform of the node, used for estimating the number
%   of link voxels need to be excluded from the line fitting. 
%   num_voxel_fit: the maximum number of voxels from each link segments are used for
%   line fitting. 
%
% Note: 
% 1. Assume the voxel that connected to the node is always the first element
% of the input cc ind vector
% 2. The minimum number of voxels used for fitting is 1(the last one) 
%
% Author: Xiang Ji
% Date: 02/15/2019
if nargin < 4
    node_dt = 3;
    num_voxel_fit = 8;
elseif nargin < 5
    num_voxel_fit = 8;    
end

if ~isvector(node_sub)
    disp('Input node subscript is not a vector. Compute mean subscript along the first axis');
    node_sub = mean(node_sub, 1);
end
% Remove the link voxel near the node. 
num_voxel_removed = round(max(node_dt/sqrt(2)));
num_cc_1 = size(cc_sub_1, 1);

kept_cc_sub_1 = cc_sub_1(min(num_voxel_removed + 1, num_cc_1):num_cc_1, :);
fit_sub_1 = kept_cc_sub_1(1 : min(num_voxel_fit, size(kept_cc_sub_1,1)), :);


num_cc_2 = size(cc_sub_2, 1);
kept_cc_sub_2 = cc_sub_2(min(num_voxel_removed + 1, num_cc_2):num_cc_2, :);
fit_sub_2 = kept_cc_sub_2(1 : min(num_voxel_fit, size(kept_cc_sub_2,1)), :);
% 3D linear regression 
fit_sub = [fit_sub_1; fit_sub_2];
fit_sub_mean = mean(fit_sub, 1);
fit_sub_c = fit_sub - fit_sub_mean;
cov_mat = (fit_sub_c' * fit_sub_c)./(size(fit_sub_c, 1) - 1);
% Diagonal the covariant matrix and get the largest eigenvalue and its
% eigenvector
[eig_vec, eig_val] = eig(cov_mat, 'vector');
[eig_val, tmp_idx] = sort(eig_val, 'descend');
eig_vec = eig_vec(:, tmp_idx);

% Closest point on the fitted line to the node position 
fit_vec = eig_vec(:,1);
sub_closest_to_node = fit_sub_mean' + ((node_sub - fit_sub_mean) * fit_vec) * fit_vec;

% Save information 
mid_point_str = struct;
mid_point_str.pos = sub_closest_to_node';
mid_point_str.sub = round(sub_closest_to_node');
mid_point_str.midpoint2node = sqrt(sum((node_sub - mid_point_str.pos).^2));
mid_point_str.link_sub_1 = kept_cc_sub_1;
mid_point_str.link_sub_2 = kept_cc_sub_2;

mid_point_str.sub_mean = fit_sub_mean;
mid_point_str.cov = cov_mat;
mid_point_str.cov_eig_vec = eig_vec;
mid_point_str.cov_eig_val = eig_val;
mid_point_str.fit_cov = fit_vec;
% Estimate the standard deviation in the radial direction 
mid_point_str.std_r = sqrt(eig_val(2) + eig_val(3));
end
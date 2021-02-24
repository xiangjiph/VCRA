function anisotropy_str = fun_analysis_get_link_anisotropy(ori_vec_list, ~, ~)
% fun_analysis_get_link_anisotropy computes the singular value
% decomposition on the link orientation vectors. The largest singluar value
% is divided by the smallest singluar value to quantify the anisotropy of
% the link orientation. 
% Input: 
%   ori_vec_list: N-by-3 vector. Can be the endpoint-to-endpoint unit
%   vector
%   undirected_Q: logical scalar.
% Output: 
%   anisotropy_str: structure with fields defined below. 
%
% Implemented by Xiang Ji on 03/22/2019
%
% Note: 
% 1. It seems that svd_1 itself has much larger variance compare to
% fractional anisotropy and min-max eigenvalue ratio for the same number of
% simulations. 
persistent str_template
    
assert(~iscolumn(ori_vec_list), 'The input orientation vector should be a N-by-3 vector');

if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'num_data', 'cov_mat', ...
        'svd_u', 'svd_value', 'svd_min2max', 'svd_value_sum', 'svd_value_norm', ...
        'svd_value_ratio', 'svd_max_vec'});
end
anisotropy_str = str_template;
is_valid_Q = all(~isnan(ori_vec_list), 2);
ori_vec_list = ori_vec_list(is_valid_Q, :);
anisotropy_str.num_data = size(ori_vec_list, 1);

if isempty(ori_vec_list) || size(ori_vec_list, 1) < 3
    return;
end

% ori_vec_list = cat(1, ori_vec_list, -ori_vec_list);
ori_vec_list = ori_vec_list + ori_vec_list .* (-2 * (ori_vec_list(:, 3) < 0));
cov_mat = (ori_vec_list' * ori_vec_list) ./ (size(ori_vec_list, 1) - 1);

[svd_U, svd_val, ~] = svd(cov_mat, 'econ');

anisotropy_str.cov_mat = cov_mat;
anisotropy_str.svd_u = svd_U;
anisotropy_str.svd_value = svd_val([1,5,9]);
anisotropy_str.svd_min2max = anisotropy_str.svd_value(3) / anisotropy_str.svd_value(1);
anisotropy_str.svd_value_sum = sum(anisotropy_str.svd_value);
anisotropy_str.svd_value_norm = sum(anisotropy_str.svd_value .^ 2);
anisotropy_str.svd_value_ratio = anisotropy_str.svd_value ./anisotropy_str.svd_value_sum;
anisotropy_str.svd_max_vec = svd_U(:, 1);
anisotropy_str.svd_value_mean = anisotropy_str.svd_value_sum / 3;
anisotropy_str.fractional_anisotropy = sqrt(3 * ...
    sum((anisotropy_str.svd_value - anisotropy_str.svd_value_mean) .^2) / ...
    (2 * anisotropy_str.svd_value_norm));
end
%% Debug
% debug_vec_list = ori_vec_list_0 + ori_vec_list_0 .* (-2 * double(ori_vec_list_0(:, 3) < 0));
% debug_cov_mat = (debug_vec_list' * debug_vec_list) ./ (size(debug_vec_list, 1) - 1);
% [debug_u, debug_v, ~] = svd(debug_cov_mat, 'econ');
% pca(debug_vec_list)
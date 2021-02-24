function anisotropy_str = fun_analysis_get_link_anisotropy_masking(ori_vec_list, masking)
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
if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'num_data', ...
        'svd_u', 'svd_value', 'svd_min2max', 'svd_value_sum', ...
        'svd_value_ratio', 'svd_max_vec', 'fa'});
end
assert(size(ori_vec_list, 2) == 3);

ori_vec_is_finite_Q = all(isfinite(ori_vec_list), 2);
if ~all(ori_vec_is_finite_Q)
    ori_vec_list = ori_vec_list(ori_vec_is_finite_Q, :);
    masking = masking(ori_vec_is_finite_Q, :);
end

anisotropy_str = str_template;
num_pts = size(ori_vec_list, 1);
if nargin < 2
    masking = true(num_pts, 1);
end
assert(size(masking, 1) == num_pts, 'Number of vectors in ori_vec_list is inconsistent with the number of rows in the masking matrix');

if isempty(ori_vec_list) || size(ori_vec_list, 1) < 3
    return;
end
num_selection = size(masking, 2);
svd_U_list = nan(3, 3, num_selection);
svd_vec_1_list = nan(3, num_selection);
svd_value_list = nan(3, num_selection);
num_selected_data = nan(1, num_selection);

for iter_selection = 1 : num_selection
    tmp_ori_vec = ori_vec_list(masking(:, iter_selection), :);
    tmp_num_pts = size(tmp_ori_vec, 1);
    if tmp_num_pts >= 3
        tmp_ori_vec = tmp_ori_vec + tmp_ori_vec .* (-2 * (tmp_ori_vec(:, 3) < 0));
        cov_mat = (tmp_ori_vec' * tmp_ori_vec) ./ (size(tmp_ori_vec, 1) - 1);
        [svd_u, svd_v, ~] = svd(cov_mat, 'econ');
        svd_U_list(:, :, iter_selection) = svd_u;
        svd_vec_1_list(:, iter_selection) = svd_u(:, 1);
        svd_value_list(:, iter_selection) = svd_v([1, 5, 9]);
        num_selected_data(iter_selection) = tmp_num_pts;
    end
end
anisotropy_str.num_data = num_selected_data;
anisotropy_str.svd_u = svd_U_list;
anisotropy_str.svd_value = svd_value_list;

anisotropy_str.svd_min2max = anisotropy_str.svd_value(3, :) ./ anisotropy_str.svd_value(1, :);
anisotropy_str.svd_value_sum = sum(svd_value_list, 1);
svd_value_norm = sum(anisotropy_str.svd_value .^ 2, 1);
anisotropy_str.svd_value_ratio = anisotropy_str.svd_value ./anisotropy_str.svd_value_sum;
anisotropy_str.svd_1 = anisotropy_str.svd_value_ratio(1, :);
anisotropy_str.svd_max_vec = svd_vec_1_list;
svd_value_mean = anisotropy_str.svd_value_sum ./ 3;
anisotropy_str.fa = sqrt(3 * ...
    sum((anisotropy_str.svd_value - svd_value_mean) .^2, 1) ./ ...
    (2 * svd_value_norm));
end








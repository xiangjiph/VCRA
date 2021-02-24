function iso_str = fun_simulation_get_weighted_isotropy_stat_with_masking(weight_vec, masking, num_simu)
% fun_simulation_get_weighted_isotropy_stat computed isotropy measure for
% orientation vector uniformly sampled from the unit sphere surface and
% weighted by the given weight. 
% This function is used for quantifying the deviation from isotropy for
% local capillary network and take the volume of the vessel into account,
% which might be more relavant for fMRI signal formation. 
% Input: 
%   weight_vec: numerical vector, whose length is the number of link for
%   simulation. 
%   num_simu: times of simulation for computing the mean and standard
%   deviation. 
% Output: 
%   iso_str: structure with fields. 
% 
persistent str_template;
if isempty(str_template)
    str_template = fun_initialized_structure_array_with_fieldname_list({'num_simu', ...
        'svd_min2max', 'fa', 'svd_1', 'data'});
end
iso_str = str_template;
num_pts = numel(weight_vec);
min_num_pts = 3;

if isrow(weight_vec)
    weight_vec = weight_vec .';
end

if nargin < 2
    masking = true(num_pts, 1);
    num_simu = 10000;
elseif nargin < 3
    num_simu = 10000;
end
iso_str.num_simu = num_simu;


if num_pts < min_num_pts
   % How to handle this case?  
   return;    
end
% Compute statistics in simulation 
% Determine the unique masking 
num_valid_data_point = sum(masking, 1);
is_valid_test_Q = num_valid_data_point >= min_num_pts;
masking(:, ~is_valid_test_Q) = 0;

[~, unique_masking_ind, broadcast_ind] = unique(masking.', 'rows', 'stable');
num_selection = numel(unique_masking_ind);

svd_list = nan(3, num_selection, num_simu);
is_valid_unique_masking = is_valid_test_Q(unique_masking_ind);
valid_unique_masking_ind = find(is_valid_unique_masking);
valid_selection_list_ind = unique_masking_ind(is_valid_unique_masking);
num_valid_selection = numel(valid_selection_list_ind);
for iter_trial = 1 : num_simu
    % Pair the weight with random orientation
    rand_ori = randn(3, num_pts);
    rand_ori = rand_ori ./ vecnorm(rand_ori);
    rand_ori = bsxfun(@times, rand_ori.', weight_vec);
    for iter_selection = 1 : num_valid_selection
        tmp_selection_ind = valid_selection_list_ind(iter_selection);
        selected_ori = rand_ori(masking(:, tmp_selection_ind), :);        
        selected_ori = selected_ori + selected_ori .* (-2 * (selected_ori(:, 3) < 0));
%         selected_ori = cat(1, selected_ori, -selected_ori);
        cov_mat = (selected_ori.' * selected_ori) ./ (size(selected_ori, 1) - 1);
        
        svd_val = svd(cov_mat, 'econ');
        
        svd_list(:, valid_unique_masking_ind(iter_selection), iter_trial) = svd_val;
    end
end
svd_list = svd_list(:, broadcast_ind, :);

svd_sum = sum(svd_list, 1);
svd_norm2 = sum(svd_list .^2, 1);

fractional_anisotropy = sqrt(3 * sum((svd_list - svd_sum ./ 3) .^2, 1) ./...
    (2 * svd_norm2));
fractional_anisotropy = permute(fractional_anisotropy, [3, 2, 1]);

svd_ratio_list = bsxfun(@rdivide, svd_list, svd_sum);
svd_ratio_list = permute(svd_ratio_list, [3, 2, 1]);
svd_min2max_list = svd_ratio_list(:, :, 3) ./ svd_ratio_list(:, :, 1);
%% Output
iso_str.is_valid_Q = is_valid_test_Q;

iso_str.svd_min2max = fun_analysis_get_basic_statistics_in_column(svd_min2max_list, 1);
iso_str.svd_1 = fun_analysis_get_basic_statistics_in_column(svd_ratio_list(:, :, 1), 1);
iso_str.fa = fun_analysis_get_basic_statistics_in_column(fractional_anisotropy, 1);

% Record simulation result for computing p-value
iso_str.data.fa = fractional_anisotropy;
iso_str.data.svd_min2max = svd_min2max_list;
iso_str.data.svd_1 = svd_ratio_list(:, :, 1);
end
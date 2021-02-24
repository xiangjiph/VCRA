function iso_str = fun_simulation_get_weighted_isotropy_stat(weight_vec, num_simu)
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
    str_template = fun_initialized_structure_array_with_fieldname_list({'num_simu', 'svd_ratio_mean', ...
        'svd_ratio_std', 'svd_min2max_mean', 'svd_min2max_std', ...
        'fractional_anisotropy_mean', 'fractional_anisotropy_std', 'data'});
end
iso_str = str_template;
if nargin < 2
    num_simu = 10000;
end
iso_str.num_simu = num_simu;

num_pts = numel(weight_vec);
if num_simu < 3
   % How to handle this case?  
   return;    
end
% Compute statistics in simulation 
svd_list = nan(3, num_simu);
for iter_trial = 1 : num_simu
    % Pair the weight with random orientation
    rand_ori = randn(num_pts, 3);
    rand_ori = rand_ori ./ vecnorm(rand_ori.').';
    rand_ori = bsxfun(@times, rand_ori, weight_vec);
    rand_ori = rand_ori + rand_ori .* (-2 * (rand_ori(:, 3) < 0));
%     rand_ori = cat(1, rand_ori, -rand_ori);
    cov_mat = (rand_ori.' * rand_ori) ./ (size(rand_ori, 1) - 1);
    
    svd_val = svd(cov_mat, 'econ');
    % Store 
    svd_list(:, iter_trial) = svd_val;
end
svd_sum = sum(svd_list, 1) .';
svd_norm2 = sum(svd_list .^2, 1) .';

svd_list = svd_list .';
svd_ratio_list = bsxfun(@rdivide, svd_list, svd_sum);
svd_min2max_list = svd_ratio_list(:, 3) ./ svd_ratio_list(:, 1);

fractional_anisotropy = sqrt(3 * sum((svd_list - svd_sum ./ 3) .^2, 2) ./ (2 * svd_norm2));

iso_str.svd_ratio_mean = mean(svd_ratio_list, 1);
iso_str.svd_ratio_std = std(svd_ratio_list, 0, 1);

iso_str.svd_min2max_mean = mean(svd_min2max_list);
iso_str.svd_min2max_std = std(svd_min2max_list);

iso_str.fractional_anisotropy_mean = mean(fractional_anisotropy);
iso_str.fractional_anisotropy_std = std(fractional_anisotropy);
% Record simulation result for computing p-value
iso_str.data.fa = fractional_anisotropy;
iso_str.data.svd_min2max = svd_min2max_list;
iso_str.data.svd_1 = svd_ratio_list(:, 1);
end
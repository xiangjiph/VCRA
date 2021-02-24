function iso_str = fun_simulation_get_isotropy_stat(num_pts, num_simu)
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
iso_str = fun_initialized_structure_array_with_fieldname_list({'svd_ratio_mean', ...
    'svd_ratio_std', 'svd_min2max_mean', 'svd_min2max_std'});
if nargin < 2
    num_simu = 50000;
end
if num_simu < 3
   % How to handle this case?  
   return;    
end
% Compute statistics in simulation 
svd_ratio_list = nan(3, num_simu);
svd_min2max_list = nan(num_simu, 1);
parfor iter_trial = 1 : num_simu
    % Pair the weight with random orientation
    rand_ori = randn(num_pts, 3);
    rand_ori = rand_ori ./ vecnorm(rand_ori.').';    
    rand_ori_flip = [rand_ori; -rand_ori];
    cov_mat = (rand_ori_flip.' * rand_ori_flip) ./ (2 * num_pts - 1);
    svd_val = svd(cov_mat, 'econ');
    svd_val_n = svd_val ./ sum(svd_val);
    % Store 
    svd_min2max_list(iter_trial) = svd_val_n(3) / svd_val_n(1);
    svd_ratio_list(:, iter_trial) = svd_val_n;
end
iso_str.svd_ratio_mean = mean(svd_ratio_list, 2)';
iso_str.svd_ratio_std = std(svd_ratio_list, 0, 2)';

iso_str.svd_min2max_mean = mean(svd_min2max_list);
iso_str.svd_min2max_std = std(svd_min2max_list);
end
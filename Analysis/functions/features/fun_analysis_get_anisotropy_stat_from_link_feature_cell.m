function info = fun_analysis_get_anisotropy_stat_from_link_feature_cell(link_feature, max_select_radius, vol_weighted_Q)
% fun_analysis_get_anisotropy_stat_from_link_feature_cell computes the
% regional vessel orientation anisotropy significance level. 
% Input: 
%   link_feature: table of link features, including dt_median, length,
%   in_bbox_volume, ep1_tp_ep2_direction_vec
%   radius_th: numerical nonnegative scalar, for selecting vessels for
%   anisotropy stat by radius. Vessel of radius not greater than this value
%   will be selected. Default value is inf.
%   vol_weighted_Q: logical scalar, scale the orientation vector by the
%   volume of the vector inside the local volume. ( well... actually it
%   should be the in-bounding-box endpoint to endpoint vecotr...). Deault
%   value is false.
%
% Output: 
%   info: structure defined below. 
%
% Implemented by Xiang Ji on 08/28/2019
   
%% Parameters
if nargin < 2
    max_select_radius = inf;
    vol_weighted_Q = false;
elseif nargin < 3
    vol_weighted_Q = false;
end
volume_wrighted_num_simulation = 1000;
% Not sure if making uni_ori_isotropy_info a persistent variable would
% cause problem in parallel computing or not
persistent uni_ori_isotropy_info
if ~vol_weighted_Q && isempty(uni_ori_isotropy_info)
    uni_ori_isotropy_info = load('./Metadata/uni_ori_isotropy.mat');
end
%%
grid_size = size(link_feature);
% Need to record both the orientation of the major SVD vector and the
% largest singular value.
[fa_cell, fa_z_cell, min2max_cell, svd_1_cell, min2max_z_cell, svd_1_z_cell, ori_vec_cell, num_data_cell] = deal(cell(grid_size(3), 1));
parfor iter_idx_3 = 1 : grid_size(3)
    maxNumCompThreads('automatic');
    tmp_tic = tic;
    fprintf('Processing grid layer %d\n', iter_idx_3);
    [tmp_fa, tmp_fa_z, tmp_min2max, tmp_svd_1, tmp_min2max_z, tmp_svd_1_z, tmp_num_data] = deal(nan(grid_size(1:2)));
    tmp_vec = nan([3, grid_size(1:2)]);
    tmp_data_cell = link_feature(:, :, iter_idx_3);
    for iter_idx_2 = 1 : grid_size(2)
        for iter_idx_1 = 1 : grid_size(1)
            tmp_link_table = tmp_data_cell{iter_idx_1, iter_idx_2};
            if ~isempty(tmp_link_table)
               if isfinite(max_select_radius)
                   tmp_radius = tmp_link_table.dt_median;
                   tmp_selected_Q = tmp_radius <= max_select_radius;
               else
                   tmp_selected_Q = true(size(tmp_link_table.dt_median));
               end
               
               if vol_weighted_Q
                   tmp_ep2ep_vec = bsxfun(@times, tmp_link_table.ep1_to_ep2_direction_vec, tmp_link_table.in_bbox_volume);
               else
                   tmp_ep2ep_vec = tmp_link_table.ep1_to_ep2_direction_vec;
               end
               
               tmp_isotropy_str = fun_analysis_get_link_anisotropy(tmp_ep2ep_vec(tmp_selected_Q, :), true);
               tmp_num_data(iter_idx_1, iter_idx_2) = tmp_isotropy_str.num_data;
               if ~isempty(tmp_isotropy_str.svd_min2max)
                   if ~vol_weighted_Q
                       tmp_uni_ori_mean = uni_ori_isotropy_info.min2max.mean(tmp_isotropy_str.num_data);
                       tmp_uni_ori_std = uni_ori_isotropy_info.min2max.std(tmp_isotropy_str.num_data);
                       
                       tmp_uni_ori_svd_mean = uni_ori_isotropy_info.svd_1.mean(tmp_isotropy_str.num_data);
                       tmp_uni_ori_svd_std = uni_ori_isotropy_info.svd_1.std(tmp_isotropy_str.num_data);
                       
                       tmp_uni_fa_mean = uni_ori_isotropy_info.fa.mean(tmp_isotropy_str.num_data);
                       tmp_uni_fa_std = uni_ori_isotropy_info.fa.std(tmp_isotropy_str.num_data);
                   else
                       tmp_uni_ori_str = fun_simulation_get_weighted_isotropy_stat(tmp_link_table.in_bbox_volume(tmp_selected_Q, :), volume_wrighted_num_simulation);
                       
                       tmp_uni_ori_mean = tmp_uni_ori_str.svd_min2max_mean;
                       tmp_uni_ori_std = tmp_uni_ori_str.svd_min2max_std;
                       
                       tmp_uni_ori_svd_mean = tmp_uni_ori_str.svd_ratio_mean(1);
                       tmp_uni_ori_svd_std = tmp_uni_ori_str.svd_ratio_std(1);
                       
                       tmp_uni_fa_mean = tmp_uni_ori_str.fractional_anisotropy_mean;
                       tmp_uni_fa_std = tmp_uni_ori_str.fractional_anisotropy_std;
                   end
                   
                   tmp_fa(iter_idx_1, iter_idx_2) = tmp_isotropy_str.fractional_anisotropy;
                   tmp_min2max(iter_idx_1, iter_idx_2) = tmp_isotropy_str.svd_min2max;
                   tmp_svd_1(iter_idx_1, iter_idx_2) = tmp_isotropy_str.svd_value_ratio(1);
                   
                   
                   tmp_min2max_z(iter_idx_1, iter_idx_2) = (tmp_isotropy_str.svd_min2max - tmp_uni_ori_mean) / tmp_uni_ori_std;
                   tmp_svd_1_z(iter_idx_1, iter_idx_2) = (tmp_isotropy_str.svd_value_ratio(1) - tmp_uni_ori_svd_mean) / tmp_uni_ori_svd_std;
                   tmp_fa_z(iter_idx_1, iter_idx_2) = (tmp_isotropy_str.fractional_anisotropy - tmp_uni_fa_mean) / tmp_uni_fa_std;
                   tmp_vec(:, iter_idx_1, iter_idx_2) = tmp_isotropy_str.svd_max_vec;
               end
            end
        end
    end
    num_data_cell{iter_idx_3} = tmp_num_data;
    fa_cell{iter_idx_3} = tmp_fa;
    min2max_cell{iter_idx_3} = tmp_min2max;
    svd_1_cell{iter_idx_3} = tmp_svd_1;
    
    min2max_z_cell{iter_idx_3} = tmp_min2max_z;
    svd_1_z_cell{iter_idx_3} = tmp_svd_1_z;
    fa_z_cell{iter_idx_3} = tmp_fa_z;
    
    ori_vec_cell{iter_idx_3} = tmp_vec;
    fprintf('Finish processing grid layer %d. Elapsed time is %f seconds\n', iter_idx_3, toc(tmp_tic));
end
info = struct;
info.grid_size = grid_size;
info.radius_th = max_select_radius;
info.weighted_Q = vol_weighted_Q;

info.num_data = cat(3, num_data_cell{:});
info.min2max = cat(3, min2max_cell{:});
info.svd_1 = cat(3, svd_1_cell{:});
info.fractional_anisotropy = cat(3, fa_cell{:});

info.min2max_z = cat(3, min2max_z_cell{:});
info.svd_1_z = cat(3, svd_1_z_cell{:});
info.fractional_anisotropy_z = cat(3, fa_z_cell{:});

info.ori_vec = cat(4, ori_vec_cell{:});

end
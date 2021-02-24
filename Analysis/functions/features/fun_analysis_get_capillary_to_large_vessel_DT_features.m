function DT_features = fun_analysis_get_capillary_to_large_vessel_DT_features(cc_ind, label_array, dt, dt_ind)
% fun_analysis_get_capillary_to_large_vessel_DT_features computes the
% noncapillary distance transform properties for link connected components
% Input: 
%   cc_ind: cell array. Each cell contains the linear indices of the
%   cneterline voxel of capillaries in the mask. 
%   label_array: 3D numerical array. Reconstucted labeled array generated
%   by fun_skeleton_reconstruction_label_aprox. 
%   [dt, dt_ind] = bwdist(label_array ~= 0);
% Output: 
%
% Implemented by Xiang Ji on 09/02/2019
num_cc = numel(cc_ind);
[DT_features.dist_to_nearest_noncapillary_ep1, DT_features.dist_to_nearest_noncapillary_ep2, ...
    DT_features.nearest_noncapillary_voxel_label_mode, ...
    DT_features.dist_to_nearest_noncapillary_mean, DT_features.dist_to_nearest_noncapillary_median] = deal(nan(num_cc, 1));
DT_features.nearest_noncapillary_voxel_label = cell(num_cc, 1);
for iter_cc = 1 : num_cc
    tmp_ind = cc_ind{iter_cc};
    tmp_dt = dt(tmp_ind);
    tmp_label = label_array(dt_ind(tmp_ind));
    DT_features.nearest_noncapillary_voxel_label{iter_cc} = tmp_label;
    DT_features.dist_to_nearest_noncapillary_ep1(iter_cc) = tmp_dt(1);
    DT_features.dist_to_nearest_noncapillary_ep2(iter_cc) = tmp_dt(end);
    DT_features.nearest_noncapillary_voxel_label_mode(iter_cc) = mode(tmp_label);
    DT_features.dist_to_nearest_noncapillary_mean(iter_cc) = mean(tmp_dt);
    DT_features.dist_to_nearest_noncapillary_median(iter_cc) = median(tmp_dt); 
end
end
function merged_stack_str = fun_analysis_radius_calibration_merge_registeration(merged_data_cell)

merged_data_cell = cat(1, merged_data_cell{:});



merge_field = {'Dist_fixed_2_moving', 'Fixed_est_corr', 'Fixed_ori_z', ...
    'Fixed_r_um', 'Moving_ori_z', 'Moving_r_um', };
merged_stack_str = struct;
for iter_field = 1 : numel(merge_field)
    tmp_field_name = merge_field{iter_field};
    merged_stack_str.voxel.(tmp_field_name) = cat(1, merged_data_cell.(tmp_field_name));
end

fixed_cc_str = cat(1, merged_data_cell.Fixed_cc_features);
merge_field = fieldnames(fixed_cc_str);
for iter_field = 1 : numel(merge_field)
    tmp_field_name = merge_field{iter_field};
    merged_stack_str.cc.Fixed.(tmp_field_name) = cat(1, fixed_cc_str.(tmp_field_name));
end

moving_cc_str = cat(1, merged_data_cell.Moving_cc_features);
merge_field = fieldnames(moving_cc_str);
for iter_field = 1 : numel(merge_field)
    tmp_field_name = merge_field{iter_field};
    merged_stack_str.cc.Moving.(tmp_field_name) = cat(1, moving_cc_str.(tmp_field_name));
end

merge_field = {'Dist_cc_mean', 'Dist_cc_median', 'Dist_cc_std', 'Dist_cc_cv'};
for iter_field = 1 : numel(merge_field)
    tmp_field_name = merge_field{iter_field};
    merged_stack_str.cc.(tmp_field_name) = cat(1, merged_data_cell.(tmp_field_name));
end
end
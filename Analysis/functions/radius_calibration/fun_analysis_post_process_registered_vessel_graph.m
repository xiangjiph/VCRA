function registered_str = fun_analysis_post_process_registered_vessel_graph(...
    registered_str, vg_fix, vg_mov)
persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
vg_fix.link.pos_r = cat(1, vg_fix.link.cc_r{:});
vg_mov.link.pos_r = cat(1, vg_mov.link.cc_r{:});
registered_str.Fixed_r_um = vg_fix.link.pos_r(registered_str.Fixed_matched_input_idx, :);
registered_str.Moving_r_um = vg_mov.link.pos_r(registered_str.Moving_matched_input_idx, :);

vg_fix.link.pos_local_ori_z = cat(1, vg_fix.link.cc_local_ori_vec_z{:});
vg_mov.link.pos_local_ori_z = cat(1, vg_mov.link.cc_local_ori_vec_z{:});
registered_str.Fixed_ori_z = abs(vg_fix.link.pos_local_ori_z(registered_str.Fixed_matched_input_idx));
registered_str.Moving_ori_z = abs(vg_mov.link.pos_local_ori_z(registered_str.Moving_matched_input_idx));

if strcmp(vg_fix.info.image_group, 'In_vivo')
   vg_fix.link.pos_est_corr = cat(1, vg_fix.link.cc_est_corr{:});
   registered_str.Fixed_est_corr = vg_fix.link.pos_est_corr(registered_str.Fixed_matched_input_idx);
end
%% Select matched pairs
[fixed_label_list_idx, registered_str.Fixed_cc_label_unique] = fun_bin_data_to_idx_list(registered_str.Fixed_cc_label);
unique_idx = cellfun(@(x) x(1), fixed_label_list_idx);
registered_str.Moving_cc_label_unique = registered_str.Moving_cc_label(unique_idx);
assert(numel(registered_str.Fixed_cc_label_unique) == numel(registered_str.Moving_cc_label_unique));

registered_str.Fixed_cc_features = fun_structure_field_slicing_by_index(vg_fix.link.features, registered_str.Fixed_cc_label_unique, 1);
registered_str.Moving_cc_features = fun_structure_field_slicing_by_index(vg_mov.link.features, registered_str.Moving_cc_label_unique, 1);
%% Distance between matched segments
registered_str.Dist_cc = cell(numel(fixed_label_list_idx), 1);
[registered_str.Dist_cc_mean, registered_str.Dist_cc_median, ...
    registered_str.Dist_cc_std] = deal(nan(numel(fixed_label_list_idx), 1));
for iter_cc = 1 : numel(fixed_label_list_idx)
    tmp_dist = registered_str.Dist_fixed_2_moving(fixed_label_list_idx{iter_cc});
    registered_str.Dist_cc{iter_cc} = tmp_dist;
    registered_str.Dist_cc_mean(iter_cc) = mean(tmp_dist);
    registered_str.Dist_cc_median(iter_cc) = median(tmp_dist);
    registered_str.Dist_cc_std(iter_cc) = std(tmp_dist);
end
registered_str.Dist_cc_cv = registered_str.Dist_cc_std ./ ...
    registered_str.Dist_cc_mean;
%% Information
registered_str.Moving_image_group = vg_mov.info.image_group;
registered_str.Fixed_image_group = vg_fix.info.image_group;

registered_str.Moving_grid_version = vg_mov.info.grid_version;
registered_str.Fixed_grid_version = vg_fix.info.grid_version;

assert(vg_mov.info.ROI_ID == vg_fix.info.ROI_ID, 'Mismatched ROI ID');
registered_str.ROI_ID = vg_mov.info.ROI_ID;
assert(strcmp(vg_mov.info.stack, vg_fix.info.stack), 'Mismatched image stack');
registered_str.stack = vg_mov.info.stack;
assert(strcmp(vg_mov.info.dataset_name, vg_fix.info.dataset_name), 'Mismatched dataset name');
registered_str.dataset_name = vg_mov.info.dataset_name;

registered_str.filepath = fullfile(DataManager.fp_analysis_data_folder(...
    registered_str.dataset_name, registered_str.stack), 'Matched_vessel', ...
    sprintf('%s_%s_%d_%s_%s_to_%s.mat', registered_str.dataset_name, ...
    registered_str.stack, registered_str.ROI_ID, registered_str.method, ...
    registered_str.Moving_image_group, registered_str.Fixed_image_group));
end
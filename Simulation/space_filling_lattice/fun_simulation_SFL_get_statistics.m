function stat_str = fun_simulation_SFL_get_statistics(valid_link_feature_table, target_link_length)

stat_str = struct;
tmp_stat_properties_name = valid_link_feature_table.Properties.VariableNames;
for iter_prop = 1 : numel(tmp_stat_properties_name)
    tmp_feature_name = tmp_stat_properties_name{iter_prop};
    tmp_data = valid_link_feature_table.(tmp_feature_name);
    stat_str.(sprintf('%s_mean', tmp_feature_name)) = mean(tmp_data, 1, 'omitnan');
    stat_str.(sprintf('%s_median', tmp_feature_name)) = median(tmp_data, 1, 'omitnan');
    stat_str.(sprintf('%s_std', tmp_feature_name)) = std(tmp_data, 1, 1, 'omitnan');
end
stat_str.tissue_dt_prctile99_mean = stat_str.tissue_dt_prctile_mean(25);
stat_str.tissue_dt_prctile99_median = stat_str.tissue_dt_prctile_median(25);

stat_str.ValidTotalLength = sum(valid_link_feature_table.length);
stat_str.ValidTotalEp2epDist = sum(valid_link_feature_table.ep2ep_dist);

stat_str.ValidTotalSpace = sum(valid_link_feature_table.nearest_tissue_volume + ...
    valid_link_feature_table.recon_mask_num_voxel);

stat_str.ValidTotalLengthDensity = stat_str.ValidTotalLength / stat_str.ValidTotalSpace;
stat_str.ValidTotalLengthDensity_m_mm3 = stat_str.ValidTotalLengthDensity .* 1e3;

stat_str.ValidTotalEp2epDistDensity = stat_str.ValidTotalEp2epDist / stat_str.ValidTotalSpace;
stat_str.ValidTotalEp2epDist_m_mm3 = stat_str.ValidTotalEp2epDistDensity .* 1e3;

stat_str.TargetTotalLength = numel(valid_link_feature_table.length) * target_link_length;
stat_str.TargetTotalLengthDensity = stat_str.TargetTotalLength / stat_str.ValidTotalSpace;
stat_str.TargetTotalLengthDist_m_mm3 = stat_str.TargetTotalLengthDensity .* 1e3;
end
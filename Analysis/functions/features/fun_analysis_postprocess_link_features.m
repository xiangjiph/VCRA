function link_feature_table = fun_analysis_postprocess_link_features(link_feature_table)

%% Add radius to the maximum tissue-vessel distance
link_feature_table.nearest_tissue_2_skl_dist_max = link_feature_table.nearest_tissue_dt_max + ...
    link_feature_table.dt_median;
% Compute the coefficient of variance for the radius 
link_feature_table.dt_cv = link_feature_table.dt_std ./ link_feature_table.dt_mean;
% Tortuosity
link_feature_table.tortuosity = 1 ./ link_feature_table.straightness;
%% Add neighbor volume sum 
cdf_order_edge = [1 : 1 : 6, inf];
num_bin = numel(cdf_order_edge) - 1;
num_link = size(link_feature_table, 1);
vol_r_in_bch_od = nan(num_bin, num_link);
for iter_link = 1 : num_link
    vol_r_in_bch_od(:, iter_link) = fun_analysis_dt_pt_sum_vol_r_by_branch_order(...
        link_feature_table.nb_lk_bch_od{iter_link}, link_feature_table.nb_lk_vol_r{iter_link}, cdf_order_edge);
end
link_feature_table.vol_r_in_bch_od = vol_r_in_bch_od';
link_feature_table.vol_r_in_bch_od_1 = link_feature_table.vol_r_in_bch_od(:, 1);
link_feature_table.vol_r_in_bch_od_gt_5 = link_feature_table.vol_r_in_bch_od(:, end);
%% Delete unused field
field_to_delete = {'link_com', 'cc_sub_pca1_vec', 'cc_sub_pca2_vec', 'cc_sub_pca3_vec', ...
    'cc_sub_cov_eig_val', 'shortest_loop_node_label', 'shortest_loop_link_label', 'nb_lk_label'};
for iter_field = 1 : numel(field_to_delete)
    tmp_field_name = field_to_delete{iter_field};
    if ismember(tmp_field_name, link_feature_table.Properties.VariableNames)
        link_feature_table.(tmp_field_name) = [];
    end
end
end
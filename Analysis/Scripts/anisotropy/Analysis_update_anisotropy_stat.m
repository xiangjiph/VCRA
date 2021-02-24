DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
% stack = 'ML_2019_01_24';
mask_version = '240_cube_recon';
skel_version = '240_cube_auto';
grid_version = '240_cube';
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
capillary_max_raidus = 3.5; % um;
im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
    'whole_brain_stat', 'internal_cubes', 'anisotropy');

registration_name = 'Allen_2017_25um_nonrigid.mat';
registration_str = DataManager.load_registration_data(dataset_name, stack, registration_name);
%% Compute the anisotropy data on the fly
fprintf('Update the anisotropy data\n');
[ai_all_cell, ai_cap_cell] = deal(cell(grid_info.num_valid_cube, 1));

[ai_all_lw_cell, ai_cap_lw_cell] = deal(cell(grid_info.num_valid_cube, 1));
p = gcp('nocreate');
delete(p);
parpool(12);
parfor iter_cube = 1 : grid_info.num_valid_cube
    % for iter_cube = 1 : grid_info.num_valid_cube
    maxNumCompThreads(4);
    tmp_tic = tic;
    tmp_grid_sub = grid_info.bbox_grid_sub_list(iter_cube, :);
    try
        tmp_skel = DataManager.load_block_skl(dataset_name, stack, skel_version, ...
            tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3));
        tmp_subgrid_graph = fun_skeleton_to_graph(tmp_skel.ind, tmp_skel.block_size);
        tmp_subgrid_graph.radius = sparse(double(tmp_skel.ind), 1, double(tmp_skel.r), prod(tmp_skel.block_size), 1);
    catch ME
        fprintf('Fail to load %s (%d, %d, %d)\n', skel_version, tmp_grid_sub);
        continue;
    end
    anisotropy_all_vw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
        [0, inf], 'volume');
    %     anisotropy_capillary_vw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
    %         [0, capillary_max_raidus], 'volume');
    anisotropy_all_lw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
        [0, inf], 'length');
    anisotropy_capillary_lw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
        [0, capillary_max_raidus], 'length');
    
    ai_all_cell{iter_cube} = anisotropy_all_vw;
    %     ai_cap_cell{iter_cube} = anisotropy_capillary_vw;
    ai_cap_lw_cell{iter_cube} = anisotropy_capillary_lw;
    ai_all_lw_cell{iter_cube} = anisotropy_all_lw;
    fprintf('Finish processing %d / %d. Elapse time is %f seconds.\n', ...
        iter_cube, grid_info.num_valid_cube, toc(tmp_tic));
end

[ai_all_array, ai_cap_array, ai_all_lw_array, ai_cap_lw_array] = deal(cell(grid_info.grid_size));
bbox_ind = sub2ind(grid_info.grid_size, grid_info.bbox_grid_sub_list(:, 1), ...
    grid_info.bbox_grid_sub_list(:, 2), grid_info.bbox_grid_sub_list(:, 3));

ai_cap_lw_array(bbox_ind) = ai_cap_lw_cell;
ai_all_lw_array(bbox_ind) = ai_all_lw_cell;
ai_all_array(bbox_ind) = ai_all_cell;
ai_cap_array(bbox_ind) = ai_cap_cell;
%% Find the 240-cubes that are completely inside the brain mask
wb_mask_ds_ratio = 16;
wb_mask = DataManager.load_data(sprintf('%s_mask.nii.gz', fullfile(DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration'),...
    sprintf('%s_%s_d16x_registration', dataset_name, stack)))) > 0;
is_internal_240_cube_ratio = fun_analysis_get_bbox_in_mask_vol_ratio(wb_mask, grid_info.bbox_xyz_mmxx_list ./ wb_mask_ds_ratio);
is_internal_240_cube_Q = is_internal_240_cube_ratio > 0;

internal_cube_ind = grid_info.bbox_grid_ind_list(is_internal_240_cube_ratio == 1 );
%% Save data
% anisotropy_data = struct;
% anisotropy_data.grid = grid_info;
% anisotropy_data.cube_grid_ind = grid_info.bbox_grid_ind_list;
% anisotropy_data.internal_cube_ratio = is_internal_240_cube_ratio;
% anisotropy_data.ai_vw_all = ai_all_array;
% anisotropy_data.ai_vw_cap = ai_cap_array;
% anisotropy_data.ai_lw_cap = ai_cap_lw_array;
% anisotropy_data.ai_lw_all = ai_all_lw_array;
% anisotropy_data.fp = fullfile(DataManager.fp_metadata_folder(dataset_name, stack), sprintf('%s_%s_volume_weighted_anisotropy_data.mat', ...
%     dataset_name, stack));
% save(anisotropy_data.fp, '-struct', 'anisotropy_data');
% fprintf('Finish saving anisotropy data to %s\n', anisotropy_data.fp);
anisotropy_data = load(fullfile(DataManager.fp_metadata_folder(dataset_name, stack), sprintf('%s_%s_volume_weighted_anisotropy_data.mat', ...
    dataset_name, stack)));
%% Capillary anisotropy - volume weighted
ai_cap_stat = struct;
ai_cap_stat.fa_p = cellfun(@(x) x.fa_p, ai_cap_array(internal_cube_ind));
ai_cap_stat.fa = cellfun(@(x) x.fractional_anisotropy, ai_cap_array(internal_cube_ind));
ai_cap_stat.svd_rmax = cellfun(@(x) x.svd_value_ratio(1), ai_cap_array(internal_cube_ind));
ai_cap_stat.fa_z = cellfun(@(x) x.fa_z, ai_cap_array(internal_cube_ind));
ai_cap_stat.svd_1_p = cellfun(@(x) x.svd_1_p, ai_cap_array(internal_cube_ind));

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1] .* 2;
ax_hdl_1 = subplot(1,3,1);
tmp_x = ai_cap_stat.fa;
histogram(ax_hdl_1, tmp_x, 'Normalization', 'pdf');
ax_hdl_1.XLabel.String = 'FA';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
box(ax_hdl_1, 'off');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)));

ax_hdl_2 = subplot(1,3,2);
tmp_x = ai_cap_stat.svd_rmax;
histogram(ax_hdl_2, tmp_x, 'Normalization', 'pdf');
ax_hdl_2.XLabel.String = 'Normalized PCV';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)));

ax_hdl_3 = subplot(1,3,3);
tmp_x = ai_cap_stat.fa_p;
tmp_x_bin = [0, 10.^(-4 : 0.5 : 0)];
histogram(ax_hdl_3, tmp_x, tmp_x_bin, 'Normalization', 'cdf');
ax_hdl_3.XScale = 'log';
ax_hdl_3.YLabel.String = 'CDF';
ax_hdl_3.XLabel.String = 'FA p-value';
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.FontSize = 14;
box(ax_hdl_3, 'off');
legend(ax_hdl_3, sprintf('# cubes: %d\np < 0.05: %d%%\np < 10^{-4}: %d%%', ...
    numel(ai_cap_stat.fa_p), round(100*nnz(ai_cap_stat.fa_p < 0.05)/numel(ai_cap_stat.fa_p)), ...
    round(100*nnz(ai_cap_stat.fa_p == 0) / numel(ai_cap_stat.fa_p))), 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'on');

fig_fp = fullfile(im_save_folder, sprintf('%s_%s_internal_cubes_volume_weighted_capillary_anisotropy_stat.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Capillary anisotropy - length weighted
is_not_empty_Q = cellfun(@(x) ~isempty(x.num_data), ai_cap_lw_array(internal_cube_ind));
internal_cube_ind = internal_cube_ind(is_not_empty_Q);
ai_cap_lw_stat = struct;
ai_cap_lw_stat.fa_p = cellfun(@(x) x.fa_p, ai_cap_lw_array(internal_cube_ind));
ai_cap_lw_stat.fa = cellfun(@(x) x.fractional_anisotropy, ai_cap_lw_array(internal_cube_ind));
ai_cap_lw_stat.svd_rmax = cellfun(@(x) x.svd_value_ratio(1), ai_cap_lw_array(internal_cube_ind));
ai_cap_lw_stat.fa_z = cellfun(@(x) x.fa_z, ai_cap_lw_array(internal_cube_ind));
ai_cap_lw_stat.svd_1_p = cellfun(@(x) x.svd_1_p, ai_cap_lw_array(internal_cube_ind));

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1] .* 2;
ax_hdl_1 = subplot(1,3,1);
tmp_x = ai_cap_lw_stat.fa;
histogram(ax_hdl_1, tmp_x, 'Normalization', 'pdf');
ax_hdl_1.XLabel.String = 'FA';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
box(ax_hdl_1, 'off');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)));

ax_hdl_2 = subplot(1,3,2);
tmp_x = ai_cap_lw_stat.svd_rmax;
histogram(ax_hdl_2, tmp_x, 'Normalization', 'pdf');
ax_hdl_2.XLabel.String = 'Normalized PCV';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)));

ax_hdl_3 = subplot(1,3,3);
tmp_x = ai_cap_lw_stat.fa_p;
tmp_x_bin = [0, 10.^(-4 : 0.5 : 0)];
histogram(ax_hdl_3, tmp_x, tmp_x_bin, 'Normalization', 'cdf');
ax_hdl_3.XScale = 'log';
ax_hdl_3.YLabel.String = 'CDF';
ax_hdl_3.XLabel.String = 'FA p-value';
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.FontSize = 14;
box(ax_hdl_3, 'off');
legend(ax_hdl_3, sprintf('# cubes: %d\np < 0.05: %d%%\np < 10^{-4}: %d%%', ...
    numel(ai_cap_lw_stat.fa_p), round(100*nnz(ai_cap_lw_stat.fa_p < 0.05)/numel(ai_cap_lw_stat.fa_p)), ...
    round(100*nnz(ai_cap_lw_stat.fa_p == 0) / numel(ai_cap_lw_stat.fa_p))), 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'on');

fig_fp = fullfile(im_save_folder, sprintf('%s_%s_internal_cubes_length_weighted_capillary_anisotropy_stat.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Vessel anisotropy - volume weighted
ai_all_stat = struct;
ai_all_stat.fa_p = cellfun(@(x) x.fa_p, ai_all_array(internal_cube_ind));
ai_all_stat.fa = cellfun(@(x) x.fractional_anisotropy, ai_all_array(internal_cube_ind));
ai_all_stat.svd_rmax = cellfun(@(x) x.svd_value_ratio(1), ai_all_array(internal_cube_ind));
ai_all_stat.fa_z = cellfun(@(x) x.fa_z, ai_all_array(internal_cube_ind));

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1] .* 2;
ax_hdl_1 = subplot(1,3,1);
tmp_x = ai_all_stat.fa;
histogram(ax_hdl_1, tmp_x, 'Normalization', 'pdf');
ax_hdl_1.XLabel.String = 'FA';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
box(ax_hdl_1, 'off');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)));

ax_hdl_2 = subplot(1,3,2);
tmp_x = ai_all_stat.svd_rmax;
histogram(ax_hdl_2, tmp_x, 'Normalization', 'pdf');
ax_hdl_2.XLabel.String = 'Normalized PCV';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)));

ax_hdl_3 = subplot(1,3,3);
tmp_x = ai_all_stat.fa_p;
tmp_x_bin = [0, 10.^(-4 : 0.5 : 0)];
histogram(ax_hdl_3, tmp_x, tmp_x_bin, 'Normalization', 'cdf');
ax_hdl_3.XScale = 'log';
ax_hdl_3.YLabel.String = 'CDF';
ax_hdl_3.XLabel.String = 'FA p-value';
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.FontSize = 14;
box(ax_hdl_3, 'off');
legend(ax_hdl_3, sprintf('# cubes: %d\np < 0.05: %d%%\np < 10^{-4}: %d%%', ...
    numel(ai_all_stat.fa_p), round(100*nnz(ai_all_stat.fa_p < 0.05)/numel(ai_all_stat.fa_p)), ...
    round(100*nnz(ai_all_stat.fa_p == 0) / numel(ai_all_stat.fa_p))), 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'on');

fig_fp = fullfile(im_save_folder, sprintf('%s_%s_internal_cubes_volume_weighted_vessel_anisotropy_stat.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Correlation between fa, fa_z
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, ai_cap_stat.fa, ai_cap_stat.fa_z, 'DisplayStyle', 'tile', 'Normalization', 'pdf');
ax_hdl.XLabel.String = 'Capillary FA';
ax_hdl.YLabel.String = 'Capillary FA z-score';
cb_hdl = colorbar(ax_hdl);
cb_hdl.Label.String = 'PDF';
ax_hdl.FontSize = 14;

hold(ax_hdl, 'on');
fit_data_x = ai_cap_stat.fa;
fit_data_y = ai_cap_stat.fa_z;
linear_fit_obj = fitlm(fit_data_x, fit_data_y);
plt_hdl = plot(ax_hdl, fit_data_x, linear_fit_obj.Coefficients.Estimate(1) + ...
    linear_fit_obj.Coefficients.Estimate(2) .* fit_data_x, 'LineWidth', 3);

ldg_hdl = legend(plt_hdl, sprintf('Slope: %.2e\\pm%.1e\nIntercept: %.2f\\pm%.1e\nR-squared: %.2f\nData size: %d', linear_fit_obj.Coefficients.Estimate(2), ...
    linear_fit_obj.Coefficients.SE(2), linear_fit_obj.Coefficients.Estimate(1), linear_fit_obj.Coefficients.SE(1), ...
    linear_fit_obj.Rsquared.Adjusted, linear_fit_obj.NumObservations), 'Location', 'northwest');

fig_fp = fullfile(im_save_folder, sprintf('%s_%s_internal_cubes_capillary_FA_z_vs_FA.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Correlation between fa_p and svd_1_p
% fig_hdl = figure;
% fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
% ax_hdl = axes(fig_hdl);
% tmp_x_bin = [0, 10.^(-4:0.5:0)];
% histogram2(ax_hdl, ai_cap_stat.fa_p, ai_cap_stat.svd_1_p, tmp_x_bin, tmp_x_bin, 'DisplayStyle', 'tile', 'Normalization', 'count');
% ax_hdl.XScale = 'log';
% ax_hdl.YScale = 'log';
% ax_hdl.XLabel.String = 'Capillary FA p-value';
% ax_hdl.YLabel.String = 'Capillary PCV p-value';
% ax_hdl.ColorScale = 'log';
% cb_hdl = colorbar(ax_hdl);
% cb_hdl.Label.String = 'PDF';
% ax_hdl.FontSize = 14;
% fig_fp = fullfile(im_save_folder, sprintf('%s_%s_internal_cubes_capillary_FA_z_vs_FA.png', dataset_name, stack));
% fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Regional anisotropy comparison
region_comparison_cell = cell(2, 0);
region_comparison_cell(:, end+1) = {'Isocortex_subregions', [500 453 1057 677 669 247 31 972 44 714 95 254 22 541 922 895]};
region_comparison_cell(:, end+1) = {'Olfactory_area_subregions', [507 151 159 589 814 961 639 631 788 566]};
region_comparison_cell(:, end+1) = {'Hippocampal_formation_subregions', [382 423 463 726 982 19 822]};
region_comparison_cell(:, end+1) = {'Striatum_subregions', [672 493 275 278]};
region_comparison_cell(:, end+1) = {'Thalamus_subregions', [864 856]};
region_comparison_cell(:, end+1) = {'Midbrain_sensory_related_subregions', [302 4 580 271 874 460]};
region_comparison_cell(:, end+1) = {'Midbrain_motor_related_subregions', [381 749 246 128 294 795 1100 616 214 35 975 115 757 231 66 75 58 615]};
region_comparison_cell(:, end+1) = {'Midbrain_behavior_state_related_subregions', [374 1052 165]};
region_comparison_cell(:, end+1) = {'Pons_subregions', [1132 987 1117]};
region_comparison_cell(:, end+1) = {'Medulla_subregions', [386 370 379]};
region_comparison_cell(:, end+1) = {'Cerebellar_cortex_subregions', [1144 1145 1143 645 1073]};
region_comparison_cell(:, end+1) = {'Major_regions', [315 698 1089 703 477 803 549 1097 313 771 354 528]};
num_region_comparison_cell = size(region_comparison_cell, 2);
extract_field_name = {'fractional_anisotropy', 'fa_z', 'fa_p', 'svd_1_z', 'svd_1_p', 'svd_value_ratio'};
for iter_plot = 3 : num_region_comparison_cell
    % vis_name = 'Isocortex_subregions';
    % analysis_atlas_id = [500 453 1057 677 669 247 31 972 44 714 95 254 22 541 922 895];
    vis_name = region_comparison_cell{1, iter_plot};
    analysis_atlas_id = region_comparison_cell{2, iter_plot};
    num_region = numel(analysis_atlas_id);
    region_stat_cell = cell(num_region, 1);
    for iter_region = 1 : num_region
        tmp_tic = tic;
        vis_atlas_id = analysis_atlas_id(iter_region);
        tmp_region_stat = struct;
        
        tmp_region_stat.region_name = registration_str.structure_table.name{registration_str.id_2_ind(vis_atlas_id)};
        tmp_region_stat.region_name_abv = registration_str.structure_table.acronym{registration_str.id_2_ind(vis_atlas_id)};
        tmp_region_stat.has_valid_cubeQ = false;
        tmp_region_stat.num_valid_cube = nan;
        
        [subregion_cc, subregion_mask] = fun_registration_get_region_cc(registration_str, vis_atlas_id);
        in_cc_volume_fraction = fun_analysis_get_bbox_in_mask_vol_ratio(subregion_mask > 0, grid_info.bbox_xyz_mmxx_list ./ ...
            [registration_str.info.voxel_size_um, registration_str.info.voxel_size_um]);
        
        valid_cube_ind = anisotropy_data.cube_grid_ind(in_cc_volume_fraction > 0.5 & is_internal_240_cube_ratio == 1);
        
        tmp_is_valid_Q = cellfun(@(x) isfield(x, 'fractional_anisotropy'), anisotropy_data.ai_vw_all(valid_cube_ind));
        tmp_is_valid_Q = tmp_is_valid_Q & cellfun(@(x) isfield(x, 'fractional_anisotropy'), anisotropy_data.ai_lw_cap(valid_cube_ind));
        valid_cube_ind = valid_cube_ind(tmp_is_valid_Q);
        
        tmp_region_stat.valid_cube_ind = valid_cube_ind;
        tmp_region_stat.num_valid_cube = numel(valid_cube_ind);
        tmp_region_stat.has_valid_cubeQ = tmp_region_stat.num_valid_cube > 0;
        
        if tmp_region_stat.has_valid_cubeQ            
            region_ai_all_stat = struct;
            region_ai_cap_stat = struct;
            for iter_field = 1 : numel(extract_field_name)
                tmp_field_name = extract_field_name{iter_field};
                region_ai_all_stat.(tmp_field_name) = cellfun(@(x) x.(tmp_field_name)(1), ...
                    anisotropy_data.ai_vw_all(valid_cube_ind));
                region_ai_cap_stat.(tmp_field_name) = cellfun(@(x) x.(tmp_field_name)(1), ...
                    anisotropy_data.ai_lw_cap(valid_cube_ind));
            end
            tmp_region_stat.ai_cap = region_ai_cap_stat;
            tmp_region_stat.ai_all = region_ai_all_stat;
        else
            fprintf('Structure %s does not have internal valid cube\n', tmp_region_stat.region_name);
        end
        region_stat_cell{iter_region} = tmp_region_stat;
        fprintf('Finish collecting data for region %s. Elapsed time is %f seconds\n', tmp_region_stat.region_name, ...
            toc(tmp_tic));
    end
    is_valid_region_Q = cellfun(@(x) x.has_valid_cubeQ, region_stat_cell);
    region_stat_cell = region_stat_cell(is_valid_region_Q);
    num_region = numel(region_stat_cell);
    
    %% Generate box plot
    region_name_list = cellfun(@(x) x.region_name, region_stat_cell, 'UniformOutput', false);
    region_name_abv_list = cellfun(@(x) x.region_name_abv, region_stat_cell, 'UniformOutput', false);
    
    regional_fa_z = cellfun(@(x) x.ai_cap.fa_z, region_stat_cell, 'UniformOutput', false);
    regional_fa = cellfun(@(x) x.ai_cap.fractional_anisotropy, region_stat_cell, 'UniformOutput', false);
    regional_svd_1 = cellfun(@(x) x.ai_cap.svd_value_ratio, region_stat_cell, 'UniformOutput', false);
    
    
    tmp_region_sort_val = cellfun(@median, regional_fa_z);
    [~, tmp_region_sort_ind] = sort(tmp_region_sort_val, 'ascend');
    regional_fa_z = regional_fa_z(tmp_region_sort_ind);
    regional_fa = regional_fa(tmp_region_sort_ind);
    regional_svd_1 = regional_svd_1(tmp_region_sort_ind);
    
    region_name_list = region_name_list(tmp_region_sort_ind);
    region_name_abv_list = region_name_abv_list(tmp_region_sort_ind);
    
    tmp_num_data = cellfun(@numel, regional_fa_z);
    tmp_label = repelem((1 : num_region)', tmp_num_data);
    
    
    fig_hdl = figure('Visible', 'off');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1,3] .* 2;
    
    ax_hdl_1 = subplot(3,1,1);
    tmp_box_data = cat(1, regional_fa{:});
    boxplot(ax_hdl_1, tmp_box_data, tmp_label);
    ax_hdl_1.XTickLabel = region_name_abv_list;
    ax_hdl_1.XTickLabelRotation = 90;
    ax_hdl_1.YLabel.String = 'Fractional Anisotropy';
    ax_hdl_1.FontSize = 14;
    ax_hdl_1.YLim = [0,1];
    
    ax_hdl_2 = subplot(3,1,2);
    tmp_box_data = cat(1, regional_fa_z{:});
    boxplot(ax_hdl_2, tmp_box_data, tmp_label);
    ax_hdl_2.XTickLabel = region_name_abv_list;
    ax_hdl_2.XTickLabelRotation = 90;
    ax_hdl_2.YLabel.String = 'FA z-score';
    ax_hdl_2.FontSize = 14;
    
    ax_hdl_3 = subplot(3,1,3);
    tmp_box_data = cat(1, regional_svd_1{:});
    boxplot(ax_hdl_3, tmp_box_data, tmp_label);
    ax_hdl_3.XTickLabel = region_name_abv_list;
    ax_hdl_3.XTickLabelRotation = 90;
    ax_hdl_3.YLabel.String = 'PCV';
    ax_hdl_3.FontSize = 14;
    ax_hdl_3.YLim = [0,1];
    
    fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_anisotropy_metric_boxplot.png', ...
        dataset_name, stack, vis_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% Statistical test if two data are from the same distribution
    fa_z_pairwise_pvalue = fun_analysis_get_pairwise_KStest_pvalue_matrix(regional_fa_z);
    fa_pairwise_pvalue = fun_analysis_get_pairwise_KStest_pvalue_matrix(regional_fa);
    svd_1_pairwise_pvalue = fun_analysis_get_pairwise_KStest_pvalue_matrix(regional_svd_1);
    
    fig_hdl = figure('Visible', 'off');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1.5];
    ax_hdl_1 = subplot(1,3,1);
    imagesc(fa_z_pairwise_pvalue);
    ax_hdl_1.DataAspectRatio = [1,1,1];
    ax_hdl_1.ColorScale = 'log';
    ax_hdl_1.XTick = 1 : numel(region_name_list);
    ax_hdl_1.XTickLabel = region_name_abv_list;
    ax_hdl_1.XTickLabelRotation = 90;
    ax_hdl_1.YTick = 1 : numel(region_name_list);
    ax_hdl_1.YTickLabel = region_name_list;
    cbar_1_hdl = colorbar(ax_hdl_1, 'Location', 'eastoutside');
    cbar_1_hdl.Label.String = 'p-Value';
    ax_hdl_1.Title.String = 'FA z-score';
    
    ax_hdl_2 = subplot(1,3,2);
    imagesc(fa_pairwise_pvalue)
    ax_hdl_2.DataAspectRatio = [1,1,1];
    ax_hdl_2.ColorScale = 'log';
    ax_hdl_2.XTick = 1 : numel(region_name_list);
    ax_hdl_2.XTickLabel = region_name_abv_list;
    ax_hdl_2.XTickLabelRotation = 90;
    ax_hdl_2.Title.String = 'FA';
    ax_hdl_2.YAxis.Visible = 'off';
    cbar_2_hdl = colorbar(ax_hdl_2, 'Location', 'eastoutside');
    cbar_2_hdl.Label.String = 'p-Value';
    
    
    ax_hdl_3 = subplot(1,3,3);
    imagesc(svd_1_pairwise_pvalue);
    ax_hdl_3.DataAspectRatio = [1,1,1];
    ax_hdl_3.ColorScale = 'log';
    ax_hdl_3.XTick = 1 : numel(region_name_list);
    ax_hdl_3.XTickLabel = region_name_abv_list;
    ax_hdl_3.XTickLabelRotation = 90;
    ax_hdl_3.Title.String = 'PCV 1';
    ax_hdl_3.YAxis.Visible = 'off';
    cbar_3_hdl = colorbar(ax_hdl_3, 'Location', 'eastoutside');
    cbar_3_hdl.Label.String = 'p-Value';
    
    fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_anisotropy_metric_distribution_kstest_pvalue_matrix.png',...
        dataset_name, stack, vis_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
end
%% Generate figures for illustrating anisotropy
% load_skel_ver = skel_version;
load_skel_ver = '240_cube_re';
vis_atlas_name = strjoin(registration_str.structure_table.name(full(registration_str.id_2_ind(vis_atlas_id))), '_n_');
vis_sub = fun_ind2sub(grid_info.grid_size, valid_cube_ind(143));
vis_recon_str = DataManager.load_block_mask(dataset_name, stack, ...
    mask_version, vis_sub(1), vis_sub(2), vis_sub(3));
vis_recon = fun_reconstruct_block_mask(vis_recon_str);

tmp_skel = DataManager.load_block_skl(dataset_name, stack, load_skel_ver, ...
    vis_sub(1), vis_sub(2), vis_sub(3));
tmp_subgrid_graph = fun_skeleton_to_graph(tmp_skel.ind, tmp_skel.block_size);
tmp_subgrid_graph.radius = sparse(double(tmp_skel.ind), 1, double(tmp_skel.r), prod(tmp_skel.block_size), 1);
anisotropy_capillary_lw = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
    [0, capillary_max_raidus], 'length');

lf.length = nan(tmp_subgrid_graph.link.num_cc, 1);
lf.ep1_to_ep2_direction_vec = nan(3, tmp_subgrid_graph.link.num_cc);
lf.dt_median = nan(tmp_subgrid_graph.link.num_cc, 1);
[ep_sub_1, ep_sub_2] = deal(nan(tmp_subgrid_graph.link.num_cc, 3));
for iter_link = 1 : tmp_subgrid_graph.link.num_cc
    if tmp_subgrid_graph.link.num_voxel_per_cc(iter_link) > 1
        tmp_ind = tmp_subgrid_graph.link.cc_ind{iter_link};
        tmp_sub = fun_ind2sub(tmp_subgrid_graph.num.mask_size, tmp_ind);
        ep_sub_1(iter_link, :) = tmp_sub(1, :);
        ep_sub_2(iter_link, :) = tmp_sub(end, :);
        tmp_ep1_to_ep2_voxel_vec = tmp_sub(end, :) - tmp_sub(1, :);
        tmp_ep1_to_ep2_vec_norm = sqrt(sum((tmp_ep1_to_ep2_voxel_vec ).^2));
        lf.ep1_to_ep2_direction_vec(:, iter_link) = tmp_ep1_to_ep2_voxel_vec ./ ...
            tmp_ep1_to_ep2_vec_norm;
        tmp_dt = full(tmp_subgrid_graph.radius(tmp_ind));
        lf.dt_median(iter_link) = median(tmp_dt(tmp_dt>0));
        lf.length(iter_link) = fun_graph_sub_to_length(tmp_sub, 1);
    end
end
ep_sub = cat(2, ep_sub_1, ep_sub_2);

lf.ep1_to_ep2_direction_vec = lf.ep1_to_ep2_direction_vec.';
selected_r_range = [0, 3.5];
selected_Q = lf.dt_median >= selected_r_range(1) & ...
    lf.dt_median <= selected_r_range(2);
weight_vec = lf.length(selected_Q);
ep2ep_vec = bsxfun(@times, lf.ep1_to_ep2_direction_vec(selected_Q, :), ...
    lf.length(selected_Q));
anisotropy_stat = fun_analysis_get_link_anisotropy(ep2ep_vec, true);
tmp_uni_ori_str = fun_simulation_get_weighted_isotropy_stat(weight_vec, 1e4);

lw_anisotropy_stat = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
    selected_r_range, 'length');
vw_anisotropy_stat = fun_analysis_get_anisotropy_stat_from_vessel_graph(tmp_subgrid_graph, ...
    selected_r_range, 'volume');
% volumeViewer(vis_recon)
vis_mask_smooth = imclose(vis_recon, strel('sphere', 3));
fv = isosurface(vis_mask_smooth);
fv = reducepatch(fv, 0.1);
%%
vis_agl = [120, 16];

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,3];

ax_hdl_1 = subplot(3,2,1);
patch_hdl = patch(ax_hdl_1, 'Vertices', fv.vertices, 'Faces', fv.faces, 'FaceColor', 'g', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.8, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
patch_hdl.FaceColor = [1, 0 , 0] .* 0.9;
ax_hdl_1.DataAspectRatio = [1,1,1];
ax_hdl_1.View = vis_agl;
camlight(ax_hdl_1);
ax_hdl_1.XLabel.String = 'X (\mum)';
ax_hdl_1.YLabel.String = 'Y (\mum)';
ax_hdl_1.ZLabel.String = 'Z (\mum)';
ax_hdl_1.Title.String = 'Reconstructed vessel network';

ax_hdl_2 = subplot(3,2,2);
for iter_link = 1 : tmp_subgrid_graph.link.num_cc
    line(ax_hdl_2, ep_sub(iter_link, [2,5]), ep_sub(iter_link, [1,4]), ep_sub(iter_link, [3,6]), 'LineWidth', 2);
    hold(ax_hdl_2, 'on');
end
ax_hdl_2.View = vis_agl;
ax_hdl_2.XLabel.String = 'X (\mum)';
ax_hdl_2.YLabel.String = 'Y (\mum)';
ax_hdl_2.ZLabel.String = 'Z (\mum)';
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.Title.String = 'Endpoint to endpoint connections';


ax_hdl_3 = subplot(3,2,3);
vis_vec = lf.ep1_to_ep2_direction_vec;
vis_vec = cat(1, vis_vec, -vis_vec);
scatter3(ax_hdl_3, vis_vec(:, 2), vis_vec(:, 1), vis_vec(:, 3), '+');
ax_hdl_3.DataAspectRatio = [1,1,1];
ax_hdl_3.Title.String = 'Vessel orientation';
ax_hdl_3.XLabel.String = 'X-component';
ax_hdl_3.YLabel.String = 'Y-component';
ax_hdl_3.ZLabel.String = 'Z-component';
ax_hdl_3.View = vis_agl;

rand_ori = randn(tmp_subgrid_graph.link.num_cc, 3);
rand_ori = rand_ori ./ vecnorm(rand_ori.').';
rand_ori_flip = [rand_ori; -rand_ori];

ax_hdl_4 = subplot(3,2,4);
scatter3(ax_hdl_4, rand_ori_flip(:, 2), rand_ori_flip(:, 1), rand_ori_flip(:, 3), '+');
ax_hdl_4.DataAspectRatio = [1,1,1];
ax_hdl_4.Title.String = 'Isotropic orientation';
ax_hdl_4.XLabel.String = 'X-component';
ax_hdl_4.YLabel.String = 'Y-component';
ax_hdl_4.ZLabel.String = 'Z-component';
ax_hdl_4.View = vis_agl;

ax_hdl_5 = subplot(3,2,5:6);
hist_hdl = histogram(ax_hdl_5, tmp_uni_ori_str.data.fa, 'Normalization', 'probability');
ax_hdl_5.XLim = [0,1];
ax_hdl_5.XLabel.String = 'Anisotropy';
ax_hdl_5.YLabel.String = 'Probability';
ax_hdl_5.Title.String = 'Anisotropy statistical test';
ax_hdl_5.YScale = 'log';
hold(ax_hdl_5, 'on');
line_hdl = line(ax_hdl_5, repelem(lw_anisotropy_stat.fractional_anisotropy, 2,1), ax_hdl_5.YLim);
line_hdl.LineWidth = 2;
line_hdl.Color = [0.8, 0, 0];
leg_hdl = legend(ax_hdl_5, 'Null hypothesis', sprintf('Capillary\nAnisotropy:%.2f\nZ-score:%.2f\n', ...
    lw_anisotropy_stat.fractional_anisotropy, lw_anisotropy_stat.fa_z));


fig_fp = fullfile(im_save_folder, sprintf('%s_%s_240_cube_%d_%d_%d_recon_in_%s.png', dataset_name, stack, ...
    vis_sub(1), vis_sub(2), vis_sub(3), strrep(vis_atlas_name, ' ', '_')));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Determine the major anatomical region
% All the major regions in the cerebrum
registration_str.structure_table.id(any(registration_str.structure_id_path_mat == 519, 2) & ...
    registration_str.structure_id_depth == 4, :)

% registration_str.structure_table(registration_str.fun_name_to_ind(registration_str.structure_table.name, 'vestibular'), :)
% registration_str.structure_table(registration_str.id_2_ind(370), :)
% 1. Do the neocortical regions have similar statistics?
%% 3D Visualization of regional anisotropy
extract_field_name = {'fractional_anisotropy', 'fa_z', 'fa_p', 'svd_1_z', 'svd_1_p', 'svd_value_ratio', 'svd_max_vec'};
vis_atlas_id = [1117];
vis_atlas_name = strjoin(registration_str.structure_table.name(full(registration_str.id_2_ind(vis_atlas_id))), '_n_');
[subregion_cc, subregion_mask] = fun_registration_get_region_cc(registration_str, vis_atlas_id);

% Determine the regional bounding box at 16 um resolution
vis_region_ind = cat(2, subregion_cc.PixelIdxList{:});
vis_region_sub = fun_ind2sub(size(subregion_mask), vis_region_ind);
vis_region_bbox_mm = min(vis_region_sub, [], 1);
vis_region_bbox_xx = max(vis_region_sub, [], 1);
vis_region_bbox_mmll = [vis_region_bbox_mm, vis_region_bbox_xx - vis_region_bbox_mm + 1];
vis_local_mask = crop_bbox3(subregion_mask, vis_region_bbox_mmll);
% 3D rendering of the mask
vis_mask_smooth = imclose(vis_local_mask > 0, strel('sphere', 3));
fv = isosurface(vis_mask_smooth);
fv = reducepatch(fv, 0.1);
% Compute the vector to be plot
in_cc_volume_fraction = fun_analysis_get_bbox_in_mask_vol_ratio(subregion_mask >0, grid_info.bbox_xyz_mmxx_list ./ ...
    [registration_str.info.voxel_size_um, registration_str.info.voxel_size_um]);
valid_cube_Q = in_cc_volume_fraction > 0.1 & is_internal_240_cube_ratio == 1;
valid_cube_ind = anisotropy_data.cube_grid_ind(valid_cube_Q);
% Compute the position of the 240 cube in the downsampled image
valid_cube_center_sub = grid_info.bbox_xyz_mmxx_list(valid_cube_Q, :);
valid_cube_center_sub = (valid_cube_center_sub(:, 1:3) + valid_cube_center_sub(:, 4:6)) ./ 2;
valid_cube_center_sub_rz = valid_cube_center_sub ./ registration_str.info.voxel_size_um;
valid_cube_ori = cellfun(@(x) x.svd_max_vec, anisotropy_data.ai_lw_cap(valid_cube_ind), 'UniformOutput', false);

region_ai_vw_all_stat = struct;
region_ai_lw_cap_stat = struct;
for iter_field = 1 : numel(extract_field_name)
    tmp_field_name = extract_field_name{iter_field};
    tmp_cell = cellfun(@(x) x.(tmp_field_name), anisotropy_data.ai_vw_all(valid_cube_ind), 'UniformOutput', false);
    tmp_cell = cat(1, tmp_cell{:});
    region_ai_vw_all_stat.(tmp_field_name) = tmp_cell;
    
    tmp_cell = cellfun(@(x) x.(tmp_field_name), anisotropy_data.ai_lw_cap(valid_cube_ind), 'UniformOutput', false);
    tmp_cell = cat(1, tmp_cell{:});
    region_ai_lw_cap_stat.(tmp_field_name) = tmp_cell;
end
region_ai_vw_all_stat.svd_max_vec = reshape(region_ai_vw_all_stat.svd_max_vec, 3, [])';
region_ai_lw_cap_stat.svd_max_vec = reshape(region_ai_lw_cap_stat.svd_max_vec, 3, [])';
region_ai_lw_cap_stat.svd_max_vec(region_ai_lw_cap_stat.svd_max_vec(:, 3) < 0, :) = ...
    - region_ai_lw_cap_stat.svd_max_vec(region_ai_lw_cap_stat.svd_max_vec(:, 3) < 0, :);
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
patch_hdl = patch(ax_hdl, 'Vertices', fv.vertices, 'Faces', fv.faces, 'FaceColor', 'g', ...
    'EdgeColor', 'none', 'FaceAlpha', 0.1, 'FaceLighting', 'gouraud', 'AmbientStrength', 0.2);
patch_hdl.FaceColor = [1, 1, 1] .* 0.9;
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.View = [30, 60];
camlight(ax_hdl);
hold(ax_hdl, 'on');
vis_cube_center_local_sub = valid_cube_center_sub_rz - vis_region_bbox_mmll(1:3);
vis_vec = region_ai_lw_cap_stat.svd_max_vec;
vis_vec = vis_vec .* region_ai_lw_cap_stat.fa_z;
q_hdl = quiver3(ax_hdl, vis_cube_center_local_sub(:, 2), vis_cube_center_local_sub(:, 1), vis_cube_center_local_sub(:, 3), ...
    vis_vec(:, 2), vis_vec(:, 1), vis_vec(:, 3));
q_hdl.LineWidth = 1;
q_hdl.ShowArrowHead = 'off';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.ZAxis.Direction = 'reverse';
ax_hdl.Title.String = strrep(vis_atlas_name, '_n_', ' \n');
ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%d'), ax_hdl.XTick .* registration_str.info.voxel_size_um(1), 'UniformOutput', false);
ax_hdl.YTickLabel = arrayfun(@(x) num2str(x, '%d'), ax_hdl.YTick .* registration_str.info.voxel_size_um(2), 'UniformOutput', false);
ax_hdl.ZTickLabel = arrayfun(@(x) num2str(x, '%d'), ax_hdl.ZTick .* registration_str.info.voxel_size_um(3), 'UniformOutput', false);

fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_lw_cap_ori_vis_3d.png', dataset_name, stack, strrep(vis_atlas_name, ' ', '_')));

% fun_print_image_in_several_formats(fig_hdl, fig_fp);
% delete(fig_hdl);
%% 3D Visualization in max projection reconstruction
recon_mask_proj_str = fun_analysis_load_whole_brain_mask_max_proj(grid_info, mask_version, grid_info.layer);
[ai_cap_lw_vec] = deal(nan([3, grid_info.grid_size]));
[ai_cap_lw_fa, ai_cap_lw_fa_z, ai_cap_lw_svd1, ai_cap_lw_num_seg] = deal(nan(grid_info.grid_size));
tmp_tic = tic;
for iter_cube = 1 : grid_info.num_valid_cube
    tmp_grid_ind = grid_info.bbox_grid_ind_list(iter_cube);
    tmp_grid_sub = grid_info.bbox_grid_sub_list(iter_cube, :);
    tmp_cap_ai_str = ai_cap_lw_cell{iter_cube};
    if ~isempty(tmp_cap_ai_str.svd_max_vec)
        ai_cap_lw_vec(:, tmp_grid_sub(1), tmp_grid_sub(2), tmp_grid_sub(3)) = ...
            tmp_cap_ai_str.svd_max_vec;
        ai_cap_lw_fa(tmp_grid_ind) = tmp_cap_ai_str.fractional_anisotropy;
        ai_cap_lw_fa_z(tmp_grid_ind) = tmp_cap_ai_str.fa_z;
        ai_cap_lw_svd1(tmp_grid_ind) = tmp_cap_ai_str.svd_value_ratio(1);
        ai_cap_lw_num_seg(tmp_grid_ind) = tmp_cap_ai_str.num_data;
    end
end
fprintf('Finish reorganizing the anisotropy data into numerical array. Elapsed time is %f seconds\n', ...
    toc(tmp_tic));

plot_info_cell = cell(4, 0);
plot_info_cell(:, end+1) = {ai_cap_lw_fa, 'Length_weighted_capillary_fractional_anisotropy', 'Fractional anisotropy', ai_cap_lw_vec};
plot_info_cell(:, end+1) = {ai_cap_lw_fa_z, 'Length_weighted_capillary_fractional_anisotropy_z_score', 'Z score', ai_cap_lw_vec};
plot_info_cell(:, end+1) = {ai_cap_lw_svd1, 'Length_weighted_capillary_PCV_1', 'PCV-1', ai_cap_lw_vec};

valid_array_Q = ai_cap_lw_num_seg > 50;
for iter_dir = 1 : 3
    fun_vis_whole_brain_local_network_orientation(plot_info_cell, recon_mask_proj_str, ...
        grid_info, valid_array_Q, iter_dir);
end
%%
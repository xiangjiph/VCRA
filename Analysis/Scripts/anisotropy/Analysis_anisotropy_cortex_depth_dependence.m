set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
% stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
stack_list = {'ML_2018_08_15', 'ML20190124'};
num_stack = numel(stack_list);
switch num_stack
    case 2
        merged_stack_name = 'BrainAB';
    case 3
        merged_stack_name = 'all_stack';
end

stack_name_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
grid_version = '240_cube';
data_folder = 'Radius_dependence_anisotropy';
skel_version = '240_cube_rec';
reconstruction_version = '240_cube_recon_sc';
linear_scaling_factor = 1.0521;

wb_stat_folder_name = 'whole_brain_stat_sc';
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_landmark.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');
vis_folder_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, merged_stack_name), ...
    'Radius_dependence_anisotropy');
%% Load data
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    tic;
    stack = stack_list{iter_stack};
    wb_region_stat_data = fun_analysis_load_whole_brain_data_for_regional_analysis(...
        dataset_name, stack, grid_version, reconstruction_version, ...
        skel_version, registration_version, wb_stat_folder_name, false);
    wb_region_stat_data.cube_stat.num_link = wb_region_stat_data.cube_stat.link_all_stat.num_data.length;
    wb_region_stat_data.cube_stat.num_cap = wb_region_stat_data.cube_stat.link_cap_stat.num_data.length;
    wb_region_stat_data.rda_data = DataManager.load_data(DataManager.fp_analysis_data_file(...
        dataset_name, stack, sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, ...
        stack, skel_version, data_folder)));
    %% Compute DT gradient and local penetrating direction 
    grid_bbox_mmxx_ds = ceil(wb_region_stat_data.grid_info.bbox_xyz_mmxx_list ./ ...
        wb_region_stat_data.brain_mask.ds_ratio);
    wb_region_stat_data.cube_stat.dt_ori_pca = fun_analysis_cube_dt_grad_pca(...
        wb_region_stat_data.brain_mask.dt .* linear_scaling_factor, grid_bbox_mmxx_ds);    
    %%
    wb_data_cell{iter_stack} = wb_region_stat_data;
    fprintf('Finish loading data in %s. Elapsed time is %f seconds\n', stack, toc);
end
%% data extraction setting 
allen_atlas_map_old_id_to_new_id = wb_data_cell{1}.registration.map_oldID_to_newID;
de_opt_str = struct;
de_opt_str.cube_stat_Q = true;
de_opt_str.node_feature_Q = false;
de_opt_str.link_feature_Q = false;
de_opt_str.vessel_graph_Q = false;
de_opt_str.depth_dependent_density_Q = false;
de_opt_str.merge_cc_Q = false;
de_opt_str.save_Q = false;
de_opt_str.save_fp = [];
if de_opt_str.merge_cc_Q
    num_cc = 1;
else
    num_cc = 2;
    cc_label = {'Left', 'Right'};
end
%% Extracted region information
parent_str_name = 'Cerebral_cortex_transition';
analysis_region_id_list = [453 247 669 985 895 31 961];

num_region = numel(analysis_region_id_list);
analysis_region_name = allen_atlas.structure_table.name(full(allen_atlas.id_2_ind(analysis_region_id_list)));
analysis_region_name_abbrv = allen_atlas.structure_table.acronym(full(allen_atlas.id_2_ind(analysis_region_id_list)));
vis_subfolder_name = 'Depth_dependence';
visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merged_stack_name), vis_subfolder_name, parent_str_name);

analysis_region_id_list = full(allen_atlas_map_old_id_to_new_id(analysis_region_id_list));
%% Collect all the region statistics in all the stack
region_data = cell(num_cc, num_stack, num_region);
collect_data_tic = tic;
for iter_region = 1 : num_region
    for iter_stack = 1 : num_stack
        tmp_tic = tic;
        tmp_region_id = analysis_region_id_list(iter_region);
        tmp_region_name = strrep(analysis_region_name{iter_region}, ' ', '_');
        tmp_wb_data = wb_data_cell{iter_stack};
        tmp_region_data = fun_analysis_get_atlas_regional_data(tmp_wb_data, ...
            tmp_region_id, tmp_region_name, de_opt_str);
        %% Process regional anisotropy data
        tmp_data_field_name = fieldnames(tmp_wb_data.rda_data.volume_weighted);
        if all(isfield(tmp_region_data, {'cc_1', 'cc_2'})) || isfield(tmp_region_data, 'cc_merged')
            for iter_cc = 1 : tmp_region_data.num_cc
                switch  tmp_region_data.num_cc
                    case 2
                        tmp_cc_data = tmp_region_data.(sprintf('cc_%d', iter_cc));
                    case 1
                        tmp_cc_data = tmp_region_data.cc_merged;
                end
                tmp_cc_data.local_cube_stat.dt_ori_pca = fun_structure_field_indexing(tmp_cc_data.local_cube_stat.dt_ori_pca, ...
                    tmp_cc_data.is_in_cc_cube_Q);
                
                tmp_cc_data.local_cube_stat.node_density_mm3 = tmp_cc_data.local_cube_stat.node_stat.num_data.degree ./ (0.24 * linear_scaling_factor)^3;
                
                tmp_cc_data.vw_rda = fun_structure_field_indexing(tmp_wb_data.rda_data.volume_weighted, ...
                    tmp_cc_data.is_in_cc_cube_Q);
                tmp_cc_data.lw_rda = fun_structure_field_indexing(tmp_wb_data.rda_data.length_weighted, ...
                    tmp_cc_data.is_in_cc_cube_Q);
                tmp_cc_data.rda_r_max = tmp_wb_data.rda_data.radius_max;
                tmp_cc_data.rda_r_min = tmp_wb_data.rda_data.radius_min;
                
                tmp_cc_data.dataset_name = tmp_region_data.dataset_name;
                tmp_cc_data.stack = tmp_region_data.stack;
                region_data{iter_cc, iter_stack, iter_region} = tmp_cc_data;
            end
        end
        toc(tmp_tic);
    end
end
fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));
%% Vessel anisotropy - collect data
min_dist_to_brain_surface_um = 100;
min_in_cc_vf = 0.5;
min_in_brain_mask_vf = 0.95;
min_depth_pen_big_vessel = 200;
max_depth_pen_big_vessel = inf;
rda_cap_idx = 10;
rda_noncap_idx = 17;

rda_cap_r_max = wb_data_cell{1}.rda_data.radius_max(rda_cap_idx);

fig_group_name = sprintf('capillary_angle_N_anisotropy_rcapmax_%.2f', rda_cap_r_max);
anisotropy_data = struct;
[anisotropy_data.mean_dist_to_surface, ...
    anisotropy_data.cos_cap_to_pen_vw, anisotropy_data.cos_cap_to_dt_ori_vw, ...
    anisotropy_data.an_cap_vw_fa, anisotropy_data.an_cap_vw_svd1, ...
    anisotropy_data.an_cap_vw_fa_p, anisotropy_data.an_cap_vw_fa_z] = deal(cell(2, num_stack, num_region));
for iter_region = 1 : num_region
    for iter_stack = 1 : num_stack
        for iter_cc = 1 : num_cc
            tmp_cc_data = region_data{iter_cc, iter_stack, iter_region};
            
            tmp_vw_cap_ori = tmp_cc_data.vw_rda.svd_max_vec(:, :, rda_cap_idx);
            tmp_vw_noncap_ori = tmp_cc_data.vw_rda.svd_max_vec(:, :, rda_noncap_idx);
            tmp_dt_ori = tmp_cc_data.local_cube_stat.dt_ori_pca.ori_eig_vec;
            % Compute the orientation of large vessels - volume weighted - remove
            % surface vessels
            tmp_cube_noncap_cap_cos = abs(sum(tmp_vw_cap_ori .* tmp_vw_noncap_ori, 2));
            tmp_cube_dt_ori_cap_cos = abs(sum(tmp_vw_cap_ori .* tmp_dt_ori, 2));
            %% Get mean of local 240-cube statistics
            tmp_gt_half_in_mask_Q = tmp_cc_data.local_cube_stat.cube_in_brain_mask_ratio > min_in_brain_mask_vf & ...
                tmp_cc_data.local_cube_stat.dt_ori_pca.avg_dist_2_surface >= min_dist_to_brain_surface_um;
            
            anisotropy_data.mean_dist_to_surface{iter_cc, iter_stack, iter_region} = tmp_cc_data.local_cube_stat.dt_ori_pca.avg_dist_2_surface(tmp_gt_half_in_mask_Q);
            anisotropy_data.cos_cap_to_pen_vw{iter_cc, iter_stack, iter_region} = tmp_cube_noncap_cap_cos(tmp_gt_half_in_mask_Q);
            anisotropy_data.cos_cap_to_dt_ori_vw{iter_cc, iter_stack, iter_region} = tmp_cube_dt_ori_cap_cos(tmp_gt_half_in_mask_Q);
            
            anisotropy_data.an_cap_vw_fa{iter_cc, iter_stack, iter_region} = tmp_cc_data.vw_rda.fa(tmp_gt_half_in_mask_Q, rda_cap_idx);
            
            anisotropy_data.an_cap_vw_svd1{iter_cc, iter_stack, iter_region} = tmp_cc_data.vw_rda.svd_1(tmp_gt_half_in_mask_Q, rda_cap_idx);
            anisotropy_data.an_cap_vw_fa_p{iter_cc, iter_stack, iter_region} = tmp_cc_data.vw_rda.fa_p(tmp_gt_half_in_mask_Q, rda_cap_idx);
            anisotropy_data.an_cap_vw_fa_z{iter_cc, iter_stack, iter_region} = tmp_cc_data.vw_rda.fa_z(tmp_gt_half_in_mask_Q, rda_cap_idx);
        end
    end
end
%%     Plot - Capillary anisitropy - merge all the cc to compute the regional curve - selected regions 
tmp_vis_region_list_ind = 1 : numel(analysis_region_id_list);
fig_y_left_axis_label = '|Cos(\theta_{cn})|';
fig_hdl = figure('Visible', 'on');
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 3;
ax_hdl_1 = subplot(2,2,1);
ax_hdl_2 = subplot(2,2,2);
ax_hdl_3 = subplot(2,2,3);
ax_hdl_4 = subplot(2,2,4);
hold(ax_hdl_1, 'on');
hold(ax_hdl_2, 'on');
hold(ax_hdl_3, 'on');
hold(ax_hdl_4, 'on');
for iter_region = 1 : num_region
    iter_stack = 1 : num_stack;
    tmp_dist = cat(1, anisotropy_data.mean_dist_to_surface{:, iter_stack, iter_region});
%     tmp_cos = cat(1, anisotropy_data.cos_cap_to_pen_vw{:, iter_stack, iter_region});
    tmp_cos = cat(1, anisotropy_data.cos_cap_to_dt_ori_vw{:, iter_stack, iter_region});
    
    tmp_fa = cat(1, anisotropy_data.an_cap_vw_fa{:, iter_stack, iter_region});
    tmp_svd1 = cat(1, anisotropy_data.an_cap_vw_svd1{:, iter_stack, iter_region});
    tmp_fa_p = cat(1, anisotropy_data.an_cap_vw_fa_p{:, iter_stack, iter_region});
    tmp_fa_p(tmp_fa_p == 0) = 5e-5;
    tmp_fa_z = cat(1, anisotropy_data.an_cap_vw_fa_z{:, iter_stack, iter_region});
    
    [tmp_ind, tmp_x] = fun_bin_data_to_idx_list_by_edges(tmp_dist, ...
        linspace(min(tmp_dist, [], 'all', 'omitnan'), max(tmp_dist, [], 'all', 'omitnan'), 20));
    tmp_selected_Q = tmp_x <= 1e3;
    tmp_ind = tmp_ind(tmp_selected_Q);
    tmp_x = tmp_x(tmp_selected_Q);
    [tmp_cos_binned] = fun_bin_data_to_cells_by_ind(abs(tmp_cos), tmp_ind);
    
    tmp_svd1_binned = fun_bin_data_to_cells_by_ind(tmp_svd1, tmp_ind);
    tmp_fa_binned = fun_bin_data_to_cells_by_ind(tmp_fa, tmp_ind);
    tmp_fa_p_binned = fun_bin_data_to_cells_by_ind(tmp_fa_p, tmp_ind);
    tmp_fa_z_binned = fun_bin_data_to_cells_by_ind(tmp_fa_z, tmp_ind);
    
    tmp_avg_cos_binned = cellfun(@(x) mean(x, 'omitnan'), tmp_cos_binned);
    tmp_avg_fa_binned = cellfun(@(x) mean(x, 'omitnan'), tmp_fa_binned);
    tmp_avg_fa_p_binned = cellfun(@(x) median(x, 'omitnan'), tmp_fa_p_binned);
    tmp_avg_svd1_binned = cellfun(@(x) mean(x, 'omitnan'), tmp_svd1_binned);
    tmp_avg_faz_binned = cellfun(@(x) mean(x, 'omitnan'), tmp_fa_z_binned);
	%% plot
    yyaxis(ax_hdl_1, 'left');
    plot(ax_hdl_1, tmp_x, tmp_avg_cos_binned, 'LineWidth', 2);
    ax_hdl_1.YLabel.String = fig_y_left_axis_label;
    ax_hdl_1.YLim = [0, 1];
    yyaxis(ax_hdl_1, 'right');
    plot(ax_hdl_1, tmp_x, tmp_avg_fa_binned, 'LineWidth', 2);
    ax_hdl_1.YLabel.String = 'Fractional anisotropy';
    ax_hdl_1.YLim = [0,1];
    
    yyaxis(ax_hdl_2, 'left');
    plot(ax_hdl_2, tmp_x, tmp_avg_cos_binned, 'LineWidth', 2);
    ax_hdl_2.YLabel.String = fig_y_left_axis_label;
    ax_hdl_2.YLim = [0, 1];
    yyaxis(ax_hdl_2, 'right');
    plot(ax_hdl_2, tmp_x, tmp_avg_fa_p_binned, 'LineWidth', 2);
    ax_hdl_2.YLabel.String = 'FA p-Value';
    ax_hdl_2.YLim = [1e-5,1];
    ax_hdl_2.YScale = 'log';

    yyaxis(ax_hdl_3, 'left');
    plot(ax_hdl_3, tmp_x, tmp_avg_cos_binned, 'LineWidth', 2);
    ax_hdl_3.YLabel.String = fig_y_left_axis_label;
    ax_hdl_3.YLim = [0, 1];
    yyaxis(ax_hdl_3, 'right');
    plot(ax_hdl_3, tmp_x, tmp_avg_svd1_binned, 'LineWidth', 2);
    ax_hdl_3.YLabel.String = 'PCV1';
    ax_hdl_3.YLim = [0,1];

    yyaxis(ax_hdl_4, 'left');
    plot(ax_hdl_4, tmp_x, tmp_avg_cos_binned, 'LineWidth', 2);
    ax_hdl_4.YLabel.String = fig_y_left_axis_label;
    ax_hdl_4.YLim = [0, 1];
    yyaxis(ax_hdl_4, 'right');
    plot(ax_hdl_4, tmp_x, tmp_avg_faz_binned, 'LineWidth', 2);
    ax_hdl_4.YLabel.String = 'FA_z';
    ax_hdl_4.YLim = [0, 6];
end
ax_hdl_1.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_1.XLim = [0, 1100];
ax_hdl_2.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_2.XLim = [0, 1100];
ax_hdl_3.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_3.XLim = [0, 1100];
ax_hdl_4.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_4.XLim = [0, 1100];

legend(ax_hdl_1, analysis_region_name{tmp_vis_region_list_ind});
legend(ax_hdl_2, analysis_region_name{tmp_vis_region_list_ind});
legend(ax_hdl_3, analysis_region_name{tmp_vis_region_list_ind});
legend(ax_hdl_4, analysis_region_name{tmp_vis_region_list_ind});
fig_name = fullfile(visualization_folder, fig_group_name, sprintf('%s_%s_%s_vs_depth_in_%s.png', dataset_name, merged_stack_name, fig_group_name, parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_name);
%%     Plot - Capillary anisotropy and pValue 2 figures
tmp_vis_region_list_ind = 1 : numel(analysis_region_id_list);
fig_hdl = figure('Visible', 'on');
ax_hdl_1 = subplot(2,1,1);
ax_hdl_2 = subplot(2,1,2);
hold(ax_hdl_1, 'on');
hold(ax_hdl_2, 'on');
for iter_region = 1 : num_region
    iter_stack = 1 : num_stack;
    tmp_dist = cat(1, anisotropy_data.mean_dist_to_surface{:, iter_stack, iter_region});
    tmp_cos = cat(1, anisotropy_data.cos_cap_to_dt_ori_vw{:, iter_stack, iter_region});
    
    tmp_fa_p = cat(1, anisotropy_data.an_cap_vw_fa_p{:, iter_stack, iter_region});
    tmp_fa_p(tmp_fa_p == 0) = 5e-5;
    
    [tmp_ind, tmp_x] = fun_bin_data_to_idx_list_by_edges(tmp_dist, ...
        linspace(min(tmp_dist, [], 'all', 'omitnan'), max(tmp_dist, [], 'all', 'omitnan'), 20));
    tmp_selected_Q = tmp_x <= 1e3;
    tmp_ind = tmp_ind(tmp_selected_Q);
    tmp_x = tmp_x(tmp_selected_Q);
    tmp_cos_binned = fun_bin_data_to_cells_by_ind(abs(tmp_cos), tmp_ind);
    tmp_fa_p_binned = fun_bin_data_to_cells_by_ind(tmp_fa_p, tmp_ind);
    
    tmp_avg_cos_binned = cellfun(@(x) mean(x, 'omitnan'), tmp_cos_binned);
    tmp_avg_fa_p_binned = cellfun(@(x) median(x, 'omitnan'), tmp_fa_p_binned);
%% Plot
    plt_hdl = plot(ax_hdl_1, tmp_x, tmp_avg_cos_binned, 'LineWidth', 1.5);
    plot(ax_hdl_2, tmp_x, tmp_avg_fa_p_binned, 'LineWidth', 1.5, 'Color', plt_hdl.Color);  
end
ax_hdl_1.YLabel.String = 'Average |Cos(\theta_{cn})|';
ax_hdl_1.YLim = [0, 1];
grid(ax_hdl_1, 'on');
box(ax_hdl_1, 'on');

ax_hdl_2.YLabel.String = 'Median FA_p';
ax_hdl_2.YLim = [1e-5,1];
ax_hdl_2.YScale = 'log';
ax_hdl_2.YTick = 10 .^ (-5 : 1 : 1);
grid(ax_hdl_2, 'on');
box(ax_hdl_2, 'on');

[ax_hdl_1.XLabel.String, ax_hdl_2.XLabel.String] = deal('Distance to the cortical surface (\mum)');
[ax_hdl_1.XLim, ax_hdl_2.XLim] = deal([0, 1100]);
legend(ax_hdl_1, analysis_region_name{tmp_vis_region_list_ind});
fig_name = fullfile(visualization_folder, fig_group_name, sprintf('%s_%s_%s_vs_depth_in_%s.png', dataset_name, merged_stack_name, 'Cos_and_FAp', parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_name);
%% Average over cc and stack curve
%% 	Extract data
anisotropy_vs_depth_cell = cell(num_region, 1);
ai_extract_field_name = fieldnames(anisotropy_data);
ai_feature_field_name = setdiff(ai_extract_field_name, {'mean_dist_to_surface', 'an_cap_vw_fa_p'}, 'sorted');
ai_feature_ylabel = {'Fractional anisotropy', 'FA_z', 'PCV_1', '|Cos(\theta_{cn})|', '|Cos(\theta_{cp})|'};
for iter_region = 1 : num_region
    % Bin the data by depth and compute the average value for each depth
    tmp_data_str = struct;
    for iter_field = 1 : numel(ai_extract_field_name)
        tmp_field_name = ai_extract_field_name{iter_field};
        tmp_data_str.(tmp_field_name) = anisotropy_data.(tmp_field_name)(:, :, iter_region);
        [tmp_data_str.legend, tmp_data_str.(sprintf('avg_%s', tmp_field_name))] = deal(cell(size(tmp_data_str.(tmp_field_name))));
    end
    tmp_data_str.data_size = numel(tmp_data_str.mean_dist_to_surface);
    for iter_ind = 1 : tmp_data_str.data_size
        [iter_cc, iter_stack] = ind2sub([2, num_stack], iter_ind);
        tmp_data_str.legend{iter_cc, iter_stack} = sprintf('%s %s', strrep(stack_list{iter_stack}, '_', ''), cc_label{iter_cc});
        
        tmp_dist_data = tmp_data_str.mean_dist_to_surface{iter_ind};        
        tmp_dist_bin = linspace(min(tmp_dist_data, [], 'all', 'omitnan'), max(tmp_dist_data, [], 'all', 'omitnan'), ...
            min(20, ceil(numel(tmp_dist_data) / 10)));
       [tmp_list_ind, tmp_dist] = fun_bin_data_to_idx_list_by_edges(tmp_dist_data, tmp_dist_bin);
       for iter_field = 1 : numel(ai_extract_field_name)
           tmp_field_name = ai_extract_field_name{iter_field};
           tmp_data_to_bin = tmp_data_str.(tmp_field_name){iter_ind};
           tmp_data_binned = fun_bin_data_to_cells_by_ind(tmp_data_to_bin, tmp_list_ind);
           tmp_data_str.(sprintf('avg_%s', tmp_field_name)){iter_ind} = cellfun(@(x) mean(x, 'omitnan'), tmp_data_binned);
       end
       tmp_field_name = 'an_cap_vw_fa_p';
       tmp_data_to_bin = tmp_data_str.(tmp_field_name){iter_ind};
       tmp_data_to_bin(tmp_data_to_bin == 0) = 5e-5;
       tmp_data_binned = fun_bin_data_to_cells_by_ind(tmp_data_to_bin, tmp_list_ind);
       tmp_data_str.(sprintf('avg_%s', tmp_field_name)){iter_ind} = cellfun(@(x) median(x, 'omitnan'), tmp_data_binned);
    end
    
    for iter_field = 1 : numel(ai_feature_field_name)
        tmp_field_name = ai_feature_field_name{iter_field};
        tmp_data_str.(sprintf('avg_itp_%s', tmp_field_name)) = fun_analysis_get_xy_curve_avgNstd_by_interpolation(...
            tmp_data_str.avg_mean_dist_to_surface, tmp_data_str.(sprintf('avg_%s', tmp_field_name)));
    end
    
    tmp_field_name = 'an_cap_vw_fa_p';
    tmp_data_str.(sprintf('avg_itp_%s', tmp_field_name)) = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
            tmp_data_str.avg_mean_dist_to_surface, tmp_data_str.(sprintf('avg_%s', tmp_field_name)));
    
    anisotropy_vs_depth_cell{iter_region} = tmp_data_str;
end
%%  Plot - double y-axis, single region
for iter_region = 1 : num_region
    tmp_data_str = anisotropy_vs_depth_cell{iter_region};
    %% Plot
    fig_hdl = figure('Visible', 'on');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2,3] .* 1.2;
    num_plot = numel(ai_feature_field_name);
    
    for iter_ax = 1 : num_plot
        tmp_ax = subplot(3,2,iter_ax);
        for iter_cell = 1 : tmp_data_str.data_size
            plot(tmp_ax, tmp_data_str.avg_mean_dist_to_surface{iter_cell}, tmp_data_str.(sprintf('avg_%s', ai_feature_field_name{iter_ax})){iter_cell}, 'LineWidth', 2, 'LineStyle', '--');
            hold(tmp_ax, 'on');
        end        
        errorbar(tmp_ax, tmp_data_str.(sprintf('avg_itp_%s', ai_feature_field_name{iter_ax})).interpolate_x, ...
            tmp_data_str.(sprintf('avg_itp_%s', ai_feature_field_name{iter_ax})).y_avg, ...
            tmp_data_str.(sprintf('avg_itp_%s', ai_feature_field_name{iter_ax})).y_std, 'LineWidth', 3);
        tmp_ax.YLabel.String = ai_feature_ylabel{iter_ax};
        tmp_ax.XLabel.String = 'Distance to the cortical surface (\mum)';
        tmp_ax.Title.String = sprintf('%s', analysis_region_name{iter_region});
        
        tmp_legend = cat(1, tmp_data_str.legend(:), {'Average \pm STD'});
        leg_hdl = legend(tmp_ax, tmp_legend{:}, 'Location', 'best');
        if ~strcmp(ai_feature_ylabel{iter_ax}, 'FA_z')
            tmp_ax.YLim = [0,1];
        else
            tmp_ax.YLim(1) = 0;
        end
    end
    tmp_ax_pv = subplot(3,2,num_plot + 1);
    for iter_cell = 1 : tmp_data_str.data_size
        plot(tmp_ax_pv, tmp_data_str.avg_mean_dist_to_surface{iter_cell}, tmp_data_str.avg_an_cap_vw_fa_p{iter_cell}, 'LineWidth', 2, 'LineStyle', '--');
        hold(tmp_ax_pv, 'on');
    end
    plot(tmp_ax_pv, tmp_data_str.avg_itp_an_cap_vw_fa_p.interpolate_x,...
        tmp_data_str.avg_itp_an_cap_vw_fa_p.prctile_value(:, 4), 'LineWidth', 3);
    tmp_ax_pv.YScale = 'log';
    tmp_legend = cat(1, tmp_data_str.legend(:), {'Median'});
    leg_hdl = legend(tmp_ax_pv, tmp_legend{:}, 'Location', 'best');
    tmp_ax_pv.XLabel.String = 'Distance to the cortical surface (\mum)';
    tmp_ax_pv.YLabel.String = 'Median FA p-Value';
    tmp_ax_pv.YLim = [1e-5, 1];
    tmp_ax_pv.YTick = [1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1];
    tmp_ax_pv.Title.String = sprintf('%s', analysis_region_name{iter_region});
    
    fig_name = fullfile(visualization_folder, fig_group_name, sprintf('%s_%s_%s_vs_depth_in_%s_cc_avg.png', ...
        dataset_name, merged_stack_name, fig_group_name, strrep(analysis_region_name{iter_region}, ' ' , '_')));
    fun_print_image_in_several_formats(fig_hdl, fig_name);
    delete(fig_hdl);
end
%%  Plot - double y-axis, Mean for selected regions 
fig_hdl = figure('Visible', 'on');
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 3;
ax_hdl_1 = subplot(2,2,1);
ax_hdl_2 = subplot(2,2,2);
ax_hdl_3 = subplot(2,2,3);
ax_hdl_4 = subplot(2,2,4);
hold(ax_hdl_1, 'on');
hold(ax_hdl_2, 'on');
hold(ax_hdl_3, 'on');
hold(ax_hdl_4, 'on');
for iter_region = 1 : num_region
    tmp_data_str = anisotropy_vs_depth_cell{iter_region};
    
    tmp_x = tmp_data_str.avg_itp_an_cap_vw_fa.interpolate_x;
    tmp_dist_selected_Q = tmp_x < 1e3;
    tmp_x = tmp_x(tmp_dist_selected_Q);
    tmp_cos = tmp_data_str.avg_itp_cos_cap_to_dt_ori_vw.y_avg(tmp_dist_selected_Q);
	% plot
    yyaxis(ax_hdl_1, 'left');
    plot(ax_hdl_1, tmp_x, tmp_cos, 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_1.YLabel.String = fig_y_left_axis_label;
    ax_hdl_1.YLim = [0, 1];
    yyaxis(ax_hdl_1, 'right');
    plot(ax_hdl_1, tmp_x, tmp_data_str.avg_itp_an_cap_vw_fa.y_avg(tmp_dist_selected_Q), 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_1.YLabel.String = 'Fractional anisotropy';
    ax_hdl_1.YLim = [0,1];
    
    yyaxis(ax_hdl_2, 'left');
    plot(ax_hdl_2, tmp_x, tmp_cos, 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_2.YLabel.String = fig_y_left_axis_label;
    ax_hdl_2.YLim = [0, 1];
    yyaxis(ax_hdl_2, 'right');
    plot(ax_hdl_2, tmp_x, tmp_data_str.avg_itp_an_cap_vw_fa_p.prctile_value(tmp_dist_selected_Q, 4), 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_2.YLabel.String = 'Median FA p-Value';
    ax_hdl_2.YLim = [1e-5,1];
    ax_hdl_2.YScale = 'log';

    yyaxis(ax_hdl_3, 'left');
    plot(ax_hdl_3, tmp_x, tmp_cos, 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_3.YLabel.String = fig_y_left_axis_label;
    ax_hdl_3.YLim = [0, 1];
    yyaxis(ax_hdl_3, 'right');
    plot(ax_hdl_3, tmp_x, tmp_data_str.avg_itp_an_cap_vw_svd1.y_avg(tmp_dist_selected_Q), 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_3.YLabel.String = 'PCV1';
    ax_hdl_3.YLim = [0,1];

    yyaxis(ax_hdl_4, 'left');
    plot(ax_hdl_4, tmp_x, tmp_cos, 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_4.YLabel.String = fig_y_left_axis_label;
    ax_hdl_4.YLim = [0, 1];
    yyaxis(ax_hdl_4, 'right');
    plot(ax_hdl_4, tmp_x, tmp_data_str.avg_itp_an_cap_vw_fa_z.y_avg(tmp_dist_selected_Q), 'LineWidth', 2, 'MarkerSize', 3);
    ax_hdl_4.YLabel.String = 'FA_z';
    ax_hdl_4.YLim = [0, 6];

end
ax_hdl_1.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_1.XLim = [0, 1100];
ax_hdl_2.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_2.XLim = [0, 1100];
ax_hdl_3.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_3.XLim = [0, 1100];
ax_hdl_4.XLabel.String = 'Distance to the cortical surface(\mum)';
ax_hdl_4.XLim = [0, 1100];

legend(ax_hdl_1, analysis_region_name{tmp_vis_region_list_ind});
legend(ax_hdl_2, analysis_region_name{tmp_vis_region_list_ind});
legend(ax_hdl_3, analysis_region_name{tmp_vis_region_list_ind});
legend(ax_hdl_4, analysis_region_name{tmp_vis_region_list_ind});

fig_name = fullfile(visualization_folder, fig_group_name, sprintf('%s_%s_%s_vs_depth_cc_avg_in_%s.png', dataset_name, merged_stack_name, fig_group_name, parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_name);
%%  Plot angle and FAp two subfigures
fig_hdl = figure('Visible', 'on');
ax_hdl_1 = subplot(2,1,1);
ax_hdl_2 = subplot(2,1,2);
hold(ax_hdl_1, 'on');
hold(ax_hdl_2, 'on');
for iter_region = 1 : num_region
    tmp_data_str = anisotropy_vs_depth_cell{iter_region};
    
    tmp_x = tmp_data_str.avg_itp_an_cap_vw_fa.interpolate_x;
    tmp_dist_selected_Q = tmp_x < 1e3;
    tmp_x = tmp_x(tmp_dist_selected_Q);
    tmp_cos = tmp_data_str.avg_itp_cos_cap_to_dt_ori_vw.y_avg(tmp_dist_selected_Q);
    tmp_cos_sd = tmp_data_str.avg_itp_cos_cap_to_dt_ori_vw.y_std(tmp_dist_selected_Q);
    tmp_fap = tmp_data_str.avg_itp_an_cap_vw_fa_p.prctile_value(tmp_dist_selected_Q, 4);
%% Plot
    plt_hdl = plot(ax_hdl_1, tmp_x, tmp_cos, 'LineWidth', 1.5);
%     [ax_hdl_1, plt_hdl, patch_hdl] = fun_vis_confidence_interval_shaded(tmp_x, tmp_cos, tmp_cos - tmp_cos_sd, ...
%         tmp_cos + tmp_cos_sd, ax_hdl_1);
    
    plot(ax_hdl_2, tmp_x, tmp_fap, 'LineWidth', 1.5, 'Color', plt_hdl.Color);  
end
ax_hdl_1.YLabel.String = 'Average |Cos(\theta_{cn})|';
ax_hdl_1.YLim = [0, 1];
grid(ax_hdl_1, 'on');
box(ax_hdl_1, 'on');

ax_hdl_2.YLabel.String = 'Median FA_p';
ax_hdl_2.YLim = [1e-5,1];
ax_hdl_2.YScale = 'log';
ax_hdl_2.YTick = 10 .^ (-5 : 1 : 1);
grid(ax_hdl_2, 'on');
box(ax_hdl_2, 'on');

[ax_hdl_1.XLabel.String, ax_hdl_2.XLabel.String] = deal('Distance to the cortical surface (\mum)');
[ax_hdl_1.XLim, ax_hdl_2.XLim] = deal([0, 1100]);
legend(ax_hdl_1, analysis_region_name{tmp_vis_region_list_ind});
fig_name = fullfile(visualization_folder, fig_group_name, sprintf('%s_%s_%s_vs_depth_cc_avg_in_%s.png', dataset_name, merged_stack_name, 'Cos_and_FAp', parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_name);
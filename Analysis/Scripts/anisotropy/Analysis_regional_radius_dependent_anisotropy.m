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
end
%%
parent_str_name = 'Major_brain_regions';
analysis_region_id_list = {453 961 1080 477 549 [302 294] 4 771 354 528 776, 1097}; % Remove hypothalamus, use hippocampus region
vis_subfolder_name = 'Region_comparsion';

tmp_vis_region_list_ind = 1 : numel(analysis_region_id_list);
num_region = numel(analysis_region_id_list);
analysis_region_name = cell(num_region, 1);
for iter_region = 1 : num_region
    tmp_atlas_id = analysis_region_id_list{iter_region};
    if isscalar(tmp_atlas_id)
        analysis_region_name{iter_region} = allen_atlas.structure_table.name{full(allen_atlas.id_2_ind(tmp_atlas_id))};
    end
    tmp_atlas_id = full(allen_atlas_map_old_id_to_new_id(tmp_atlas_id));
    analysis_region_id_list{iter_region} = tmp_atlas_id;
end
analysis_region_name{6} = 'Superior colliculus';
%% Collect all the region statistics in all the stack
region_data = cell(num_cc, num_stack, num_region);
collect_data_tic = tic;
for iter_region = 1 : num_region
    for iter_stack = 1 : num_stack
        tmp_tic = tic;
        tmp_region_id = analysis_region_id_list{iter_region};
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
                tmp_cc_data.local_cube_stat.node_density_mm3 = tmp_cc_data.local_cube_stat.node_stat.num_data.degree ./ (0.24 * linear_scaling_factor)^3;
                
                tmp_cc_data.vw_rda = fun_structure_field_indexing(tmp_wb_data.rda_data.volume_weighted, ...
                    tmp_cc_data.is_in_cc_cube_Q);
                tmp_cc_data.lw_rda = fun_structure_field_indexing(tmp_wb_data.rda_data.length_weighted, ...
                    tmp_cc_data.is_in_cc_cube_Q);
                tmp_cc_data.rda_r_max = tmp_wb_data.rda_data.radius_max;
                tmp_cc_data.rda_r_min = tmp_wb_data.rda_data.radius_min;
                region_data{iter_cc, iter_stack, iter_region} = tmp_cc_data;
            end
        end
        toc(tmp_tic);
    end
end
fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));
%% Merge all the data inside a region
% Radius dependence anisotropy in all the region
min_in_cc_vf = 0.75;
min_in_brain_mask_vf = 1;
if de_opt_str.merge_cc_Q
    merged_cc_data = cell(num_region, 1);
    rda_field_name = {'vw_rda', 'lw_rda'};
    for iter_region = 1 : num_region
        tmp_region_cell = region_data(:, :, iter_region);
        tmp_mc_str = struct;
        %% Extract radius-dependent anisotropy data
        for iter_rda_field = 1 : numel(rda_field_name)
            tmp_field_name = rda_field_name{iter_rda_field};
            tmp_extract_subfield_names = fieldnames(tmp_region_cell{1}.(tmp_field_name));
            tmp_mc_str.r_max = tmp_region_cell{1}.rda_r_max;
            tmp_mc_str.r_min = tmp_region_cell{1}.rda_r_min;
            % Selection
            for iter_field = 1 : numel(tmp_extract_subfield_names)
                tmp_fn = tmp_extract_subfield_names{iter_field};
                tmp_in_cc_vol_f = cellfun(@(x) x.local_cube_stat.in_cc_vol_f, tmp_region_cell, 'UniformOutput', false);
                tmp_in_cc_vol_f = cat(1, tmp_in_cc_vol_f{:});
                
                tmp_in_brain_mask_vol_f = cellfun(@(x) x.local_cube_stat.cube_in_brain_mask_ratio, tmp_region_cell, 'UniformOutput', false);
                tmp_in_brain_mask_vol_f = cat(1, tmp_in_brain_mask_vol_f{:});
                
                tmp_selected_Q = (tmp_in_cc_vol_f >= min_in_cc_vf) & ...
                    (tmp_in_brain_mask_vol_f >= min_in_brain_mask_vf);
                tmp_cell = cellfun(@(x) x.(tmp_field_name).(tmp_fn), tmp_region_cell, 'UniformOutput', false);
                tmp_cell = cat(1, tmp_cell{:});
                if ~isreal(tmp_cell)
                    warning('%s is not real array. Take the real part\n', tmp_fn);
                    tmp_cell = real(tmp_cell);
                end
                tmp_cell = tmp_cell(tmp_selected_Q, :);
                % Set the minimum value of p-value to be 5e-5
                if any(strcmp(tmp_fn, {'fa_p', 'svd_p'}))
                    tmp_is_zero_Q = (tmp_cell == 0);
                    tmp_cell(tmp_is_zero_Q) = 5e-5;
                end
                tmp_mc_str.(tmp_field_name).data.(tmp_fn) = tmp_cell;
                tmp_mc_str.(tmp_field_name).stat.(tmp_fn) = fun_analysis_get_basic_statistics_in_column(tmp_cell, 1);
            end
        end
        %% Extract cube statistics?
        %%
        merged_cc_data{iter_region} = tmp_mc_str;
    end
    %% RDA FA-median and FA_p median for selected regions
    selected_region_idx = [1, 4, 6, 7, 8, 11];
    vis_feature_list = {'fa', 'fa_p'};
    num_vis_features = numel(vis_feature_list);
    vis_feature_stat = 'median';
    vis_feature_disp_name_list = {'FA', 'FA_p'};
    vis_y_scale = {'linear', 'log', 'linear'};
    vis_feature_disp_name_list = cellfun(@(x) sprintf('%s %s', vis_feature_stat, x), ...
        vis_feature_disp_name_list, 'UniformOutput', false);
    vix_y_range = {[0, 1], [1e-5, 1], [0, 1]};
    
    tmp_x = merged_cc_data{1}.r_max;
    tmp_x_selected_idx = 1 : (numel(tmp_x) - 2);
    tmp_x = tmp_x(tmp_x_selected_idx);
    
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 1.3];
    for iter_feature = 1 : num_vis_features
        leg_target_hdl = [];
        tmp_feature_name = vis_feature_list{iter_feature};
        tmp_y = cellfun(@(x) x.vw_rda.stat.(tmp_feature_name).(vis_feature_stat)(tmp_x_selected_idx), merged_cc_data, 'UniformOutput', false);
        tmp_y_up = cellfun(@(x) cat(1, x.vw_rda.stat.(tmp_feature_name).('prctile_val'){:}), merged_cc_data, 'UniformOutput', false);
        tmp_y_low = cellfun(@(x) x(tmp_x_selected_idx, 7), tmp_y_up, 'UniformOutput', false);
        tmp_y_up = cellfun(@(x) x(tmp_x_selected_idx, 9), tmp_y_up, 'UniformOutput', false);
        tmp_ax = subplot(num_vis_features, 1, iter_feature);
        for iter_plt = selected_region_idx
            tmp_plt_hdl = plot(tmp_ax, tmp_x, tmp_y{iter_plt}, 'LineWidth', 2);
            hold(tmp_ax, 'on');
            %         tmp_hdl = plot(tmp_ax, tmp_x, tmp_y_low{iter_plt}, 'LineWidth', 1.5, 'LineStyle', '-.', 'Color', tmp_plt_hdl.Color);
            %         tmp_hdl = plot(tmp_ax, tmp_x, tmp_y_up{iter_plt}, 'LineWidth', 1.5, 'LineStyle', '--', 'Color', tmp_plt_hdl.Color);
            %         [tmp_ax, tmp_plt_hdl, tmp_patch_hdl] = fun_vis_confidence_interval_shaded(tmp_x, tmp_y{iter_plt}, ...
            %             tmp_y_low{iter_plt}, tmp_y_up{iter_plt}, tmp_ax);
            leg_target_hdl(end+1) = tmp_plt_hdl;
        end
        tmp_ax.XScale = 'log';
        tmp_ax.XTick = 2 : 10;
        grid(tmp_ax, 'on');
        tmp_ax.XLabel.String = 'Maximum radius (\mum)';
        tmp_ax.YLabel.String = vis_feature_disp_name_list{iter_feature};
        tmp_ax.YScale = vis_y_scale{iter_feature};
        tmp_ax.YLim = vix_y_range{iter_feature};
        if iter_feature == 1
            leg_hdl = legend(leg_target_hdl, analysis_region_name{selected_region_idx}, 'Location', 'northwest');
        end
    end
    fig_fp = fullfile(vis_folder_fp, sprintf('%s_%s_merge_cc_RDA_%s_FA_FAp_selected_region_cc.png', dataset_name, merged_stack_name, vis_feature_stat));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
else
    %% Average over cc
    avg_cc_data = cell(num_region, 1);
    rda_field_name = {'vw_rda', 'lw_rda'};
    for iter_region = 1 : num_region
        tmp_region_cell = region_data(:, :, iter_region);
        tmp_avg_str = struct;
        %% Extract radius-dependent anisotropy data
        for iter_rda_field = 1 : numel(rda_field_name)
            tmp_field_name = rda_field_name{iter_rda_field};
            tmp_extract_subfield_names = fieldnames(tmp_region_cell{1}.(tmp_field_name));
            tmp_avg_str.r_max = tmp_region_cell{1}.rda_r_max;
            tmp_avg_str.r_min = tmp_region_cell{1}.rda_r_min;
            for iter_field = 1 : numel(tmp_extract_subfield_names)
                tmp_fn = tmp_extract_subfield_names{iter_field};
                tmp_cell = cellfun(@(x) x.(tmp_field_name).(tmp_fn), tmp_region_cell, 'UniformOutput', false);
                
                tmp_in_brain_mask_vol_f = cellfun(@(x) x.local_cube_stat.cube_in_brain_mask_ratio, tmp_region_cell, 'UniformOutput', false);
                tmp_in_cc_vol_f = cellfun(@(x) x.local_cube_stat.in_cc_vol_f, tmp_region_cell, 'UniformOutput', false);
               
                tmp_avg_stat_cell = cell(size(tmp_cell));
                for iter_cell = 1 : numel(tmp_cell)
                    tmp_data = tmp_cell{iter_cell};
                    if ~isreal(tmp_data)
                        warning('%s is not real array. Take the real part\n', tmp_fn);
                        tmp_data = real(tmp_data);
                    end
                    tmp_selected_Q = (tmp_in_cc_vol_f{iter_cell} >= min_in_cc_vf) & ...
                        (tmp_in_brain_mask_vol_f{iter_cell} >= min_in_brain_mask_vf);
                    tmp_data = tmp_data(tmp_selected_Q, :);
                    if any(strcmp(tmp_fn, {'fa_p', 'svd_p'}))
                        tmp_is_zero_Q = (tmp_data == 0);
                        tmp_data(tmp_is_zero_Q) = 5e-5;
                    end                    
                    tmp_avg_stat_cell{iter_cell} = fun_analysis_get_basic_statistics_in_column(tmp_data, 1);
                end
                tmp_avg_stat_cell = cat(1, tmp_avg_stat_cell{:});
                tmp_avg_str.(tmp_field_name).(tmp_fn) = tmp_avg_stat_cell;
            end
        end
        %%
        avg_cc_data{iter_region} = tmp_avg_str;
    end
    %% RDA FA-median and FA_p median for selected regions
    selected_region_idx = [1, 4, 5, 6, 7, 8];
    vis_feature_list = {'fa', 'fa_p'};
    num_vis_features = numel(vis_feature_list);
    vis_feature_disp_name_list = {'Fractional anisotropy', 'FA_p'};
    vis_y_scale = {'linear', 'log', 'linear'};
    vis_feature_stat = 'median';
    vix_y_range = {[0, 1], [1e-5, 1], [0, 1]};
    
    tmp_x = avg_cc_data{1}.r_max;
    tmp_x_selected_idx = 1 : (numel(tmp_x) - 2);
    tmp_x = tmp_x(tmp_x_selected_idx);
    
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 1.5];
    for iter_feature = 1 : num_vis_features
        leg_target_hdl = [];
        tmp_feature_name = vis_feature_list{iter_feature};
        % Organize data
        tmp_y_cell = cellfun(@(x) {x.vw_rda.(tmp_feature_name).(vis_feature_stat)}, avg_cc_data, 'UniformOutput', false);
        tmp_x_cell = arrayfun(@(x) tmp_x, ones(num_cc * num_stack, 1), 'UniformOutput', false);
        % Plot
        tmp_ax = subplot(num_vis_features, 1, iter_feature);
        for iter_plt = selected_region_idx
            tmp_y = tmp_y_cell{iter_plt};
            tmp_y = cellfun(@(x) x(tmp_x_selected_idx), tmp_y, 'UniformOutput', false);
            tmp_avg_str = fun_analysis_get_xy_curve_stat_curves_by_interpolation(...
                tmp_x_cell, tmp_y, tmp_x);
            tmp_plt_hdl = errorbar(tmp_ax, tmp_avg_str.interpolate_x, ...
                tmp_avg_str.y_avg, tmp_avg_str.y_std, 'LineWidth', 1.5);
%             [tmp_ax, tmp_plt_hdl, tmp_patch_hdl] = fun_vis_errorbar_shaded(tmp_avg_str.interpolate_x, ...
%                 tmp_avg_str.y_avg, tmp_avg_str.y_std, tmp_ax);
            hold(tmp_ax, 'on');
            leg_target_hdl(end+1) = tmp_plt_hdl;
        end
        tmp_ax.XScale = 'log';
        tmp_ax.XTick = 1 : 10;
        grid(tmp_ax, 'on');
        tmp_ax.XLabel.String = 'Maximum radius (\mum)';
        tmp_ax.YLabel.String = vis_feature_disp_name_list{iter_feature};
        tmp_ax.YScale = vis_y_scale{iter_feature};
        tmp_ax.YLim = vix_y_range{iter_feature};
        if iter_feature == 1
            leg_hdl = legend(leg_target_hdl, analysis_region_name{selected_region_idx}, 'Location', 'southeast');
        end
    end   
end
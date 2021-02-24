clc;clear;
set_env
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
skl_grid_name = '240_cube_rec';
grid_c_version = '240_cube_combined_5_o_2';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
if ~isfield(grid_c, 'internal_subgrid_label_array')
    grid_c = fun_grid_get_internal_subgrid(grid_c);
    DataManager.write_grid_info(grid_c, grid_c.dataset_name, grid_c.stack, grid_c.version);
end
task_function_name = 'fun_task_simulation_pO2';
task_name = sprintf('%s_%s', 'Simulation_pO2_SA', datestr(now, 'YYYYmmDD'));

vis_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Tissue_pO2_simulation');
%% Basic information
task_str = struct;
task_str.dataset_name = dataset_name;
task_str.stack = stack;
task_str.DataManager = DataManager;
task_str.grid_c_name = grid_c_version;
task_str.grid_c_info = DataManager.load_grid(task_str.dataset_name, ...
    task_str.stack, task_str.grid_c_name);
task_str.task_function_name = task_function_name;

task_str.overwrite_Q = false;
%% Task specific parameters
    %% Simulation setting
    task_specific_opt = struct;
    task_specific_opt.recon_max_error_rate = 0.1; % Reconstruction 
    task_specific_opt.recon_pad_half_length = 80; % Before downsampling
    task_specific_opt.local_extrema_window_size = [8 : 4 : 32, 40 : 10 : 70, 80 : 20 : 120]; % Before downsampling 
    task_specific_opt.recon_downsample_rate = 2; % Downsample the reconstruction to reduce computation load
    task_specific_opt.min_recon_vessel_radius_um = 0; % The minimum radius of the vessel used for reconstruction. Not sure if 2 is a more appropriate number. 
    % Poisson equation 
    task_specific_opt.inhomogeneous_term = 1; % Also before downsampling
    % Krogh model
    task_specific_opt.krogh_vessel_r_um = 2;
    task_specific_opt.krogh_est_corr_coeff = 0.5;
    
    task_specific_opt.use_gpu_Q = true;
    task_specific_opt.save_folder_name = 'pO2_SA';
    
    task_str.skl_grid_name = skl_grid_name;

    task_str.task_name = task_name;
    task_str.opt = task_specific_opt;
%% 
grid_c_label_list = cell(2, 0);
grid_c_label_list(:, end+1) = {'Neocortex', [779 854 940 1014]};
grid_c_label_list(:, end+1) = {'Inferior colliculus', [2149 2137]};
grid_c_label_list(:, end+1) = {'Piriform cortex', [1042 935]};
% grid_c_label_list(:, end+1) = {'Corpus callsoum', [791 855]};
num_regions = size(grid_c_label_list, 2);
fig_legend = grid_c_label_list(1, :);
region_pO2_SA_data = cell(num_regions, 1);
for iter_region = 1 : num_regions
    region_grid_c_label = grid_c_label_list{2, iter_region};
    tmp_num_grid_c = numel(region_grid_c_label);
    tmp_region_pO2_data = cell(tmp_num_grid_c, 1);
    for iter_cube = 1 : tmp_num_grid_c
        tmp_grid_c_label = region_grid_c_label(iter_cube);
        [exit_code, int_cube_stat_cell] = fun_simulation_OT_internal_subgrid_SA(grid_c, ...
            skl_grid_name, tmp_grid_c_label, task_str.opt);
        tmp_region_pO2_data{iter_cube} = int_cube_stat_cell;
    end
    region_pO2_SA_data{iter_region} = tmp_region_pO2_data;
end
%% Post processing
lm_window_size_list = region_pO2_SA_data{1}{1}{1}.local_extrema_window_size;
diff_window_size = diff(lm_window_size_list);
tmp_x = movmean(lm_window_size_list, 2, 'Endpoints', 'discard');
leg_hdl = [];
fig_hdl = figure;
ax_hdl = subplot(2,1,1);
ax_hdl_2 = subplot(2,1,2);
for iter_region = 1 : 3
    tmp_region_data = region_pO2_SA_data{iter_region};
    tmp_region_data = cat(1, tmp_region_data{:});    
    for iter_cell = 1 : numel(tmp_region_data)
        tmp_cube_data = tmp_region_data{iter_cell};
        [tmp_cube_data.dt_lm_stat.dt_v_mean, tmp_cube_data.dt_lm_stat.dt_v_median]...
            = fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_data.dt_lm.dt_v);
        [tmp_cube_data.dt_lm_stat.pO2_v_mean, tmp_cube_data.dt_lm_stat.pO2_v_median]...
            = fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_data.dt_lm.pO2_v);
        
        [tmp_cube_data.pO2_lm_stat.dt_v_mean, tmp_cube_data.pO2_lm_stat.dt_v_median]...
            = fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_data.pO2_lm.dt_v);
        [tmp_cube_data.pO2_lm_stat.pO2_v_mean, tmp_cube_data.pO2_lm_stat.pO2_v_median]...
            = fun_analysis_get_mean_N_median_in_each_cell(tmp_cube_data.pO2_lm.pO2_v);
        tmp_region_data{iter_cell} = tmp_cube_data;
    end
    
    tmp_region_dt_v_mean = cellfun(@(x) x.dt_lm_stat.dt_v_mean', tmp_region_data, ...
        'UniformOutput', false);
    tmp_region_dt_v_mean = cat(1, tmp_region_dt_v_mean{:});
    
    tmp_region_dt_v_mean_stat = fun_analysis_get_basic_statistics_in_column(tmp_region_dt_v_mean);
    %%    
    tmp_diff_dt_v_mean = diff(tmp_region_dt_v_mean, 1, 2);
    tmp_diff_dt_v_mean = tmp_diff_dt_v_mean ./ diff_window_size;
    tmp_diff_stat = fun_analysis_get_basic_statistics_in_column(tmp_diff_dt_v_mean);
    %%
    [ax_hdl, plt_hdl, ~] = fun_vis_confidence_interval_shaded(lm_window_size_list, ...
        tmp_region_dt_v_mean_stat.median, tmp_region_dt_v_mean_stat.prctile_val(:, 7), ...
        tmp_region_dt_v_mean_stat.prctile_val(:, 9), ax_hdl);
%     plot(ax_hdl, lm_window_size_list, tmp_region_dt_v_mean_stat.mean, 'LineWidth', 1);
    hold(ax_hdl, 'on');    
%     [ax_hdl_2, plt_hdl, patch_hdl] = fun_vis_confidence_interval_shaded(tmp_x, ...
%         tmp_diff_stat.median, tmp_diff_stat.prctile_val(:, 7), ...
%         tmp_diff_stat.prctile_val(:, 9), ax_hdl_2);
    plt_hdl = plot(ax_hdl_2, tmp_x, tmp_diff_stat.mean, 'LineWidth', 1, 'Color', plt_hdl.Color);
    hold(ax_hdl_2, 'on');
    
    leg_hdl(iter_region) = plt_hdl;
end
ax_hdl_2.XLim = ax_hdl.XLim;
legend(ax_hdl_2, leg_hdl, fig_legend{1:3});

ax_hdl_2.XLabel.String = 'Window size (\mum)';
ax_hdl_2.YLabel.String = 'Derivative';
ax_hdl.YLabel.String = 'd_{lm} (\mum)';
ax_hdl_2.YScale = 'log';
fig_fp = fullfile(vis_folder, sprintf('%s_%s_dt_lm_vs_window_size_SS_IC_PI.png', ...
    dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
stack_name_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
grid_version = '240_cube';
save_folder_name = 'pO2_SA_li_sc';

skel_version = '240_cube_rec';
reconstruction_version = '240_cube_recon_sc';
wb_stat_folder_name = 'whole_brain_stat_sc';
length_density_file_postfix = 'cube_length_density_sc';

Allen_atlas_id = load('Allen_atlas_id.mat');
merge_stack_name = 'all_stack';

vis_folder_name = 'pO2_dt_SA_scaling';
save_im_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, merge_stack_name),...
    vis_folder_name);
if ~isfolder(save_im_folder)
    mkdir(save_im_folder);
end

pO2_SA = struct;
pO2_SA.dataset_name = dataset_name;
pO2_SA.pO2_version = save_folder_name;
pO2_SA.skel_version = skel_version;
pO2_SA.reconstruction_version = reconstruction_version;
pO2_SA.stack_list = stack_list;
pO2_SA.filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, ...
    merge_stack_name), 'pO2_lm_SA_result_ptl.mat');

output_dtlm_vs_dulm_figure_Q = true;
output_lm_vs_ws_figure_Q = true;
%% Load preprocessed data
num_stack = numel(stack_list);
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : numel(stack_list)
    stack = stack_list{iter_stack};
    tmp_pO2_stat_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack), ...
        sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, stack, skel_version, save_folder_name));
    
    tmp_len_den_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack), ...
        sprintf('%s_%s_%s_%s.mat', dataset_name, stack, skel_version, length_density_file_postfix));
    
    tic_start = tic;
    wb_region_stat_data.dataset_name = dataset_name;
    wb_region_stat_data.stack = stack;
    wb_region_stat_data.grid_version = grid_version;
    wb_region_stat_data.skeleton_version = skel_version;
    wb_region_stat_data.cube_stat = DataManager.load_analysis_data(dataset_name, ...
        stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
        dataset_name, stack, reconstruction_version), ...
        wb_stat_folder_name);
    wb_region_stat_data.pO2_data = DataManager.load_data(tmp_pO2_stat_filepath);
    wb_region_stat_data.length_density_data = DataManager.load_data(tmp_len_den_filepath);
    wb_data_cell{iter_stack} = wb_region_stat_data;
    fprintf('Finish loading whole brain pO2 simulation result in %s. Elapsed time is %f seconds.\n', ...
        stack, toc(tic_start));
end
clearvars wb_region_stat_data;

lm_wz_um = wb_data_cell{1}.pO2_data.local_extrema_window_size;
lm_wz_um_mean = movmean(lm_wz_um, 2, 'Endpoints', 'discard');
lm_wz_um_diff = diff(lm_wz_um);
num_window_size = numel(lm_wz_um);
%% Use the same criteria for selecting the 240-cube for all the analysis 
min_cap2vsl_vol_fraction = 0.5;
is_selected_cube_cell = cell(num_stack, 1);
is_internal_cube_cell = cell(num_stack, 1);

upper_prctile = 99.5;
lower_prctile = 0.5;
for iter_stack = 1 : num_stack
    tmp_stack_data = wb_data_cell{iter_stack}.cube_stat;
    tmp_is_internal_Q = tmp_stack_data.cube_in_brain_mask_ratio == 1;
    % Focus on cubes dominanted by capillaries 
    tmp_primary_selection_Q = tmp_is_internal_Q & ...    
        tmp_stack_data.cap2vsl_vol_fraction >= min_cap2vsl_vol_fraction;
    % Further selection by removing outliers
    tmp_secondary_selection_Q = ...
        fun_analysis_is_in_percentile_interval_Q(tmp_stack_data.link_all_stat.mean.nearest_tissue_dt_max, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
        fun_analysis_is_in_percentile_interval_Q(tmp_stack_data.link_all_stat.num_data.length, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
        fun_analysis_is_in_percentile_interval_Q(tmp_stack_data.link_length_density_m_mm3, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
        fun_analysis_is_in_percentile_interval_Q(wb_data_cell{iter_stack}.pO2_data.local_dt_stat.mean, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
        fun_analysis_is_in_percentile_interval_Q(wb_data_cell{iter_stack}.pO2_data.local_pO2_stat.mean, lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
        fun_analysis_is_in_percentile_interval_Q(wb_data_cell{iter_stack}.pO2_data.pO2_lm.pO2_mean(:, 11), lower_prctile, upper_prctile, tmp_primary_selection_Q) & ...
        fun_analysis_is_in_percentile_interval_Q(wb_data_cell{iter_stack}.pO2_data.dt_lm.dt_mean(:, 11), lower_prctile, upper_prctile, tmp_primary_selection_Q);
    is_selected_cube_cell{iter_stack} = tmp_secondary_selection_Q;
    is_internal_cube_cell{iter_stack} = tmp_is_internal_Q;
end
is_selected_cube_Q = cat(1, is_selected_cube_cell{:});
is_internal_cube_Q = cat(1, is_internal_cube_cell{:});

vsl_len_den = cellfun(@(x) x.cube_stat.link_length_density_m_mm3, wb_data_cell, 'UniformOutput', false);
vsl_len_den = cat(1, vsl_len_den{:});
histogram(vsl_len_den(is_selected_cube_Q));

pO2_SA.is_internal_cube_Q = is_internal_cube_cell;
pO2_SA.is_selected_cube_Q = is_selected_cube_cell;
pO2_SA.min_cap2vsl_vol_fraction = min_cap2vsl_vol_fraction;
pO2_SA.selection_upper_prctile = upper_prctile;
pO2_SA.selection_lower_prctile = lower_prctile;

avg_r = cellfun(@(x) x.cube_stat.link_cap_stat.mean.dt_median, wb_data_cell, 'UniformOutput', false);
avg_r = cat(1, avg_r{:});
pO2_SA.avg_cap_r_um = mean(avg_r(is_selected_cube_Q), 'omitnan');
%% Questions to address:
% 1. Whick window size should be used to show the result?
% 2. How does the fitting coefficient, and eventually pO2 drop depdends on
% the window size?
% 3. Why does the local extrema work while the attempt to correct the
% entire curve fail?
% 4.
%% Overall local extrema vs window size
lm_vs_ws_info = struct;
vis_fn_l1 = {'dt_lm', 'pO2_lm'};
vis_info_cell = cell(2, 0);
vis_info_cell(:, end+1) = {'dt_mean', 'Average d_{lm} (\mum)'};
vis_info_cell(:, end+1) = {'dt_mean_prctile', '1 - CDF(Average d_{lm})'};
vis_info_cell(:, end+1) = {'dt_median', 'Median d_{lm} (\mum)'};
vis_info_cell(:, end+1) = {'dt_median_prctile', '1 - CDF(Median d_{lm})'};
vis_info_cell(:, end+1) = {'pO2_mean', 'Average u_{lm}'};
vis_info_cell(:, end+1) = {'pO2_mean_prctile', 'CDF(Average u_{lm})'};
vis_info_cell(:, end+1) = {'pO2_median', 'Median u_{lm}'};
vis_info_cell(:, end+1) = {'pO2_median_prctile', 'CDF(Median u_{lm})'};
for iter_fn = 1 : numel(vis_fn_l1)
    tmp_fn = vis_fn_l1{iter_fn};
    for iter_vis = 1 : size(vis_info_cell, 2)
        tmp_vis_data_n = vis_info_cell{1, iter_vis};
        tmp_vis_data_yl = vis_info_cell{2, iter_vis};
        
        tmp_data = cellfun(@(x) x.pO2_data.(tmp_fn).(tmp_vis_data_n), wb_data_cell, 'UniformOutput', false);
        tmp_data = cat(1, tmp_data{:});
        tmp_data = tmp_data(is_selected_cube_Q, :);
        if contains(tmp_vis_data_n, 'prctile') && contains(tmp_vis_data_n, 'dt')
            tmp_data = 1 - tmp_data;
        end
        tmp_data_diff = diff(tmp_data, 1, 2) ./ lm_wz_um_diff;
        tmp_data_stat = fun_analysis_get_basic_statistics_in_column(tmp_data);
        tmp_data_diff_stat = fun_analysis_get_basic_statistics_in_column(tmp_data_diff);
        
        lm_vs_ws_info.(tmp_fn).(tmp_vis_data_n).stat = tmp_data_stat;
        lm_vs_ws_info.(tmp_fn).(tmp_vis_data_n).derivative_stat = tmp_data_diff_stat;
        
        vis_prctile_low = 7;
        vis_prctile_high = 9;
        %% Generate figures
        if output_lm_vs_ws_figure_Q
            fig_hdl = figure;
            ax_hdl_1 = axes(fig_hdl);
            yyaxis('left');
            [ax_hdl_1, plt_hdl_1, ~] = fun_vis_confidence_interval_shaded(lm_wz_um, ...
                tmp_data_stat.median, tmp_data_stat.prctile_val(:, vis_prctile_low), tmp_data_stat.prctile_val(:, vis_prctile_high), ax_hdl_1);
            ax_hdl_1.YLabel.String = tmp_vis_data_yl;
            if contains(tmp_vis_data_n, 'dt') && ~contains(tmp_vis_data_n, 'prctile')
                ax_hdl_1.YLim(1) = 0;
            elseif contains(tmp_vis_data_n, 'pO2') && ~contains(tmp_vis_data_n, 'prctile')
                ax_hdl_1.YLim(2) = 0;
            end
            
            yyaxis('right');
            [ax_hdl_1, plt_hdl_2, ~] = fun_vis_confidence_interval_shaded(lm_wz_um_mean, ...
                tmp_data_diff_stat.median, tmp_data_diff_stat.prctile_val(:, vis_prctile_low), tmp_data_diff_stat.prctile_val(:, vis_prctile_high), ax_hdl_1);
            ax_hdl_1.YLabel.String = 'Derivative';
            ax_hdl_1.XLabel.String = 'Window size, w (\mum)';
            % ax_hdl_1.YScale = 'log';
            grid(ax_hdl_1, 'on');
            box(ax_hdl_1, 'on');
            ax_hdl_1.FontSize = 14;
            ax_hdl_1.Title.String = upper(strrep(tmp_fn, '_', ' '));
            fig_fp = fullfile(save_im_folder, 'lm_vs_wz', sprintf('%s_%s_%s_%s_vs_lm_window_size.png', ...
                dataset_name, merge_stack_name, tmp_fn, tmp_vis_data_n));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            delete(fig_hdl);
        end
    end
end
pO2_SA.lm_vs_wz = lm_vs_ws_info;
%% Result 1 - determination of length scale for metabolism plot
% 1. pO2LM DT increase slower than DTLM DT - is there a window size where
% two DT stat cross? Yes, 50 um. Do know why though.
% 2. 50 um seems to be the scale where the derivative of DTLM DT flats out.
% But pO2LM seems does not have any length scale where derivative flats
% out. - wired.
% 3. So, use 50 um.
if output_lm_vs_ws_figure_Q
    fig_hdl = figure;
    ax_hdl = subplot(2, 1, 1);
    yyaxis(ax_hdl, 'left');
    plot(ax_hdl, lm_wz_um, pO2_SA.lm_vs_wz.dt_lm.dt_mean.stat.mean, 'LineWidth', 1);
    ax_hdl.YLabel.String = 'Average d_{lm} (\mum)';
    ax_hdl.YLim = [20, 40];
    hold(ax_hdl, 'on');
    yyaxis(ax_hdl, 'right');
    plot(ax_hdl, lm_wz_um, pO2_SA.lm_vs_wz.pO2_lm.dt_mean.stat.mean, 'LineWidth', 1);
    ax_hdl.YLabel.String = 'Average d_{ulm} (\mum)';
    ax_hdl.YLim = [20, 40];
    ax_hdl.XLabel.String = 'Window size (\mum)';
    grid(ax_hdl, 'on');
    ax_hdl = subplot(2, 1, 2);
    yyaxis(ax_hdl, 'left');
    plot(ax_hdl, lm_wz_um, pO2_SA.lm_vs_wz.dt_lm.dt_mean.stat.mean - pO2_SA.lm_vs_wz.pO2_lm.dt_mean.stat.mean, 'LineWidth', 1);
    ax_hdl.YLabel.String = 'Average d_{lm} - Average d_{ulm} (\mum)';
    hold(ax_hdl, 'on');
    yyaxis(ax_hdl, 'right');
    plot(ax_hdl, lm_wz_um, abs(pO2_SA.lm_vs_wz.dt_lm.dt_mean.stat.mean - pO2_SA.lm_vs_wz.pO2_lm.dt_mean.stat.mean) ./ pO2_SA.lm_vs_wz.pO2_lm.dt_mean.stat.mean, 'LineWidth', 1);
    ax_hdl.YLabel.String = '|d_{lm} - d_{ulm}|/d_{ulm}';
    grid(ax_hdl, 'on');
    ax_hdl.YScale = 'log';
    ax_hdl.XLabel.String = 'Window size (\mum)';
    fig_fp = fullfile(save_im_folder, 'lm_vs_wz', sprintf('%s_%s_all_avg_diff_d_lm_d_ulm_vs_lm_window_size.png', ...
        dataset_name, merge_stack_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
end
%% Correlation between DTLM DT and pO2LM DT
% As we are not really using the slope for pO2 calculation, here we use all
% the internal 240-cubes for fitting - to justify the assumption /
% approximation roughly hold for the whole brain. 
% Another justification for choosing window size 50 um
have_intercept_Q = false;
plot_info_cell = cell(6, 0);
plot_info_cell(:, end+1) = {'dt_lm.dt_median', 'pO2_lm.dt_median', 'Median $d_{lm}\;(\mu m)$', 'Median $d_{ulm}\;(\mu m)$', 'med_lm_dt', 'med_pO2_lm_dt'};
plot_info_cell(:, end+1) = {'dt_lm.dt_mean', 'pO2_lm.dt_mean', 'Mean $d_{lm}\;(\mu m)$', 'Mean $d_{ulm}\;(\mu m)$', 'mean_lm_dt', 'mean_pO2_lm_dt'};

num_plot = size(plot_info_cell, 2);
dulm_vs_ulm_fit_info = struct;
%%  Generate figures
for iter_plot = 1 : num_plot
    tmp_x_label = plot_info_cell{3, iter_plot};
    tmp_y_label = plot_info_cell{4, iter_plot};
    tmp_x_name = plot_info_cell{5, iter_plot};
    tmp_y_name = plot_info_cell{6, iter_plot};
    
    [tmp_x_data, tmp_y_data] = deal(cell(num_stack, 1));
    for iter_stack = 1 : num_stack
        tmp_x_data{iter_stack} = fun_getfield(wb_data_cell{iter_stack}.pO2_data, plot_info_cell{1, iter_plot});
        tmp_y_data{iter_stack} = fun_getfield(wb_data_cell{iter_stack}.pO2_data, plot_info_cell{2, iter_plot});
    end
    tmp_x_data = cat(1, tmp_x_data{:});
    tmp_y_data = cat(1, tmp_y_data{:});
    %%
    tmp_fit_info_cell = cell(num_window_size, 1);
    for vis_wz_idx = 1 : num_window_size
        vis_wz_um = lm_wz_um(vis_wz_idx);
        %%
        tmp_vis_x_data = tmp_x_data(:, vis_wz_idx);
        tmp_vis_y_data = tmp_y_data(:, vis_wz_idx);
        tmp_valid_data_Q = isfinite(tmp_vis_x_data) & isfinite(tmp_vis_y_data) & is_internal_cube_Q;
        tmp_num_data_before_selection = nnz(tmp_valid_data_Q);
        tmp_valid_data_Q = tmp_valid_data_Q & is_selected_cube_Q;
        tmp_vis_x_data = tmp_vis_x_data(tmp_valid_data_Q);
        tmp_vis_y_data = tmp_vis_y_data(tmp_valid_data_Q);
        tmp_num_data_after_selection = nnz(tmp_valid_data_Q);
        if have_intercept_Q
            linear_fit_hdl = fitlm(tmp_vis_x_data, tmp_vis_y_data);
        else
            linear_fit_hdl = fitlm(tmp_vis_x_data, tmp_vis_y_data, 'Intercept', false);
        end
        
        tmp_fit_info_str = struct;
        tmp_fit_info_str.window_size_um = vis_wz_um;
        tmp_fit_info_str.num_data_point = tmp_num_data_after_selection;
        tmp_fit_info_str.fraction_of_valid_cube = tmp_num_data_after_selection / tmp_num_data_before_selection;
        tmp_fit_info_str.Estimate = linear_fit_hdl.Coefficients.Estimate;
        tmp_fit_info_str.SE = linear_fit_hdl.Coefficients.SE;
        tmp_fit_info_str.RSquaredAdj = linear_fit_hdl.Rsquared.Adjusted;
        tmp_fit_info_cell{vis_wz_idx} = tmp_fit_info_str;
        fit_data_x = linspace(0, max(tmp_vis_x_data, [], 'all'), 30);
        %% Visualization
        if output_dtlm_vs_dulm_figure_Q
            fig_hdl = figure;
            ax_hdl = axes(fig_hdl);
            histogram2(ax_hdl, tmp_vis_x_data, tmp_vis_y_data, 'DisplayStyle', 'tile');
            if contains(tmp_x_label, '$')
                ax_hdl.XLabel.Interpreter = 'latex';
            end
            if contains(tmp_y_label, '$')
                ax_hdl.YLabel.Interpreter = 'latex';
            end
            hold(ax_hdl, 'on');
            ax_hdl.XLabel.String = tmp_x_label;
            ax_hdl.YLabel.String = tmp_y_label;
            ax_hdl.Title.String = sprintf('Window size %d \\mum', vis_wz_um);
            ax_hdl.FontSize = 14;
            
            if have_intercept_Q
                fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
                    'LineWidth', 3);
            else
                fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(1), ...
                    'LineWidth', 3);
            end
            cbar_hdl = colorbar(ax_hdl);
            cbar_hdl.Label.String = 'Number of data points';
            
            leg_string = sprintf('Slope: %.3e \\pm %.1e\nR^2-Adjusted: %.2f\nData size: %d (%.4f)', ...
                tmp_fit_info_str.Estimate(1), tmp_fit_info_str.SE(1), ...
                tmp_fit_info_str.RSquaredAdj, tmp_fit_info_str.num_data_point, ...
                tmp_fit_info_str.fraction_of_valid_cube);
            ax_hdl.ColorScale = 'log';
            ax_hdl.CLim(1) = 1;
            leg_hdl = legend(ax_hdl, fit_plt_hdl, leg_string, 'Location', 'southeast');
            
            fig_fp = fullfile(save_im_folder, 'lm_dt_fit_vs_wz', sprintf('%s_%s_%s_vs_%s_wz_%d_um.png', ...
                dataset_name, merge_stack_name, tmp_x_name, tmp_y_name, vis_wz_um));
            fun_print_image_in_several_formats(fig_hdl, fig_fp);
            delete(fig_hdl);
        end
    end
    dulm_vs_ulm_fit_info.(tmp_x_name) = cat(1, tmp_fit_info_cell{:});
end
pO2_SA.dulm_vs_dlm = dulm_vs_ulm_fit_info;
%% Slope vs window size
if output_dtlm_vs_dulm_figure_Q
    tmp_slope = [pO2_SA.dulm_vs_dlm.mean_lm_dt.Estimate];
    tmp_dsdw = diff(tmp_slope) ./  lm_wz_um_diff;
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
    yyaxis(ax_hdl, 'left');
    plot(ax_hdl, lm_wz_um, tmp_slope, 'LineWidth', 1);
    ax_hdl.XLabel.String = 'Window size (\mum)';
    ax_hdl.YLabel.String = 'Slope';
    hold(ax_hdl, 'on');
    yyaxis(ax_hdl, 'right');
    plot(ax_hdl, lm_wz_um, [pO2_SA.dulm_vs_dlm.mean_lm_dt.RSquaredAdj], 'LineWidth', 1);
    ax_hdl.YLabel.String = 'R^2-Adjusted';
    grid(ax_hdl, 'on');
    ax_hdl.YLim = [0, 1];
    ax_hdl.FontSize = 14;
    fig_fp = fullfile(save_im_folder, 'lm_dt_fit_vs_wz', sprintf('%s_%s_mean_dt_lm_fit_vs_wd_sz.png', ...
        dataset_name, merge_stack_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
end
%% Scaling between local DT properties and length density
% The question here is how to choice the 240-cubes, and if we want to use
% the fitting intercept as the vessel radius in the Krogh model? 
lattice_space_scaling_str = load('./Simulation/space_filling_lattice/lattice_space_filling_scaling.mat');
have_intercept_Q = true;
num_lattice = numel(lattice_space_scaling_str.lattice_name);

plot_info_cell = cell(4, 0);
% plot_info_cell(:, end+1) = {'link_cap_stat.mean.nearest_tissue_dt_max', 'nearest_tissue_dt_max_mean', 'Average vessel $d_{max}\;(\mu m)$', 'avg_tissue_dt_max'};
% plot_info_cell(:, end+1) = {'link_cap_stat.median.nearest_tissue_dt_max', 'nearest_tissue_dt_max_median', 'Median vessel $d_{max}\;(\mu m)$', 'med_tissue_dt_max'};
%
% plot_info_cell(:, end+1) = {'link_cap_stat.mean.nearest_tissue_dt_mean', 'nearest_tissue_dt_mean_mean', 'Average vessel $d_{avg}\;(\mu m)$', 'avg_tissue_dt_max'};
% plot_info_cell(:, end+1) = {'link_cap_stat.median.nearest_tissue_dt_mean', 'nearest_tissue_dt_mean_median', 'Median vessel $d_{avg}\;(\mu m)$', 'med_tissue_dt_max'};
%
% plot_info_cell(:, end+1) = {'link_cap_stat.mean.nearest_tissue_dt_median', 'nearest_tissue_dt_median_mean', 'Average vessel $d_{med}\;(\mu m)$', 'avg_tissue_dt_med'};
% plot_info_cell(:, end+1) = {'link_cap_stat.median.nearest_tissue_dt_median', 'nearest_tissue_dt_median_median', 'Median vessel $d_{med}\;(\mu m)$', 'med_tissue_dt_med'};

plot_info_cell(:, end+1) = {'dt_lm.dt_mean', 'nearest_tissue_dt_max_mean', 'Average $d_{lm}\;(\mu m)$', 'avg_dlm'};
plot_info_cell(:, end+1) = {'dt_lm.dt_median', 'nearest_tissue_dt_max_median', 'Median $d_{lm}\;(\mu m)$', 'med_dlm'};

plot_info_cell(:, end+1) = {'pO2_lm.dt_mean', 'nearest_tissue_dt_max_mean', 'Average $d_{ulm}\;(\mu m)$', 'avg_dulm'};
plot_info_cell(:, end+1) = {'pO2_lm.dt_median', 'nearest_tissue_dt_max_median', 'Median $d_{ulm}\;(\mu m)$', 'med_dulm'};

plot_info_cell(:, end+1) = {'local_dt_stat.mean', 'nearest_tissue_dt_mean_mean', '$\bar{d}\;(\mu m)$', 'avg_d'};
plot_info_cell(:, end+1) = {'local_dt_stat.median', 'nearest_tissue_dt_median_mean', 'Median $d \;(\mu m)$', 'med_d'};

% plot_x_field_name = 'capillary_length_density_m_mm3';
% plot_x_label_name = 'Capillary $\rho_l^{-1/2}\; (\mu m)$';
% plot_output_var_name_x = 'isqrt_cap_length_density';
% str_field_name = 'cap_den_2_d';

plot_x_field_name = 'cube_stat.link_length_density_m_mm3';
% plot_x_field_name = 'length_density_data.all_w_eph';
% plot_x_field_name = 'length_density_data.all';
plot_x_label_name = 'Vessel $\rho_l^{-1/2}\; (\mu m)$';
plot_output_var_name_x = 'isqrt_link_length_density';
str_field_name = 'vsl_den_2_d';
rho_to_d_fit = struct;
rho_to_d_fit.x_label_name = plot_x_label_name;

num_fitting = size(plot_info_cell, 2);
    %% Generate figures
for iter_fitting = 1 : num_fitting
    %% Extract data
    tmp_y_field_name = plot_info_cell{1, iter_fitting};
    tmp_y_name = plot_info_cell{2, iter_fitting};
    tmp_y_label = plot_info_cell{3, iter_fitting};
    tmp_output_var_name = plot_info_cell{4, iter_fitting};
    
    tmp_y_fitting_data = lattice_space_scaling_str.(tmp_y_name);
    
    tmp_y_slope = tmp_y_fitting_data.slope;
    [tmp_y_slope, tmp_sorted_idx] = sort(tmp_y_slope, 'descend');
    tmp_lattice_list = lattice_space_scaling_str.lattice_name;
    tmp_lattice_list = tmp_lattice_list(tmp_sorted_idx);
    
    [tmp_x_data, tmp_y_data, ] = deal(cell(num_stack, 1));
    for iter_stack = 1 : num_stack
        tmp_stack_y = fun_getfield(wb_data_cell{iter_stack}.pO2_data, tmp_y_field_name);        
        tmp_stack_data = wb_data_cell{iter_stack}.cube_stat;
        tmp_stack_x = fun_getfield(wb_data_cell{iter_stack}, plot_x_field_name);
        tmp_stack_x = (tmp_stack_x * 1e-3).^ (-1/2);
        tmp_y_data{iter_stack} = tmp_stack_y;
        tmp_x_data{iter_stack} = tmp_stack_x;
    end    
    tmp_x_data = cat(1, tmp_x_data{:});
    tmp_y_data = cat(1, tmp_y_data{:});
    tmp_num_y_column = size(tmp_y_data, 2);
    %% 
    tmp_fit_info_cell = cell(num_window_size, 1);
    for vis_wz_idx = 1 : tmp_num_y_column
        vis_wz_um = lm_wz_um(vis_wz_idx);

        tmp_vis_x = tmp_x_data;
        if tmp_num_y_column == 1
            tmp_vis_y = tmp_y_data;
        else
            tmp_vis_y = tmp_y_data(:, vis_wz_idx);
        end        
        tmp_valid_data_Q = isfinite(tmp_vis_x) & isfinite(tmp_vis_y) & is_internal_cube_Q;
        tmp_num_data_before_selection = nnz(tmp_valid_data_Q);
        tmp_valid_data_Q = tmp_valid_data_Q & is_selected_cube_Q;
        
        tmp_vis_x = tmp_vis_x(tmp_valid_data_Q);
        tmp_vis_y = tmp_vis_y(tmp_valid_data_Q);
        tmp_num_data_after_selection = nnz(tmp_valid_data_Q);
        
        if have_intercept_Q
            linear_fit_hdl = fitlm(tmp_vis_x, tmp_vis_y);
        else
            linear_fit_hdl = fitlm(tmp_x_data, tmp_vis_y, 'Intercept', false);
        end        
        tmp_fit_info_str = struct;
        tmp_fit_info_str.window_size_um = vis_wz_um;
        tmp_fit_info_str.num_data_point = tmp_num_data_after_selection;
        tmp_fit_info_str.fraction_of_valid_cube = tmp_num_data_after_selection / tmp_num_data_before_selection;
        tmp_fit_info_str.Estimate = linear_fit_hdl.Coefficients.Estimate.';
        tmp_fit_info_str.SE = linear_fit_hdl.Coefficients.SE.';
        tmp_fit_info_str.RSquaredAdj = linear_fit_hdl.Rsquared.Adjusted;
        tmp_fit_info_cell{vis_wz_idx} = tmp_fit_info_str;
        %% Visualization
        fig_hdl = figure('Visible', 'on');
        fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.4;
        ax_hdl = axes(fig_hdl);
        histogram2(ax_hdl, tmp_vis_x, tmp_vis_y, 'DisplayStyle', 'tile');
        ax_hdl.XLim(1) = 0;
%         ax_hdl.YLim(2) = round(ax_hdl.YLim(2) / 5) * 5;
%         ax_hdl.XLim(2) = round(max(tmp_vis_x) / 5) * 5;
        hold(ax_hdl, 'on');
        fit_data_x = linspace(0, max(tmp_vis_x, [], 'all'), 30);
        leg_hdl_array = [];
        if have_intercept_Q
            fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
                'LineWidth', 3, 'Color', 'k');
            ax_hdl.YLim(1) = linear_fit_hdl.Coefficients.Estimate(1);
        else
            fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(1), ...
                'LineWidth', 3, 'Color', 'k');
            ax_hdl.YLim(1) = 0;
        end
        leg_hdl_array(end+1) = fit_plt_hdl;
       
        for iter_lattice = 1 : num_lattice
            if have_intercept_Q
                tmp_lattice_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* tmp_y_slope(iter_lattice) + ...
                    linear_fit_hdl.Coefficients.Estimate(1), 'LineWidth', 1.5, 'LineStyle', '-');
            else
                tmp_lattice_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* tmp_y_slope(iter_lattice), 'LineWidth', 3, 'LineStyle', '-.');
            end
            leg_hdl_array(end+1) = tmp_lattice_hdl;
        end
        if have_intercept_Q
            fit_legend = sprintf('Slope: %.3f \\pm %.3f\nIntercept: %.2f \\pm %.2f\nR^2-Adjusted: %.2f\nData size: %d (%d%%)', ...
                linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
                linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
                linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/tmp_num_data_before_selection));
        else
            fit_legend = sprintf('Slope: %.3f \\pm %.3f\nR^2-Adjusted: %.2f\nData size: %d (%d%%)', ...
                linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
                linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/tmp_num_data_before_selection));
        end
        leg_str_array = cat(1, {fit_legend}, cellfun(@(x1, x2) sprintf('%s lattice: %.3f', x1, x2), ...
            tmp_lattice_list', num2cell(tmp_y_slope), 'UniformOutput', false));
        
        leg_hdl = legend(leg_hdl_array, leg_str_array, 'Location', 'northwest', 'Box', true);
        
        cbar_hdl = colorbar;
        cbar_hdl.Label.String = 'Number of data points';
        ax_hdl.FontSize = 14;
        ax_hdl.XLabel.Interpreter = 'latex';
        ax_hdl.XLabel.String = plot_x_label_name;
        ax_hdl.YLabel.Interpreter = 'latex';
        ax_hdl.YLabel.String = tmp_y_label;
        ax_hdl.XLabel.FontSize = 18;
        ax_hdl.YLabel.FontSize = 18;
        ax_hdl.ColorScale = 'log';
        ax_hdl.CLim(1) = 1;
        if tmp_num_y_column > 1
            ax_hdl.Title.String = sprintf('Window size %d \\mum', vis_wz_um);
            tmp_fp = fullfile(save_im_folder, 'd_vs_rho', tmp_output_var_name, sprintf('%s_%s_%s_vs_%s_mc2v_vr%.2f_wd_sz_%d_um.png', ...
                dataset_name, merge_stack_name, tmp_output_var_name,...
                plot_output_var_name_x, min_cap2vsl_vol_fraction, vis_wz_um));
        else
            tmp_fp = fullfile(save_im_folder, 'd_vs_rho', sprintf('%s_%s_%s_vs_%s_mc2v_vr%.2f.png', ...
                dataset_name, merge_stack_name, tmp_output_var_name,...
                plot_output_var_name_x, min_cap2vsl_vol_fraction));            
        end        
        %%
        fun_print_image_in_several_formats(fig_hdl, tmp_fp);
        delete(fig_hdl);
    end
    rho_to_d_fit.(tmp_output_var_name) = cat(1, tmp_fit_info_cell{:});
end
pO2_SA.(str_field_name) = rho_to_d_fit;
%% 
r_cap = pO2_SA.avg_cap_r_um;
krogh_coeff = 1/2;
%%
fit_info_cell = cell(5, 0);
fit_info_cell(:, end+1) = {'pO2_data.pO2_lm.dt_mean', 'pO2_data.pO2_lm.pO2_mean', 'Average d_{ulm} (\mum)', 'Average u_{lm}', 'dulm_2_ulm'};
fit_info_cell(:, end+1) = {'pO2_data.dt_lm.dt_mean', 'pO2_data.pO2_lm.pO2_mean', 'Average d_{lm} (\mum)', 'Average u_{lm}', 'dlm_2_ulm'};
fit_info_cell(:, end+1) = {'pO2_data.dt_lm.dt_mean', 'pO2_data.dt_lm.pO2_mean', 'Average d_{lm} (\mum)', 'Average u_{dlm}', 'dlm_2_udlm'};
num_fit = size(fit_info_cell, 2);
%% Try to fit Krogh to local pO2 minimum vs its corresponding dt
scaled_krogh_str = struct;
for iter_fit = 1 : num_fit
    tmp_x_field = fit_info_cell{1, iter_fit};
    tmp_y_field = fit_info_cell{2, iter_fit};
    tmp_x_label = fit_info_cell{3, iter_fit};
    tmp_y_label = fit_info_cell{4, iter_fit};
    tmp_folder_name = fit_info_cell{5, iter_fit};
    
    tmp_x = cellfun(@(x) fun_getfield(x, tmp_x_field), wb_data_cell, 'UniformOutput', false);
    tmp_x = cat(1, tmp_x{:});
    tmp_y = cellfun(@(x) fun_getfield(x, tmp_y_field), wb_data_cell, 'UniformOutput', false);
    tmp_y = cat(1, tmp_y{:});
    tmp_fit_info_cell = cell(num_window_size, 1);
    for iter_wd_sz = 1 : num_window_size
        % Visualization for each window size
        tmp_wd_sz = lm_wz_um(iter_wd_sz);
        
        vis_x = tmp_x(:, iter_wd_sz);
        vis_y = tmp_y(:, iter_wd_sz);
        tmp_valid_Q = isfinite(vis_x) & isfinite(vis_y) & is_internal_cube_Q;
        tmp_num_data_before_selection = nnz(tmp_valid_Q);
        tmp_valid_Q = is_selected_cube_Q & tmp_valid_Q;
        vis_x = vis_x(tmp_valid_Q);
        vis_y = vis_y(tmp_valid_Q);
        fit_x_fun = @(x) x + r_cap;
        fit_x = vis_x + r_cap;
        fit_x_T = - krogh_coeff * (fit_x .^ 2 .* (log(fit_x / r_cap) - 1/2) + r_cap^2/2);
        linear_fit_hdl = fitlm(fit_x_T, vis_y, 'Intercept', false);
        
        tmp_fit_info_str = struct;
        tmp_fit_info_str.window_size_um = tmp_wd_sz;
        tmp_fit_info_str.num_data_point = nnz(tmp_valid_Q);
        tmp_fit_info_str.fraction_of_valid_cube = tmp_fit_info_str.num_data_point / tmp_num_data_before_selection;
        tmp_fit_info_str.Estimate = linear_fit_hdl.Coefficients.Estimate.';
        tmp_fit_info_str.SE = linear_fit_hdl.Coefficients.SE.';
        tmp_fit_info_str.RSquaredAdj = linear_fit_hdl.Rsquared.Adjusted;
        tmp_fit_info_cell{iter_wd_sz} = tmp_fit_info_str;
        %% Quadratic decay: 
%         fun_getfields(linear_fit_hdl, {'Coefficients.Estimate', 'Coefficients.SE', 'Rsquared.Adjusted'})
%         tmp_quadric_fun = @(x) x.^2;
%         fit_x_T = tmp_quadric_fun(fit_x);
%         tmp_quadric_fit = fitlm(fit_x_T, vis_y, 'Intercept', false);
%         fig_hdl = figure;
%         ax_hdl = axes(fig_hdl);
%         histogram2(ax_hdl, fit_x_T, vis_y, 'DisplayStyle', 'tile');
%         hold(ax_hdl, 'on');
%         tmp_plot_x = tmp_quadric_fun(0 : max(fit_x));
%         plot(ax_hdl, tmp_plot_x, tmp_plot_x .* tmp_quadric_fit.Coefficients.Estimate(1));
        %% Pure phenomelogical fitting - nonlinear fit
        % Directly on ulm vs dulm 
        % Optimized coefficient: [0.1290, -1.7702]
%         tmp_ratio = vis_y ./ fit_x.^2;
%         tmp_ratio_in_dmax_bin = fun_analysis_get_y_stat_in_x_bin(fit_x, tmp_ratio, [2 : 0.5 : 60]);
%         tmp_is_ol_Q = tmp_ratio < -0.7 | fit_x > 45;
% 
%         tmp_nl_fun = @(a, x) - a(1) .* x .^ 2 .* ( log(x ./ r_cap) + a(2) .* (1 - r_cap^2 ./ (x.^2)));
%         tmp_est_ab = [0.26, -0.5];
%         tmp_opt_coeff_full = nlinfit(fit_x, vis_y, tmp_nl_fun, tmp_est_ab);
%         tmp_opt_cov = corr(tmp_nl_fun(tmp_opt_coeff_full, fit_x), vis_y);
%         % Fit the dimensionless function
%         % nlinfit gives the same fitting result as linear regression. 
%         tmp_nldl_fun = @(a, x) - a(1) .* ( log(x ./ r_cap) + a(2) .* (1 - r_cap.^2 ./ (x .^ 2)));
%         tmp_opt_coeff_dl = nlinfit(fit_x(~tmp_is_ol_Q), tmp_ratio(~tmp_is_ol_Q), tmp_nldl_fun, tmp_est_ab);
%         tmp_opt_dl_cov = corr(tmp_nldl_fun(tmp_opt_coeff_dl, fit_x) .* fit_x.^2, vis_y);
        % Fit the dimensionless function with linear regression
%         tmp_lm_dl = fitlm(log(fit_x(~tmp_is_ol_Q)./r_cap), tmp_ratio(~tmp_is_ol_Q), 'Intercept', true);
%         tmp_lm_dl_fun = @(x) tmp_lm_dl.Coefficients.Estimate(2) .* log(x ./ r_cap) + tmp_lm_dl.Coefficients.Estimate(1);
%         tmp_lm_dl_cov = corr(tmp_lm_dl_fun(fit_x), tmp_ratio);
%         tmp_opt_sk_cov = corr(tmp_nldl_fun(tmp_est_ab, fit_x) .* fit_x.^2, vis_y);
        
        
%         tmp_fig = figure; 
%         tmp_fig.Position(3:4) = tmp_fig.Position(3:4) .* [2,1];
%         ax_hdl = subplot(1,2,1);
%         histogram2(ax_hdl, fit_x, vis_y, 'DisplayStyle', 'tile');
%         tmp_plot_x = r_cap : ax_hdl.XLim(2);
%         hold(ax_hdl, 'on');
%         tmp_1_plt_hdl = plot(ax_hdl, tmp_plot_x, tmp_nl_fun(tmp_opt_coeff_full, tmp_plot_x));
%         tmp_1_plt_hdl_2 = plot(ax_hdl, tmp_plot_x, tmp_plot_x .^ 2 .* tmp_nldl_fun(tmp_opt_coeff_dl, tmp_plot_x));
%         tmp_1_plt_hdl_3 = plot(ax_hdl, tmp_plot_x, tmp_nl_fun(tmp_est_ab, tmp_plot_x));
%         ax_hdl_2 = subplot(1,2,2);
%         histogram2(ax_hdl_2, fit_x, tmp_ratio, 'DisplayStyle', 'tile');
%         ax_hdl_2.ColorScale = 'log';
%         hold(ax_hdl_2, 'on');
%         tmp_plt_hdl_2 = plot(ax_hdl_2, tmp_plot_x, tmp_nl_fun(tmp_opt_coeff_full, tmp_plot_x) ./ tmp_plot_x .^ 2);
%         tmp_plt_hdl_3 = plot(ax_hdl_2, tmp_plot_x, tmp_nldl_fun(tmp_opt_coeff_dl, tmp_plot_x));
%         tmp_plt_hdl_4 = plot(ax_hdl_2, tmp_plot_x, tmp_nldl_fun(tmp_est_ab, tmp_plot_x));
%         tmp_ymed_in_x_hdl = plot(ax_hdl_2, tmp_ratio_in_dmax_bin.x_bin_val, tmp_ratio_in_dmax_bin.y_median, 'LineWidth', 2);
%         tmp_ymmean_in_x_hdl = plot(ax_hdl_2, tmp_ratio_in_dmax_bin.x_bin_val, tmp_ratio_in_dmax_bin.y_mean, 'LineWidth', 2);
%         leg_hdl = legend(ax_hdl_2, [tmp_plt_hdl_2, tmp_plt_hdl_3, tmp_plt_hdl_4, tmp_ymed_in_x_hdl, tmp_ymmean_in_x_hdl], ...
%             'u_{lm} wrt d_{ulm}', 'u_{lm}/d_{ulm}^2 wrt d_{ulm}', 'Scaled Krogh', 'Median in bin', 'Mean in bin');
%         ax_hdl_2.YLim(1) = -1;
        %% Try a better fit for the deviation from quadratic decay
%         tmp_x = fit_x.^2;
%         tmp_y = vis_y;
%         tmp_ratio = tmp_y ./ tmp_x;
% %         tmp_quadric_fit = fitlm(tmp_x, tmp_y, 'Intercept', false);
% %         tmp_fit_fun = @(x) (- krogh_coeff) * (log(x / r_cap) - 1/2 + r_cap^2 ./ (2 * x .^2));
%         tmp_fit_fun = @(x) log(x / r_cap);
% %         tmp_fit_fun = @(x) (- krogh_coeff) * (log(x / r_cap) - 0.5);
%         tmp_model_coeff = tmp_fit_fun(vis_x);
%         tmp_fit_2 = fitlm(tmp_model_coeff, tmp_ratio);
%         tmp_fit_str = sprintf('a: %.3f\nb: %.3f\nR^2-Adj: %.3f', ...
%             tmp_fit_2.Coefficients.Estimate(2), tmp_fit_2.Coefficients.Estimate(1),...
%             tmp_fit_2.Rsquared.Adjusted);
%         %         tmp_fit_2 = fitlm(tmp_model_coeff, tmp_ratio, 'Intercept', false);
% %         tmp_fit_str = sprintf('\\lambda: %.3f\nR^2-Adj: %.3f', tmp_fit_2.Coefficients.Estimate(1), tmp_fit_2.Rsquared.Adjusted);
%         
%         tmp_fig = figure;
%         tmp_ax = axes(tmp_fig);
%         histogram2(tmp_ax, tmp_model_coeff, tmp_ratio, 'DisplayStyle', 'tile');        
%         hold(tmp_ax, 'on');
%         tmp_plot_x = r_cap : tmp_ax.XLim(2);
%         tmp_plt = plot(tmp_ax, tmp_plot_x, tmp_fit_2.Coefficients.Estimate(2) .* tmp_plot_x ...
%             + tmp_fit_2.Coefficients.Estimate(1));        
% %         tmp_plt = plot(tmp_ax, tmp_plot_x, tmp_fit_2.Coefficients.Estimate(1) .* tmp_fit_fun(tmp_plot_x));
% %         tmp_plt = plot(tmp_ax, tmp_plot_x, tmp_fit_2.Coefficients.Estimate(2) .* tmp_fit_fun(tmp_plot_x) + tmp_fit_2.Coefficients.Estimate(1));
%         legend(tmp_ax, tmp_plt, tmp_fit_str);
% %         tmp_ax.XLabel.String = 'r_{ulm} (\mum)';
% %         tmp_ax.YLabel.String = 'u_{lm}/r_{ulm}^2';
%         tmp_ax.YLim(1) = -1;
%         tmp_ax.FontSize = 14;
%         tmp_cbar = colorbar(tmp_ax);
%         tmp_cbar.Label.String = 'Number of datapoints';
%         tmp_fig_fp = fullfile(save_im_folder, tmp_folder_name, sprintf('%s_%s_%s_corr_coeff_min_c2v_vr_%.2f_wd_sz_%d_um_dimless_part.png', ...
%             dataset_name, merge_stack_name, tmp_folder_name, ...
%             min_cap2vsl_vol_fraction, tmp_wd_sz));
%                 fun_print_image_in_several_formats(tmp_fig, tmp_fig_fp);
%         delete(tmp_fig);
        %% Output figures
        fig_hdl = figure('Visible', 'on');
        ax_hdl = axes(fig_hdl);
        histogram2(ax_hdl, vis_x, vis_y, 'DisplayStyle', 'tile');
        ax_hdl.XLim(1) = 0;
        ax_hdl.YLim(2) = 0;
        hold(ax_hdl, 'on');
        plot_x = 0 : ax_hdl.XLim(2);
        plot_fit_y = - krogh_coeff * ((plot_x + r_cap) .^ 2 .* (log((plot_x + r_cap) / r_cap) - 1/2) + r_cap^2/2);
        line_hdl = plot(ax_hdl, plot_x, plot_fit_y .* linear_fit_hdl.Coefficients.Estimate(1), ...
            'LineWidth', 2);
        ax_hdl.ColorScale = 'log';
        cbar_hdl = colorbar(ax_hdl);
        cbar_hdl.Limits(1) = 1;
        cbar_hdl.Label.String = 'Number of cubes';
        ax_hdl.FontSize = 14;
        ax_hdl.XLabel.String = tmp_x_label;
        ax_hdl.YLabel.String = tmp_y_label;
        ax_hdl.Title.String = sprintf('Window size %d \\mum', tmp_wd_sz);
        
        leg_str = sprintf('Scaled Krogh model\n\\lambda: %.4f \\pm %.4f \nR^2-Adjusted: %.3f\nData size: %d (%.3f)', ...
            linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1), ...
            linear_fit_hdl.Rsquared.Adjusted, ...
            tmp_fit_info_str.num_data_point, tmp_fit_info_str.fraction_of_valid_cube);
        leg_hdl = legend(ax_hdl, line_hdl, leg_str, 'Location', 'southwest');
        
        tmp_fp = fullfile(save_im_folder, tmp_folder_name, sprintf('%s_%s_%s_corr_coeff_min_c2v_vr_%.2f_wd_sz_%d_um.png', ...
            dataset_name, merge_stack_name, tmp_folder_name, ...
            min_cap2vsl_vol_fraction, tmp_wd_sz));
        %%
        fun_print_image_in_several_formats(fig_hdl, tmp_fp);
        delete(fig_hdl);
    end
    scaled_krogh_str.(tmp_folder_name) = cat(1, tmp_fit_info_cell{:});
end
%%
scaled_krogh_str.r_cap = r_cap;
scaled_krogh_str.krogh_coeff = krogh_coeff;
pO2_SA.scaled_krogh_coeff = scaled_krogh_str;
DataManager.write_data(pO2_SA.filepath, pO2_SA);
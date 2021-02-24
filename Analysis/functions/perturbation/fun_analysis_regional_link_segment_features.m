function region_dtpt_str = fun_analysis_regional_link_segment_features(link_features_T, info_str)

persistent DataManager;
if isempty(DataManager)
    DataManager = FileManager;
end
%% Parse input structure
dataset_name = info_str.dataset_name;
stack  = info_str.stack;
vis_group_name = info_str.vis_group_name;
tmp_region_name = info_str.structure_name;
visible_Q = info_str.visibleQ;
subfolder_name = info_str.vis_subgroup_name;
x_label_string = info_str.x_label_string;

if ~isempty(info_str.cc_name)
    vis_folder_name = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
        stack), vis_group_name, strrep(tmp_region_name, ' ', '_'), info_str.cc_name, subfolder_name);
    tmp_region_name = sprintf('%s_%s', tmp_region_name, info_str.cc_name);
else
    vis_folder_name = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
        stack), vis_group_name, strrep(tmp_region_name, ' ', '_'), subfolder_name);
end
%% Estimation
region_dtpt_str = info_str;
region_dtpt_str.filepath = fullfile(vis_folder_name, ...
    sprintf('%s_%s_%s_%s_region_stat_data.mat', dataset_name, stack, ...
    vis_group_name, tmp_region_name));

regional_dt_max_med = median(link_features_T.nearest_tissue_dt_max, 'omitnan');
regional_dt_mean_med = median(link_features_T.nearest_tissue_dt_mean, 'omitnan');

regional_mean_straightness = mean(link_features_T.straightness, 'omitnan');
regional_est_termination_length = 2 * regional_dt_max_med / regional_mean_straightness;

region_dtpt_str.nearest_dt_max_median = regional_dt_max_med;
region_dtpt_str.nearest_dt_mean_median = regional_dt_mean_med;
region_dtpt_str.straightness_mean = regional_mean_straightness;
region_dtpt_str.est_tlength = regional_est_termination_length;
%% Add derivative fields
link_features_T.vol_r_in_bch_od_le_5 = 1 - link_features_T.vol_r_in_bch_od_gt_5;
link_features_T.diff_p_up_2_p_dt_mean = 1 - link_features_T.up_ts_2_p_ts_dt_mean;
%% Quantify the dependence of vessel features on vessel length
bin_feature_info = cell(2, 0);
bin_feature_info(:, end+1) = {'dist_rm_lk_2_nlm', '$dist(lk, nlm)\;(\mu m)$'};
bin_feature_info(:, end+1) = {'dist_rm_lk_2_pt_max_dt', '$dist(lk, d^{(p)}_{max})\;(\mu m)$'};
bin_feature_info(:, end+1) = {'lk_dt_af_rm_max', '$max{d^{(p)}(\mathbf{x})|\mathbf{x}\in lk\} \;(\mu m)$'};
bin_feature_info(:, end+1) = {'lk_dt_af_rm_mean', '$\bar{d}^{(p)}_{lk} \;(\mu m)$'};
bin_feature_info(:, end+1) = {'lk_dt_af_rm_median', '$dist(lk, nlm)\;(\mu m)$'};
bin_feature_info(:, end+1) = {'nb_lk_vol_r_max', 'RVF_{max}'};
bin_feature_info(:, end+1) = {'nb_lk_vol_r_max_bch_od', 'BO(RVF_{max})'};
bin_feature_info(:, end+1) = {'nearest_tissue_dt_max', '$d_{max}$'};
bin_feature_info(:, end+1) = {'nearest_tissue_dt_mean', '$\bar{d}_{max}$'};
bin_feature_info(:, end+1) = {'nlm_v_af_rm', '$d^{(p)}_{nlm}$'};
bin_feature_info(:, end+1) = {'num_nb_lk', 'N_{nlk}'};
bin_feature_info(:, end+1) = {'pt_vol_dt_max', '$d^{(p)}_{max}\;(\mu m)$'};
bin_feature_info(:, end+1) = {'pt_vol_dt_mean', '$\bar{d}^{(p)}_{max}$'};
bin_feature_info(:, end+1) = {'shortest_loop_geodesic_length', 'N_{sl}'};
bin_feature_info(:, end+1) = {'straightness', 'S'};
bin_feature_info(:, end+1) = {'up_ts_2_p_ts_dt_mean', '$\langle d/d^{(p)}\rangle$'};
bin_feature_info(:, end+1) = {'diff_p_up_2_p_dt_mean', '$\langle 1 - d/d^{(p)}\rangle$'};
bin_feature_info(:, end+1) = {'vol_r_in_bch_od_1', 'RVF(BO = 1)'};
bin_feature_info(:, end+1) = {'vol_r_in_bch_od_gt_5', 'RVF(BO > 5)'};
bin_feature_info(:, end+1) = {'vol_r_in_bch_od_le_5', '$RVF(BO \le 5)$'};
bin_feature_info(:, end+1) = {'diff_dp_d_mean', '$\langle d^{(p)} - d\rangle \;(\mu m)$'};
num_bin_feature = size(bin_feature_info, 2);

plot_x_data = link_features_T.length;
plot_x_edge = [0 : 2.5 : 25, 30 : 5 : 50, 60 : 10 : 200];

for iter_feature = 1 : num_bin_feature
    tmp_x_data = plot_x_data;
    tmp_y_field_name = bin_feature_info{1, iter_feature};
    tmp_y_data = fun_getfield(link_features_T, tmp_y_field_name);
    y_label_string = bin_feature_info{2, iter_feature};
    
    [tmp_y_binned_str, ~, ~] = fun_analysis_get_y_stat_in_x_bin(tmp_x_data, ...
        tmp_y_data, plot_x_edge);
    tmp_y_binned_str.x_label_string = x_label_string;
    tmp_y_binned_str.y_label_string = y_label_string;
    % Record for cross-region comparison
    region_dtpt_str.(tmp_y_field_name) = tmp_y_binned_str;
end

%% Cumulate length distribution
tmp_length = sort(link_features_T.length, 'ascend');
total_cap_length = sum(tmp_length);
tmp_cum_length = cumsum(tmp_length);
[tmp_unique_length, tmp_idx, ~] = unique(tmp_length, 'stable');
tmp_idx = [tmp_idx(2:end) - 1; numel(tmp_length)];
tmp_cum_length = tmp_cum_length(tmp_idx);
cum_cap_length_itp = griddedInterpolant(tmp_unique_length, tmp_cum_length);
region_dtpt_str.cum_cap_len_frc_itp = cum_cap_length_itp;
region_dtpt_str.total_cap_length_um = total_cap_length;
region_dtpt_str.total_cap_length_in_cc_um = sum(link_features_T.length .* ...
    link_features_T.in_cc_ratio, 'omitnan');
% Find the length at which Median N_RVF_max turns to 2
region_dtpt_str.min_cap_length_NRVFmax_2 = region_dtpt_str.nb_lk_vol_r_max_bch_od.x_bin_val(...
    find(region_dtpt_str.nb_lk_vol_r_max_bch_od.y_median == 2, 1, 'first'));
%% Maximum distance increment vs original maximum distance
tmp_y_data = link_features_T.pt_vol_dt_max - link_features_T.nearest_tissue_dt_max;
[tmp_y_binned_str, ~, ~] = ...
    fun_analysis_get_y_stat_in_x_bin(link_features_T.nearest_tissue_dt_max, ...
    tmp_y_data, 0 : 2.5 : ceil(prctile(link_features_T.nearest_tissue_dt_max, 99.5)));
tmp_y_binned_str.x_label_string = '$d_{max}\;(\mu m)$';
tmp_y_binned_str.y_label_string = '$(d^{(p)}_{max} - d_{max})\;(\mu m)$';
region_dtpt_str.delta_dmax = tmp_y_binned_str;
%% Maximum DT passed by the removed link vs dist to nearest nlm
tmp_x_data = link_features_T.lk_dt_af_rm_max;
tmp_y_data = link_features_T.dist_rm_lk_2_nlm;
[tmp_y_binned_str, ~, ~] = ...
    fun_analysis_get_y_stat_in_x_bin(tmp_x_data, ...
    tmp_y_data, 0 : 2.5 : ceil(prctile(tmp_x_data, 99.5)));
tmp_y_binned_str.x_label_string = '$max\{d^{(p)}(\mathbf{x})|\mathbf{x}\in lk\}\;(\mu m)$';
tmp_y_binned_str.y_label_string = '$dist(lk, d^{(p)}_{nlm})(\mu m)$';
region_dtpt_str.maxDTp_vs_dist2nml = tmp_y_binned_str;

if ~isfolder(vis_folder_name)
    mkdir(vis_folder_name);
end
save(region_dtpt_str.filepath, '-struct', 'region_dtpt_str');
%% Visualization
if info_str.visQ
    %% 2D histogram - unperturbed space, segment length dependence of d_avg and d_max
    pt_prop_info = cell(3, 0);
    % New DT maximum in the perturbed space
    pt_prop_info(:, end+1) = {'nearest_tissue_dt_mean', 0 : 2.5 : ceil(regional_dt_mean_med * 2 / 5) * 5, '$d_{avg}\;(\mu m)$'};
    % Distance between the removed link and the new DT maximum
    pt_prop_info(:, end+1) = {'nearest_tissue_dt_max', 0 : 5 : ceil(regional_dt_max_med * 2 / 5) * 5, '$d_{max}\;(\mu m)$'};
    num_subplot = size(pt_prop_info, 2);
    
    plot_x_data = link_features_T.length;
    plot_x_edge = [1 : 1 : 9, 10 .^ linspace(log10(10), log10(200), 30)];
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1.5, 1.5];
    for iter_plot = 1 : num_subplot
        tmp_x_data = plot_x_data;
        
        tmp_y_field_name = pt_prop_info{1, iter_plot};
        tmp_y_data = fun_getfield(link_features_T, tmp_y_field_name);
        tmp_y_edge = pt_prop_info{2, iter_plot};
        y_label_string = pt_prop_info{3, iter_plot};
        
        [tmp_y_binned_str, tmp_x_data, tmp_y_data] = ...
            fun_analysis_get_y_stat_in_x_bin(tmp_x_data, ...
            tmp_y_data, plot_x_edge);
        
        tmp_log_fit_Q = tmp_y_binned_str.x_bin_val < 2 * regional_dt_max_med;
        linear_fit_hdl = fitlm(log(tmp_y_binned_str.x_bin_val(tmp_log_fit_Q)), tmp_y_binned_str.y_median(tmp_log_fit_Q));
        % Visualization
        ax_hdl_1 = subplot(num_subplot, 1, iter_plot);
        histogram2(ax_hdl_1, tmp_x_data, tmp_y_data, plot_x_edge, tmp_y_edge, 'DisplayStyle', 'tile');
        if any(strfind(y_label_string, '$'))
            ax_hdl_1.YLabel.Interpreter = 'latex';
        end
        ax_hdl_1.YLabel.String = y_label_string;
        
        tmp_cbar = colorbar;
        tmp_cbar.Label.String = 'Number of segments';
        tmp_cbar.Limits(1) = 1;
        tmp_cbar.Ticks = 10 .^ [1 : round(log10(tmp_cbar.Limits(2)))];
        
        ax_hdl_1.ColorScale = 'log';
        
        hold(ax_hdl_1, 'on');
        plot_med_hdl = plot(ax_hdl_1, tmp_y_binned_str.x_bin_val, tmp_y_binned_str.y_median, 'Color', 'k', 'LineWidth', 4);
        
        lfit_hdl = plot(ax_hdl_1, tmp_y_binned_str.x_bin_val, ...
            log(tmp_y_binned_str.x_bin_val) * linear_fit_hdl.Coefficients.Estimate(2) + ...
            linear_fit_hdl.Coefficients.Estimate(1), 'Color', 'r', 'LineStyle', ':', 'LineWidth', 4);
        
        est_terminate_l_hdl = line(ax_hdl_1, [regional_est_termination_length, regional_est_termination_length], ...
            [0, tmp_y_edge(end)], 'Color', 'g', 'LineWidth', 4, 'LineStyle', '-.');
        
        legend(ax_hdl_1, [plot_med_hdl, lfit_hdl, est_terminate_l_hdl], {'Median of y binned by x', ...
            sprintf('$%.2f\\;\\log(l) + %.2f$\n$R^2$-Adjusted: %.3f', linear_fit_hdl.Coefficients.Estimate(2), ...
            linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Rsquared.Adjusted), '2 Median $d_{max}/\bar{S_c}$'}, 'Location', 'southeast', 'Interpreter', 'latex');
        box(ax_hdl_1, 'off');
        ax_hdl_1.XLabel.String = x_label_string;
        if iter_plot ~= num_subplot
            ax_hdl_1.XAxis.Visible = 'off';
        end
        ax_hdl_1.FontSize = 14;
    end
    
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_up_dt_avg_n_max_vs_%s.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_'), subfolder_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% 2D Histograms - Perturbed space properties vs capillary length
    pt_prop_info = cell(3, 0);
    % New DT maximum in the perturbed space
    pt_prop_info(:, end+1) = {'pt_vol_dt_max', 0 : 5 : ceil(regional_dt_max_med * 2 / 5) * 5, '$d^{(p)}_{max}\; (\mu m)$'};
    % Distance between the removed link and the new DT maximum
    pt_prop_info(:, end+1) = {'dist_rm_lk_2_pt_max_dt', 0 : 5 : 40, '$dist(lk, d^{(p)}_{max})\;(\mu m)$'};
    % Local DT maximum (with window size 10) - not necessarily exist.
    pt_prop_info(:, end+1) = {'nlm_v_af_rm', 0 : 5 : ceil(regional_dt_max_med * 2 / 5) * 5, '$d^{(p)}_{l,\;max}\; (\mu m)$'};
    pt_prop_info(:, end+1) = {'dist_rm_lk_2_nlm', 0 : 5 : 40, '$dist(lk, nlm)\;(\mu m)$'};
    
    num_subplot = size(pt_prop_info, 2);
    plot_x_data = link_features_T.length;
    plot_x_edge = [0 : 2.5 : 25, 30 : 5 : 50, 60 : 10 : 200];
    
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 2] .* 1.5;
    for iter_plot = 1 : num_subplot
        tmp_x_data = plot_x_data;
        tmp_y_field_name = pt_prop_info{1, iter_plot};
        tmp_y_data = fun_getfield(link_features_T, tmp_y_field_name);
        tmp_y_edge = pt_prop_info{2, iter_plot};
        y_label_string = pt_prop_info{3, iter_plot};
        
        [tmp_y_binned_str, tmp_x_data, tmp_y_data] = ...
            fun_analysis_get_y_stat_in_x_bin(tmp_x_data, ...
            tmp_y_data, plot_x_edge);
        
        %% Visualization
        ax_hdl = subplot(num_subplot, 1, iter_plot);
        histogram2(ax_hdl, tmp_x_data, tmp_y_data, plot_x_edge, tmp_y_edge, 'DisplayStyle', 'tile');
        hold(ax_hdl, 'on');
        plot_med_hdl = plot(ax_hdl, tmp_y_binned_str.x_bin_val, tmp_y_binned_str.y_median, 'k', 'LineWidth', 3);
        est_terminate_l_hdl = line(ax_hdl, [regional_est_termination_length, regional_est_termination_length], ...
            [0, tmp_y_edge(end)], 'Color', 'g', 'LineWidth', 4, 'LineStyle', '-.');
        if any(strfind(y_label_string, '$'))
            ax_hdl.YLabel.Interpreter = 'latex';
        end
        ax_hdl.YLabel.String = y_label_string;
        % Colorbar
        tmp_cbar = colorbar;
        tmp_cbar.Limits(1) = 1;
        tmp_cbar.Label.String = 'Counts';
        tmp_cbar.Ticks = 10 .^ [1 : round(log10(tmp_cbar.Limits(2)))];
        ax_hdl.ColorScale = 'log';
        ax_hdl.FontSize = 14;
        % Bin the y data by x data and plot the box plot
        legend(ax_hdl,  [plot_med_hdl, est_terminate_l_hdl], {'Median of y binned by x', '2 Median $d_{max} /\bar{S_c}$'}, 'Interpreter', 'latex');
        box(ax_hdl, 'off');
        % Add x axis
        if iter_plot == num_subplot
            ax_hdl.XLabel.String = x_label_string;
        else
            ax_hdl.XAxis.Visible = 'off';
        end
    end
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_d_max_n_dist_2_max_vs_%s.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_'), subfolder_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% 2D Histograms - Perturbed space properties vs capillary length
    pt_prop_info = cell(3, 0);%
    pt_prop_info(:, end+1) = {'up_ts_2_p_ts_dt_mean', 0.5 : 0.05 : 1, '$\langle d/d^{(p)}\rangle$'};
    pt_prop_info(:, end+1) = {'nb_lk_vol_r_max_bch_od', 0.5 : 1 : 9.5, '$BO_{RVF_{max}}$'};
    pt_prop_info(:, end+1) = {'vol_r_in_bch_od_1', 0 : 0.1 : 1, '$RVF(N = 1)$'};
    pt_prop_info(:, end+1) = {'vol_r_in_bch_od_gt_5', 0 : 0.1 : 1, '$RVF(N > 5)$'};
    
    % pt_prop_info(:, end+1) = {'straightness', 0.5 : 0.1 : 1, '$S_c$'};
    num_subplot = size(pt_prop_info, 2);
    plot_x_data = link_features_T.length;
    plot_x_edge = [0 : 2.5 : 25, 30 : 5 : 50, 60 : 10 : 200];
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 2] .* 1.5;
    for iter_plot = 1 : num_subplot
        ax_hdl = subplot(num_subplot, 1, iter_plot);
        %     ax_hdl = nexttile;
        tmp_x_data = plot_x_data;
        tmp_y_field_name = pt_prop_info{1, iter_plot};
        tmp_y_data = fun_getfield(link_features_T, tmp_y_field_name);
        tmp_y_edge = pt_prop_info{2, iter_plot};
        y_label_string = pt_prop_info{3, iter_plot};
        
        [tmp_y_binned_str, tmp_x_data, tmp_y_data] = ...
            fun_analysis_get_y_stat_in_x_bin(tmp_x_data, ...
            tmp_y_data, plot_x_edge);
       
        %% Visualization
        histogram2(ax_hdl, tmp_x_data, tmp_y_data, plot_x_edge, tmp_y_edge, 'DisplayStyle', 'tile');
        hold(ax_hdl, 'on');
        plot_med_hdl = plot(ax_hdl, tmp_y_binned_str.x_bin_val, tmp_y_binned_str.y_median, 'k', 'LineWidth', 3);
        est_terminate_l_hdl = line(ax_hdl, [regional_est_termination_length, regional_est_termination_length], ...
            [tmp_y_edge(1), tmp_y_edge(end)], 'Color', 'g', 'LineWidth', 4, 'LineStyle', '-.');
        if any(strfind(y_label_string, '$'))
            ax_hdl.YLabel.Interpreter = 'latex';
        end
        ax_hdl.YLabel.String = y_label_string;
        % Colorbar
        tmp_cbar = colorbar;
        tmp_cbar.Limits(1) = 1;
        tmp_cbar.Label.String = 'Counts';
        tmp_cbar.Ticks = 10 .^ [1 : round(log10(tmp_cbar.Limits(2)))];
        ax_hdl.ColorScale = 'log';
        ax_hdl.FontSize = 14;
        % Bin the y data by x data and plot the box plot
        legend(ax_hdl,  [plot_med_hdl, est_terminate_l_hdl], {'Median of y binned by x', '2 Median $d_{max} /\bar{S_c}$'}, 'Interpreter', 'latex');
        box(ax_hdl, 'off');
        % Add x axis
        if iter_plot == num_subplot
            ax_hdl.XLabel.String = x_label_string;
        else
            ax_hdl.XAxis.Visible = 'off';
        end
    end
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_DT_reduction_n_RVF_vs_%s.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_'), subfolder_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% 2D Histogram - Local maximum DT before vs after - useful.
    tmp_x_data = link_features_T.nearest_tissue_dt_max;
    tmp_y_data = link_features_T.pt_vol_dt_max;
    
    tmp_is_valid_Q = ~isnan(tmp_x_data) & ~isnan(tmp_y_data);
    tmp_x_data = tmp_x_data(tmp_is_valid_Q);
    tmp_y_data = tmp_y_data(tmp_is_valid_Q);
    
    plot_x_edge = 0 : 2.5 : 45;
    tmp_y_edge = 0 : 2.5 : 45;
    
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 1.5];
    ax_hdl = subplot(4, 3, 1:9);
    histogram2(ax_hdl, tmp_x_data, tmp_y_data, plot_x_edge, tmp_y_edge, 'DisplayStyle', 'tile');
    ax_hdl.XAxis.Visible = 'off';
    ax_hdl.YLabel.Interpreter = 'latex';
    ax_hdl.YLabel.String = '$d^{(p)}_{max}\;(\mu m)$';
    tmp_cbar = colorbar;
    tmp_cbar.Label.String = 'Counts';
    tmp_cbar.Location = 'northoutside';
    ax_hdl.ColorScale = 'linear';
    ax_hdl.FontSize = 14;
    
    tmp_list_ind = fun_bin_data_to_idx_list_by_edges(tmp_x_data, plot_x_edge, true);
    tmp_y_data_binned = fun_bin_data_to_cells_by_ind(tmp_y_data, tmp_list_ind);
    tmp_x_data_binned = fun_bin_data_to_cells_by_ind(tmp_x_data, tmp_list_ind);
    tmp_y_data_avg = cellfun(@mean, tmp_y_data_binned);
    tmp_x_data_avg = cellfun(@mean, tmp_x_data_binned);
    hold(ax_hdl, 'on');
    box(ax_hdl, 'off');

    plot_med_hdl = plot(ax_hdl, tmp_x_data_avg, tmp_y_data_avg, 'LineWidth', 3, 'Color', 'r');
    legend(plot_med_hdl, 'Mean of y binned by x', 'Location', 'northwest');
    %
    ax_hdl_2 = subplot(4,3, 10:12);
    delta_dmax = tmp_y_data_avg - tmp_x_data_avg;
    plot(ax_hdl_2, tmp_x_data_avg, delta_dmax, 'LineWidth', 2, 'Color', 'r');
    ax_hdl_2.XLim = ax_hdl.XLim;
    ax_hdl_2.YLim(2) = ceil(ax_hdl_2.YLim(2));
    box(ax_hdl_2, 'off');
    grid(ax_hdl_2, 'on');
    ax_hdl_2.XLabel.Interpreter = 'latex';
    ax_hdl_2.XLabel.String = '$d_{max}\;(\mu m)$';
    ax_hdl_2.YLabel.Interpreter = 'latex';
    ax_hdl_2.YLabel.String = '$d^{(p)}_{max} - d_{max} \;(\mu m)$';
    ax_hdl_2.FontSize = 14;
    
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_d_max_after_vs_before.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_')));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% Redistributed fraction vs order CDF for vessels group by length
    length_group = [0, 20; 20, 40; 40, 60; 60, 80; 80, 100; 100, inf; 0, inf];
    num_length_group = size(length_group, 1);
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3 3] .* 1.2;
    for iter_group = 1 : num_length_group
        tmp_length_min = length_group(iter_group, 1);
        tmp_length_max = length_group(iter_group, 2);
        tmp_valid_Q =  (link_features_T.length >= tmp_length_min) & ...
            (link_features_T.length <= tmp_length_max);
        vol_r_binned = link_features_T.vol_r_in_bch_od(tmp_valid_Q, :);
        
        ax_hdl = subplot(3, 3, iter_group);
        boxplot(ax_hdl, vol_r_binned, 'Labels', {'1', '2', '3', '4', '5', '> 5'});
        ax_hdl.XLabel.String = 'Braching order';
        ax_hdl.YLabel.String = 'Volume fraction';
        ax_hdl.FontSize = 14;
        ax_hdl.Title.String = sprintf('Capillary length range [%d, %d] \\mum', ...
            tmp_length_min, tmp_length_max);
    end
    
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_redistributed_volume_fraction_vs_branching_order.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_')));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% 2D Histogram - Maximum DT passed by the removed link vs distance to nearest DT local maximum
    tmp_x_data = link_features_T.lk_dt_af_rm_max;
    tmp_y_data = link_features_T.dist_rm_lk_2_nlm;
    
    tmp_is_valid_Q = ~isnan(tmp_x_data) & ~isnan(tmp_y_data);
    tmp_x_data = tmp_x_data(tmp_is_valid_Q);
    tmp_y_data = tmp_y_data(tmp_is_valid_Q);
    
    plot_x_edge = 0 : 5 : 45;
    tmp_y_edge = 0 : 5 : 45;
    
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    ax_hdl = axes(fig_hdl);
    histogram2(ax_hdl, tmp_x_data, tmp_y_data, plot_x_edge, tmp_y_edge, 'DisplayStyle', 'tile');
    ax_hdl.YLabel.Interpreter = 'latex';
    ax_hdl.XLabel.Interpreter = 'latex';
    ax_hdl.YLabel.String = '$dist(lk, d^{(p)}_{nlm})(\mu m)$';
    ax_hdl.XLabel.String = '$max\{d^{(p)}(\mathbf{x})|\mathbf{x}\in lk\}\;(\mu m)$';
    tmp_cbar = colorbar;
    tmp_cbar.Label.String = 'Counts';
    tmp_cbar.Limits(1) = 1;
    ax_hdl.ColorScale = 'log';
    ax_hdl.FontSize = 14;
    
    % Bin the y data by x data and plot the box plot
    tmp_list_ind = fun_bin_data_to_idx_list_by_edges(tmp_x_data, plot_x_edge, true);
    tmp_y_data_binned = fun_bin_data_to_cells_by_ind(tmp_y_data, tmp_list_ind);
    tmp_y_data_med = cellfun(@median, tmp_y_data_binned);
    tmp_x_edge_avg = movmean(plot_x_edge, 2, 'Endpoints', 'discard');
    hold(ax_hdl, 'on');
    plot_med_hdl = plot(ax_hdl, tmp_x_edge_avg, tmp_y_data_med, 'LineWidth', 3, 'Color', 'r');
    
    % Average DT Max in this block
    leg_hdl = legend(plot_med_hdl, 'Median of y binned by x');
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_dist2nlm_vs_max_dt_passed_by_lk.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_')));
    
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    %% length cumulated distribution
    tmp_plot_length = 0.5 : 5: 300;
    tmp_plot_cdf = cum_cap_length_itp(tmp_plot_length) ./ total_cap_length;
    if visible_Q
        fig_hdl = figure;
    else
        fig_hdl = figure('Visible', 'off');
    end
    ax_hdl = axes(fig_hdl);
    yyaxis(ax_hdl, 'left');
    ldf_hdl = plot(ax_hdl, movmean(tmp_plot_length, 2, 'Endpoints', 'discard'), diff(tmp_plot_cdf), 'LineWidth', 2);
    ax_hdl.YLabel.String = 'Length fraction';
    yyaxis(ax_hdl, 'right');
    cldf_hdl = plot(ax_hdl, tmp_plot_length, tmp_plot_cdf, 'LineWidth', 2);
    ax_hdl.YLabel.String = 'Cumulative length fraction';
    ax_hdl.XScale = 'log';
    ax_hdl.XLim = [2.5, tmp_plot_length(end)];
    ax_hdl.XLabel.String = 'Capillary segment length (\mum)';
    ax_hdl.FontSize = 14;
    legend(ax_hdl, cldf_hdl, sprintf('Total length: %.2f m', total_cap_length /1e6), 'Location', 'northwest');
    ax_hdl.Title.String = strrep(tmp_region_name, '_', ' ');
    fig_fp = fullfile(vis_folder_name, sprintf('%s_%s_%s_%s_length_fraction_n_clf_vs_%s.png', ...
        dataset_name, stack, vis_group_name, strrep(tmp_region_name, ' ', '_'), subfolder_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
end
end
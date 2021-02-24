%% Load data
DataManager = FileManager;
recon_r = 0.25;
% result_fp = DataManager.fp_analysis_data_file('WholeBrain', 'ML_2018_08_15', ...
%     'Lattice_space_filling_link_properties_table.mat');
% lattice_valid_data_cell = DataManager.load_data(result_fp, lattice_valid_data_cell);
stat_fp = DataManager.fp_analysis_data_file('WholeBrain', 'ML_2018_08_15', ...
    sprintf('Lattice_space_filling_stat_table_min_r_%.2f_um.mat', recon_r));
lattice_stat_cell = DataManager.load_data(stat_fp);
%% Preprocessing
has_valid_structure_Q = ~cellfun(@isempty, lattice_stat_cell);
lattice_stat_table = lattice_stat_cell(has_valid_structure_Q);
% selected_data_group = 'dt_lm_link_stat';
selected_data_group = 'all_link_stat';
for iter_cell = 1 : numel(lattice_stat_table)
    tmp_data = lattice_stat_table{iter_cell};
    tmp_str = struct;
    tmp_str.lattice_name = tmp_data.lattice_name;
    tmp_str.set_length_pxl = tmp_data.set_length_pxl;
    tmp_str.system_size = tmp_data.system_size;
    tmp_str = fun_copy_fields_from_str2_to_str1(tmp_str, tmp_data.(selected_data_group));
    lattice_stat_table{iter_cell} = struct2table(tmp_str, 'AsArray', true);
end
tmp_str = [];
lattice_stat_table = cat(1, lattice_stat_table{:});
% table_fp = './Simulation/lattice_space_filling_stat_05_ds_r2_v2.csv';
% writetable(lattice_stat_table, table_fp);
%%
lattice_stat_table.ValidTotalLengthDensity_isqrt = (lattice_stat_table.ValidTotalLengthDensity) .^ (-1/2);
lattice_stat_table.ValidTotalEp2epDistDensity_isqrt = (lattice_stat_table.ValidTotalEp2epDistDensity) .^ (-1/2);
lattice_stat_table.TargetTotalLengthDensity_isqrt = (lattice_stat_table.TargetTotalLengthDensity) .^ (-1/2);
%%
lattice_name_list = {'(10, 3)-a', '(10, 3)-b', '(10, 3)-c', '(8, 3)-a', 'cubic'};
num_lattice = numel(lattice_name_list);
scaling_exp_str.lattice_name = lattice_name_list;
%% Setting plotting parameters
plot_setting = cell(2, 0);
plot_setting(:, end+1) = {'length_mean', 'Average l (\mum)'};
plot_setting(:, end+1) = {'length_median', 'Median l (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_dt_max_mean', 'Average d_{max} (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_dt_max_median', 'Median d_{max} (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_dt_mean_mean', 'Average d_{avg} (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_dt_mean_median', 'Median d_{avg} (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_dt_median_mean', 'Average d_{med} (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_dt_median_median', 'Median d_{med} (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_radius_mean', 'Average d_r (\mum)'};
plot_setting(:, end+1) = {'nearest_tissue_radius_median', 'Median d_r (\mum)'};
plot_setting(:, end+1) = {'tissue_dt_prctile99_mean', 'Average CDF_{d}(99) (\mum)'};
plot_setting(:, end+1) = {'tissue_dt_prctile99_median', 'Median CDF_{d}(99) (\mum)'};
plot_setting(:, end+1) = {'delta_dt_max_mean', 'Average (d_{max, p} - d_{max}) (\mum)'};
plot_setting(:, end+1) = {'delta_dt_max_median', 'Median (d_{max, p} - d_{max}) (\mum)'};

plot_setting(:, end+1) = {'up_ts_2_p_ts_dt_mean_mean', 'Average <d_p/d>'};
plot_setting(:, end+1) = {'up_ts_2_p_ts_dt_mean_median', 'Median <d_p/d>'};

plot_setting(:, end+1) = {'dist_rm_lk_2_pt_max_dt_mean', 'Average dist(lk, d_{max}) (\mum)'};
plot_setting(:, end+1) = {'dist_rm_lk_2_pt_max_dt_median', 'Median dist(lk, d_{max}) (\mum)'};

plot_setting(:, end+1) = {'pt_vol_dt_max_mean', 'Average d_{max, p} (\mum)'};
plot_setting(:, end+1) = {'pt_vol_dt_max_median', 'Median d_{max, p} (\mum)'};

plot_setting(:, end+1) = {'dist_rm_lk_2_nlm_mean', 'Average dist(lk, d_{nlm, p}) (\mum)'};
plot_setting(:, end+1) = {'dist_rm_lk_2_nlm_median', 'Median dist(lk, d_{nlm, p}) (\mum)'};

plot_setting(:, end+1) = {'nb_lk_bch_od_1_vol_r_mean', 'Average RVF(BO = 1)'};
plot_setting(:, end+1) = {'nb_lk_bch_od_1_vol_r_median', 'Median RVF(BO = 1)'};

num_features = size(plot_setting, 2);
%% Compute the scaling relation for each type of lattices
tmp_x_name = 'TargetTotalLengthDensity_isqrt';
tmp_x_label = '\rho_l^{-1/2} (\mum)';
% tmp_x_name = 'ValidTotalEp2epDistDensity_isqrt';
% tmp_x_label = '$\rho_{lep}^{-1/2} (\mu m)$';
%%
for iter_y = 1 : num_features
    tmp_y_name = plot_setting{1, iter_y};
    tmp_y_label = plot_setting{2, iter_y};    
    [tmp_str.slope, tmp_str.intercept, tmp_str.R2A, tmp_str.R2, ...
        tmp_str.slope_se, tmp_str.intercept_se] = deal(nan(num_lattice, 1));
    
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
    ax_hdl = axes(fig_hdl);
    leg_hdl_array = [];
    leg_str_array = {};
    for iter_lattice = 1 : num_lattice
        lattice_name = lattice_name_list{iter_lattice};
        tmp_table_selected_Q = strcmp(lattice_stat_table.lattice_name, lattice_name);
        
        tmp_x = lattice_stat_table.(tmp_x_name);
        tmp_y = lattice_stat_table.(tmp_y_name);
        tmp_x = tmp_x(tmp_table_selected_Q);
        tmp_y = tmp_y(tmp_table_selected_Q);
        %%
        dot_hdl = scatter(ax_hdl, tmp_x, tmp_y);
        hold(ax_hdl, 'on');
        ax_hdl.XLim(1) = 0;
        fit_data_x = linspace(0, max(tmp_x, [], 'all'), 30);
        linear_fit_hdl = fitlm(tmp_x, tmp_y);
        fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
            'LineWidth', 2, 'LineStyle', '-', 'Color', dot_hdl.CData);
        leg_str_array = [leg_str_array, ...
            {sprintf('Lattice: %s\nIntercept: %.4f \\pm %.4f\nSlope: %.4f \\pm %.4f\nR^2-Adjusted: %.4e', ...
            lattice_name,...
            linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
            linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
            linear_fit_hdl.Rsquared.Adjusted)}];
        leg_hdl_array = [leg_hdl_array, fit_plt_hdl];
        %%
        tmp_str.slope(iter_lattice) = linear_fit_hdl.Coefficients.Estimate(2);
        tmp_str.slope_se(iter_lattice) = linear_fit_hdl.Coefficients.SE(2);
        tmp_str.intercept(iter_lattice) = linear_fit_hdl.Coefficients.Estimate(1);
        tmp_str.intercept_se(iter_lattice) = linear_fit_hdl.Coefficients.SE(1);
        tmp_str.R2A(iter_lattice) = linear_fit_hdl.Rsquared.Adjusted;
        tmp_str.R2(iter_lattice) = linear_fit_hdl.Rsquared.Ordinary;
    end
    leg_hdl = legend(leg_hdl_array, leg_str_array, 'Location', 'southeast');
    leg_hdl.FontSize = 12;
    ax_hdl.FontSize = 18;
    ax_hdl.XLabel.String = tmp_x_label;
    ax_hdl.YLabel.String = tmp_y_label;
    %%
    fig_fp = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'all_stack'), ...
        'Space_filling_lattice', selected_data_group, tmp_x_name, sprintf('Space_filling_lattices_%s_vs_%s_05_ds_r2.png', tmp_y_name, tmp_x_name));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
    
    scaling_exp_str.(tmp_y_name) = tmp_str;
end
scaling_exp_str.y_name = plot_setting(1, :)';
scaling_exp_str.y_label_name = plot_setting(2, :)';
save('./Simulation/space_filling_lattice/lattice_space_filling_scaling.mat', '-struct', 'scaling_exp_str');
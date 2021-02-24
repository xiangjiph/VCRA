%%
plot_marker_style = {'none', '+', 'o', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
plot_line_style = {'-', '--', ':', '-.'};
num_marker_style = numel(plot_marker_style);
num_line_style = numel(plot_line_style);
[~, tmp_line_style_ind, tmp_marker_style_ind] = ndgrid(1 : 6, 1 : num_line_style, 1 : num_marker_style);
%% Compute the bond percolation threshold for vascular networks
fig_group_name = 'regional_network_percolation';
cc_label = {'Left', 'Right'};
percolation_region_data = cell(2, num_stack, num_region);
for iter_region = 1 : num_region
    for iter_stack = 1 : num_stack
        test_region_data = region_data{iter_region, iter_stack};
        for iter_cc = 1 : 2
            tmp_cc_name = sprintf('cc_%d', iter_cc);
            test_graph = test_region_data.(tmp_cc_name).vessel_graph;
            graph_ud_str = fun_analysis_get_connectivity_graph(test_graph);
            graph_ud = graph_ud_str.graph_w;
            %%
            tic
            percolation_p =  sort(unique([0.05 : 0.05 : 0.95, 0.45 : 0.0005 : 0.65]), 'ascend');
            percolation_str = fun_analysis_percolation_transition_by_bond_removal(graph_ud, percolation_p,  ...
                500, true);
            toc
            %% Visualization
            fig_hdl = figure;
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1];
            ax_hdl = subplot(1,3,1);
            plot(ax_hdl, percolation_str.bond_occupancy_p , percolation_str.avg.largest_cc_fraction, 'LineWidth', 2);
            ax_hdl.YLabel.String = 'Largest cluster fraction (s)';
            ax_hdl.XLabel.String = 'Occupation probability (p)';
            ax_hdl.XLim = [0,1];
            ax_hdl.YLim = [0,1];
%             yyaxis(ax_hdl, 'right');
%             plot(ax_hdl, percolation_str.bond_occupancy_p , percolation_str.avg.largest_cc_num_site);
%             ax_hdl.YLabel.String = 'Number of site in the largest cluster';
%             ax_hdl.XLabel.String = 'Occupation probability (p)';
%             ax_hdl.XLim = [0,1];
%             ax_hdl.YLim = [0, percolation_str.num_sites];
            % Find the value of remove fraction s.t. largest cc ratio = 0.5 by linear
            % interpolation
            tmp_pc = fun_analysis_percolation_compute_pc(percolation_str.bond_occupancy_p, ...
                percolation_str.avg.largest_cc_fraction);
            [int_x, tmp_ind] = sort(percolation_str.avg.largest_cc_fraction, 'ascend');
            int_y = percolation_str.bond_occupancy_p(tmp_ind);
            cc_ratio_2_p = griddedInterpolant(int_x, int_y);
            
            ratio_half_th = cc_ratio_2_p(0.5);
            ratio_095_th = cc_ratio_2_p(0.95);
            ratio_005_th = cc_ratio_2_p(0.05);
            ratio_001_th = cc_ratio_2_p(0.01);
            legend(ax_hdl, sprintf('Number of sites = %d\nNumber of bonds = %d\np_c = %.4f\np(s = 0.05) = %.4f\np(s = 0.50) = %.4f\np(s = 0.95) = %.4f', ...
                percolation_str.num_sites, percolation_str.num_bonds, tmp_pc, ratio_005_th, ratio_half_th, ratio_095_th), 'Location', 'northwest');
            
            ax_hdl_2 = subplot(1,3,2);
            s_range = [0.05, 0.3];
            p_range = cc_ratio_2_p(s_range);
            
            loglog_fit_Q = (percolation_p <= p_range(2) & percolation_p >= p_range(1));
            fit_x = percolation_p(loglog_fit_Q) - p_range(1);
            fit_y = log10(percolation_str.avg.largest_cc_fraction(loglog_fit_Q));
            scatter(ax_hdl_2, fit_x, fit_y);
            lin_log_fit = fitlm(fit_x, fit_y);
            hold(ax_hdl_2, 'on');
            plt_hdl = plot(ax_hdl_2, fit_x, fit_x * lin_log_fit.Coefficients.Estimate(2) + lin_log_fit.Coefficients.Estimate(1), 'LineWidth', 2);
            % ax_hdl_2.YLim = [0,1];
            ax_hdl_2.YLabel.String = 'log_{10}(s)';
            ax_hdl_2.XLabel.String = sprintf('p - p_c(s = %.3f)', s_range(1));
            legend(plt_hdl, sprintf('p_c: %.4f\nSlope: %.2f \\pm %.2f\nIntercept: %.2f \\pm %.2f\nR-squared: %.3f\nData size: %d', ...
                p_range(1), ...
                lin_log_fit.Coefficients.Estimate(2), lin_log_fit.Coefficients.SE(2), ...
                lin_log_fit.Coefficients.Estimate(1), lin_log_fit.Coefficients.SE(1), ...
                lin_log_fit.Rsquared.Adjusted, lin_log_fit.NumObservations), 'Location', 'northwest');
            
            ax_hdl_3 = subplot(1,3,3);
            scatter(ax_hdl_3, percolation_p(loglog_fit_Q) - p_range(1), ...
                percolation_str.avg.largest_cc_fraction(loglog_fit_Q));
            ax_hdl_3.XScale = 'log';
            ax_hdl_3.YScale = 'log';
            ax_hdl_3.XLabel.String = sprintf('p - p_c(s = %.3f)', s_range(1));
            ax_hdl_3.YLabel.String = 's';
            %% Save
            percolation_str.cc_ratio_2_p = cc_ratio_2_p;
            percolation_str.ratio_050_th = ratio_half_th;
            percolation_str.ratio_095_th = ratio_095_th;
            percolation_str.ratio_005_th = ratio_005_th;
            percolation_str.p_c_derivative = tmp_pc;
            
                        
            percolation_str.dataset_name = test_region_data.dataset_name;
            percolation_str.stack = test_region_data.stack;
            percolation_str.cc_ind = iter_cc;
            percolation_str.cc_label = cc_label{iter_cc};
            percolation_str.structure_name = test_region_data.structure_name;
            percolation_region_data{iter_cc, iter_stack, iter_region} = percolation_str;            
            
            tmp_fp = fullfile(DataManager.fp_visualization_folder(percolation_str.dataset_name, ...
                percolation_str.stack), fig_group_name, ...
                sprintf('%s_%s_%s_%s_%s.png', percolation_str.dataset_name, ...
                percolation_str.stack, fig_group_name, percolation_str.structure_name, ...
                percolation_str.cc_label));
            fun_print_image_in_several_formats(fig_hdl, tmp_fp);
            delete(fig_hdl);
        end
    end
end
percolation_data_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), fig_group_name, parent_str_name, ...
    sprintf('%s_%s_%s_%s_percolation_data.mat', dataset_name, ...
    merge_stack_name, fig_group_name, parent_str_name));
DataManager.write_data(percolation_data_fp, percolation_region_data);    
%% Plot percolation transition in selective region
dataset_name = 'WholeBrain';
merge_stack_name = 'all_stack';
fig_group_name = 'regional_network_percolation';
parent_str_name = 'Major_brain_regions';
percolation_data_fp = fullfile(...
    DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), fig_group_name, parent_str_name, ...
    sprintf('%s_%s_%s_%s_percolation_data.mat', dataset_name, ...
    merge_stack_name, fig_group_name, parent_str_name));
percolation_region_data = DataManager.load_data(percolation_data_fp);    



tmp_vis_region_list_ind = [1, 3, 4, 5, 6, 7, 8, 12];
percolation_vis_data = cell(numel(tmp_vis_region_list_ind), 1);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
ax_hdl = axes(fig_hdl);
hold(ax_hdl, 'on');
leg_hdl_array = [];
for iter_region = 1 : numel(percolation_vis_data)
    tmp_data = percolation_region_data(:, :, tmp_vis_region_list_ind(iter_region));
    tmp_x = cellfun(@(x) (1 - x.bond_occupancy_p), tmp_data, 'UniformOutput', false);
    tmp_y = cellfun(@(x) x.avg.largest_cc_fraction, tmp_data, 'UniformOutput', false);
    tmp_pc = cellfun(@(x) fun_analysis_percolation_compute_pc(x.bond_occupancy_p, ...
        x.avg.largest_cc_fraction), tmp_data);
    
    tmp_num_bond = cellfun(@(x) x.num_bonds, tmp_data);
    tmp_num_site = cellfun(@(x) x.num_sites, tmp_data);
    
    tmp_record = struct;
    tmp_record.avg_interpolation = fun_analysis_get_xy_curve_avgNstd_by_interpolation(...
        tmp_x, tmp_y);
    tmp_record.structure_name = tmp_data{1}.structure_name;
    tmp_record.half_th_mean = mean(tmp_pc(:));
    tmp_record.half_th_std = std(tmp_pc(:));
    tmp_record.num_bond_mean = mean(tmp_num_bond(:));
    tmp_record.num_bond_std = std(tmp_num_bond(:));
    tmp_record.num_site_mean = mean(tmp_num_site(:));
    tmp_record.num_site_std = std(tmp_num_site(:));   
    
    tmp_x = 1 - tmp_data{1}.bond_occupancy_p;
    tmp_y = tmp_record.avg_interpolation.y_avg_interpolation(tmp_x);
    
    tmp_dydx = diff(tmp_y) ./ diff(tmp_x);
    [tmp_dydx_max, tmp_max_ind] = max(tmp_dydx);
    tmp_dydx_max_x = mean(tmp_x([0, 1] + tmp_max_ind));
    tmp_y_err = tmp_record.avg_interpolation.y_std_interpolation(tmp_y);
    fprintf('Maximum gradient occupancy probability: %f\n', tmp_dydx_max_x);
%     errorbar(ax_hdl, tmp_x, tmp_y, tmp_y_err);
%     [ax_hdl, plot_hdl, patch_hld] = fun_vis_errorbar_shaded(...
%         tmp_x, tmp_y, tmp_y_err, ax_hdl);
%     plot_hdl.LineWidth = 2;
    plot_hdl = plot(ax_hdl, tmp_x, tmp_y, 'LineWidth', 2);
    plot_hdl.LineStyle = plot_line_style{tmp_line_style_ind(iter_region)};
    plot_hdl.Marker = plot_marker_style{tmp_marker_style_ind(iter_region)};
    leg_hdl_array(end+1) = plot_hdl;
%     plot(ax_hdl, tmp_x, tmp_y, 'LineWidth', 2);
    percolation_vis_data{iter_region} = tmp_record;
end

tmp_legend = cellfun(@(x) sprintf('%s :\n(# site, # bond): (%d, %d)\np_c = %.3e \\pm %.1e', ...
    strrep(x.structure_name, '_', ' '), round(x.num_site_mean), round(x.num_bond_mean), ...
    1 - x.half_th_mean, x.half_th_std), percolation_vis_data, 'UniformOutput', false);
% tmp_legend = cellfun(@(x) sprintf('%s :\n(%.1e \\pm %.1e, %.1e \\pm %.1e)\np_c = %.3e \\pm %.1e', ...
%     strrep(x.structure_name, '_', ' '), x.num_site_mean, x.num_site_std, ...
%     x.num_bond_mean, x.num_bond_std, ...
%     x.half_th_mean, x.half_th_std), percolation_vis_data, 'UniformOutput', false);

leg_hdl = legend(ax_hdl, leg_hdl_array, tmp_legend{:}, 'Location', 'best', 'NumColumns', 2);
leg_hdl.Box = 'off';
ax_hdl.YLabel.String = 'Largest cluster fraction (s)';
ax_hdl.XLabel.String = 'Link removal probability (p)';
ax_hdl.XLim = [0,1];
ax_hdl.YLim = [0,1];
ax_hdl.FontSize = 14;

fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), fig_group_name, parent_str_name, ...
    sprintf('%s_%s_%s_%s_selective_region_w_numNE.png', dataset_name, ...
    merge_stack_name, fig_group_name, parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
tmp_vis_region_list_ind = 1 : 12;
percolation_vis_data = cell(numel(tmp_vis_region_list_ind), 1);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
ax_hdl = axes(fig_hdl);
hold(ax_hdl, 'on');
leg_hdl_array = [];
for iter_region = 1 : numel(percolation_vis_data)
    tmp_data = percolation_region_data(:, :, tmp_vis_region_list_ind(iter_region));
    tmp_x = cellfun(@(x) x.bond_occupancy_p, tmp_data, 'UniformOutput', false);
    tmp_y = cellfun(@(x) x.avg.largest_cc_fraction, tmp_data, 'UniformOutput', false);
    tmp_pc = cellfun(@(x) fun_analysis_percolation_compute_pc(x.bond_occupancy_p, ...
        x.avg.largest_cc_fraction), tmp_data);
    
    tmp_num_bond = cellfun(@(x) x.num_bonds, tmp_data);
    tmp_num_site = cellfun(@(x) x.num_sites, tmp_data);
    
    tmp_record = struct;
    tmp_record.avg_interpolation = fun_analysis_get_xy_curve_avgNstd_by_interpolation(...
        tmp_x, tmp_y);
    tmp_record.structure_name = tmp_data{1}.structure_name;
    tmp_record.half_th_mean = mean(tmp_pc(:));
    tmp_record.half_th_std = std(tmp_pc(:));
    tmp_record.num_bond_mean = mean(tmp_num_bond(:));
    tmp_record.num_bond_std = std(tmp_num_bond(:));
    tmp_record.num_site_mean = mean(tmp_num_site(:));
    tmp_record.num_site_std = std(tmp_num_site(:));   
    
    tmp_x = tmp_data{1}.bond_occupancy_p;
    tmp_y = tmp_record.avg_interpolation.y_avg_interpolation(tmp_x);
    
    tmp_dydx = diff(tmp_y) ./ diff(tmp_x);
    [tmp_dydx_max, tmp_max_ind] = max(tmp_dydx);
    tmp_dydx_max_x = mean(tmp_x([0, 1] + tmp_max_ind));
    tmp_y_err = tmp_record.avg_interpolation.y_std_interpolation(tmp_y);
    fprintf('Maximum gradient occupancy probability: %f\n', tmp_dydx_max_x);
%     errorbar(ax_hdl, tmp_x, tmp_y, tmp_y_err);
%     [ax_hdl, plot_hdl, patch_hld] = fun_vis_errorbar_shaded(...
%         tmp_x, tmp_y, tmp_y_err, ax_hdl);
%     plot_hdl.LineWidth = 2;
    plot_hdl = plot(ax_hdl, tmp_x, tmp_y, 'LineWidth', 2);
    plot_hdl.LineStyle = plot_line_style{tmp_line_style_ind(iter_region)};
    plot_hdl.Marker = plot_marker_style{tmp_marker_style_ind(iter_region)};
    leg_hdl_array(end+1) = plot_hdl;
%     plot(ax_hdl, tmp_x, tmp_y, 'LineWidth', 2);
    percolation_vis_data{iter_region} = tmp_record;
end

tmp_legend = cellfun(@(x) sprintf('%s :\np_c = %.3f \\pm %.3f', ...
    strrep(x.structure_name, '_', ' '),  ...
    x.half_th_mean, x.half_th_std), percolation_vis_data, 'UniformOutput', false);
% tmp_legend = cellfun(@(x) sprintf('%s :\n(%.1e \\pm %.1e, %.1e \\pm %.1e)\np_c = %.3e \\pm %.1e', ...
%     strrep(x.structure_name, '_', ' '), x.num_site_mean, x.num_site_std, ...
%     x.num_bond_mean, x.num_bond_std, ...
%     x.half_th_mean, x.half_th_std), percolation_vis_data, 'UniformOutput', false);

leg_hdl = legend(ax_hdl, leg_hdl_array, tmp_legend{:}, 'Location', 'best', 'NumColumns', 2);
leg_hdl.Box = 'off';
ax_hdl.YLabel.String = 'Largest cluster fraction (s)';
ax_hdl.XLabel.String = 'Occupation probability (p)';
ax_hdl.XLim = [0,1];
ax_hdl.YLim = [0,1];
ax_hdl.FontSize = 14;

fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), fig_group_name, parent_str_name, ...
    sprintf('%s_%s_%s_%s_selective_region_wo_numNE.png', dataset_name, ...
    merge_stack_name, fig_group_name, parent_str_name));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
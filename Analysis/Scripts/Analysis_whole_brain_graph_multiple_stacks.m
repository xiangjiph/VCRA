set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
skl_grid_name = '240_cube_re';
output_graph_grid_name = sprintf('%s_analyzed', skl_grid_name);
grid_c_version = '240_cube_combined_5_o_2';
save_folder_name = 'whole_brain_stat';
stack_graph_cell = cell(size(stack_list));
capillary_max_radius_um = 3.5;
merge_stack_name = 'all_stack';
im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, merge_stack_name), save_folder_name);
if ~isfolder(im_save_folder)
    mkdir(im_save_folder);
end
linear_scaling_ratio = 1.0521;
radius_scaling_ratio = 1;
%%
num_stack = numel(stack_list);
stat_str_cell = cell(num_stack, 1);
for iter_stack = 1 : numel(stack_list)
    tmp_fp = fullfile(DataManager.fp_analysis_data_folder(dataset_name, stack_list{iter_stack}), save_folder_name, ...
        sprintf('%s_%s_whole_brain_stat.mat', dataset_name, stack_list{iter_stack}));
    wb_stat_str = DataManager.load_data(tmp_fp);
    wb_stat_str = fun_analysis_sc_whole_brain_stat_str(wb_stat_str, linear_scaling_ratio, radius_scaling_ratio);
    stat_str_cell{iter_stack} = wb_stat_str;
end
%% Convert whole brain statistics structure to xlsx table
stat_table_cell = cell(num_stack, 1);
scalar_field_name = {'num', 'link', 'capillary', 'graph'};
stat_field_name = {'mean', 'median', 'std', 'min', 'max', 'range', 'sum'};
for iter_stack = 1 : num_stack
   tmp_str = stat_str_cell{iter_stack}; 
   tmp_str_dist = stat_str_cell{iter_stack}.dist;
   tmp_str = rmfield(tmp_str, {'dataset_name', 'fp', 'skel_version', 'dist', 'sample_size_mm'});
   tmp_str.stack = strrep(tmp_str.stack, '_', '');
   for iter_field = 1 : numel(scalar_field_name)
       tmp_field_name = scalar_field_name{iter_field};
       tmp_substr = tmp_str.(tmp_field_name);
       tmp_subfield_name_list = fieldnames(tmp_substr);
       for iter_subfield = 1 : numel(tmp_subfield_name_list)
           tmp_subfield_name = tmp_subfield_name_list{iter_subfield};
           tmp_subfield_data = tmp_substr.(tmp_subfield_name);
           if isnumeric(tmp_subfield_data) && isscalar(tmp_subfield_data)
               tmp_new_field_name = sprintf('%s_%s', tmp_field_name, tmp_subfield_name);
               tmp_str.(tmp_new_field_name) = tmp_subfield_data;
           end
       end
       tmp_str = rmfield(tmp_str, tmp_field_name);
   end  
   tmp_str.fraction_link_wo_ep = tmp_str.num_link_wo_ep / tmp_str.num_link;
   tmp_str.num_link_2_num_node = tmp_str.num_link / tmp_str.num_node;
   % Add other   
   dist_field_name = fieldnames(tmp_str_dist);
   for iter_dist_f = 1 : numel(dist_field_name)
       tmp_dist_field_name = dist_field_name{iter_dist_f};
       for iter_stat_field = 1 : numel(stat_field_name)
           tmp_stat_field_n = stat_field_name{iter_stat_field};
           tmp_stat_field_v = tmp_str_dist.(tmp_dist_field_name).(tmp_stat_field_n);
           tmp_new_field_name = sprintf('%s_%s', tmp_dist_field_name, tmp_stat_field_n);
           tmp_str.(tmp_new_field_name) = tmp_stat_field_v;
       end
   end   
   stat_table_cell{iter_stack} = tmp_str;
end

mean_str = struct;
mean_str.stack = 'Mean';
std_str = struct;
std_str.stack = 'STD';
std_n_str = struct;
std_n_str.stack = 'CoeffVar';
compute_stat_field_list = setdiff(fieldnames(stat_table_cell{1}), {'stack'});
for iter_field = 1 : numel(compute_stat_field_list)
    tmp_field_name = compute_stat_field_list{iter_field};
    tmp_data = cellfun(@(x) x.(tmp_field_name), stat_table_cell, 'UniformOutput', false);
    tmp_data = cat(1, tmp_data{:});
    mean_str.(tmp_field_name) = mean(tmp_data, 1);
    std_str.(tmp_field_name) = std(tmp_data, 0, 1);
    std_n_str.(tmp_field_name) = std_str.(tmp_field_name) ./ mean_str.(tmp_field_name);
end
stat_table_cell = cat(1, stat_table_cell, {mean_str}, {std_str}, {std_n_str});
stat_table_cell = cellfun(@(x) struct2table(x, 'AsArray', true), stat_table_cell, 'UniformOutput', false);
stat_table = cat(1, stat_table_cell{:});

tmp_table = stat_table;
tmp_table.Properties.RowNames = tmp_table.stack;
tmp_table.stack = [];
tmp_table = rows2vars(tmp_table);
tmp_table.Properties.RowNames = tmp_table.OriginalVariableNames;
tmp_table.OriginalVariableNames = [];

merge_table_fp = fullfile(im_save_folder,...
    sprintf('%s_%s_whole_brain_stat_table.csv', dataset_name, merge_stack_name));
writetable(tmp_table, merge_table_fp, 'WriteRowNames', true);
%% Branching orders
tmp_stat_str_cell = cellfun(@(x) x.dist.cap_branching_order, stat_str_cell, 'UniformOutput', false);
tmp_hist_plot_opt = struct;
tmp_hist_plot_opt.int_x = 1 : 10;
tmp_hist_plot_opt.HistQ = true;
tmp_hist_plot_opt.ErrorBarQ = true;
tmp_hist_plot_opt.XLim = [0, 10.5];
tmp_hist_plot_opt.XLabelString = 'Capillary branch order';
tmp_hist_plot_opt.YLabelString = 'PDF';
tmp_hist_plot_opt.LegendString = cat(2, cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false), {'Average \pm STD'});
tmp_plot_data = fun_vis_plot_multi_stack_histograms(tmp_stat_str_cell, tmp_hist_plot_opt);
tmp_plot_data.ax_hdl.XTick = tmp_hist_plot_opt.int_x;
tmp_plot_data.ax_hdl.FontSize = 14;
fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_capillary_branch_order.png', ...
    dataset_name, merge_stack_name, save_folder_name));
fun_print_image_in_several_formats(tmp_plot_data.fig_hdl, fig_fp);
%% Node degree
tmp_stat_str_cell = cellfun(@(x) x.dist.node_degree, stat_str_cell, 'UniformOutput', false);
tmp_hist_plot_opt = struct;
tmp_hist_plot_opt.HistQ = true;
tmp_hist_plot_opt.ErrorBarQ = true;
tmp_hist_plot_opt.int_x = 3 : 5;
tmp_hist_plot_opt.XLim = [2.5, 5.5];
tmp_hist_plot_opt.YLim = [1e-4, 1];
tmp_hist_plot_opt.YScale = 'log';
tmp_hist_plot_opt.XLabelString = 'Node degree';
tmp_hist_plot_opt.YLabelString = 'PDF';
tmp_hist_plot_opt.LegendString = cat(2, cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false), {'Average \pm STD'});
tmp_plot_data = fun_vis_plot_multi_stack_histograms(tmp_stat_str_cell, tmp_hist_plot_opt);
tmp_plot_data.ax_hdl.XTick = 3 : 5;
tmp_plot_data.ax_hdl.FontSize = 14;
fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_node_degree.png', dataset_name, merge_stack_name, save_folder_name));
fun_print_image_in_several_formats(tmp_plot_data.fig_hdl, fig_fp);
%% Link length distribution
tmp_stat_str_cell = cellfun(@(x) x.dist.link_cap_length, stat_str_cell, 'UniformOutput', false);
tmp_hist_plot_opt = struct;
tmp_hist_plot_opt.rebin_edge = 0 : 10 : 1000;
% tmp_hist_plot_opt.int_x = 0:2.5:1000;
tmp_hist_plot_opt.cdf_limit = [1e-8, 1];
tmp_hist_plot_opt.XLim = [0, 1000];
% tmp_hist_plot_opt.YLim = [1e-8, 1];
tmp_hist_plot_opt.HistQ = true;
tmp_hist_plot_opt.ErrorBarQ = true;
tmp_hist_plot_opt.YScale = 'log';
tmp_hist_plot_opt.XScale = 'linear';
tmp_hist_plot_opt.XLabelString = 'Capillary segment length (\mum)';
tmp_hist_plot_opt.YLabelString = 'PDF';
tmp_hist_plot_opt.LegendString = cat(2, cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false), {'Average \pm STD'});
tmp_plot_data = fun_vis_plot_multi_stack_histograms(tmp_stat_str_cell, tmp_hist_plot_opt);
tmp_plot_data.ax_hdl.FontSize = 14;
tmp_plot_data.legend_hdl.Location = 'southwest';
fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_capillary_segment_length_dist.png', dataset_name, merge_stack_name, save_folder_name));
fun_print_image_in_several_formats(tmp_plot_data.fig_hdl, fig_fp);

tmp_plot_data.ax_hdl.XScale = 'log'; 
fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_capillary_segment_length_dist_loglog.png', dataset_name, merge_stack_name, save_folder_name));
fun_print_image_in_several_formats(tmp_plot_data.fig_hdl, fig_fp);
%% Radius distribution 
tmp_stat_str_cell = cellfun(@(x) x.dist.link_radius, stat_str_cell, 'UniformOutput', false);
tmp_hist_plot_opt = struct;
tmp_hist_plot_opt.rebin_edge = [1.75 : 0.5 : 10, 10 .^ (linspace(log10(10.5), log10(150), 30))];
tmp_hist_plot_opt.int_x = [2 : 0.5 : 10, 10 .^ (linspace(log10(10.5), log10(150), 30))];
% tmp_hist_plot_opt.rebin_edge = [1.8 : 0.2 : 4, 4.5 : 0.5 : 10, 10 .^ (linspace(log10(10.5), log10(150), 30))];
% tmp_hist_plot_opt.int_x = [1.9 : 0.2 : 3.9, 4.25 : 0.5 : 10, 10 .^ (linspace(log10(10.5), log10(150), 30))];

% tmp_hist_plot_opt.int_x = tmp_hist_plot_opt.rebin_edge;
tmp_hist_plot_opt.cdf_limit = [1e-8, 1];
tmp_hist_plot_opt.XLim = [0, 150];
% tmp_hist_plot_opt.YLim = [1e-8, 1];
tmp_hist_plot_opt.HistQ = true;
tmp_hist_plot_opt.ErrorBarQ = true;
tmp_hist_plot_opt.YScale = 'log';
tmp_hist_plot_opt.XScale = 'log';
tmp_hist_plot_opt.XLabelString = 'Vessel radius (\mum)';
tmp_hist_plot_opt.YLabelString = 'PDF';
tmp_hist_plot_opt.LegendString = cat(2, cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false), {'Average \pm STD'});
tmp_plot_data = fun_vis_plot_multi_stack_histograms(tmp_stat_str_cell, tmp_hist_plot_opt);
% Add radius fitting 
fit_r_list = 10 .^ (linspace(log10(3.5), log10(50), 50));
fit_prob_list = cellfun(@(x) x(fit_r_list), ...
    tmp_plot_data.interpolation_str.y_interpolation, 'UniformOutput', false);
fit_prob_list = cat(2, fit_prob_list{:});
fit_r_list = repmat(fit_r_list, 1, num_stack);
% fit_prob_list = tmp_plot_data.interpolation_str.y_avg_interpolation(fit_r_list);
fit_r_list = log10(fit_r_list);
fit_prob_list = log10(fit_prob_list);

loglog_fit_hdl = fitlm(fit_r_list, fit_prob_list);
hold(tmp_plot_data.ax_hdl, 'on');
tmp_plot_data.loglog_fit_hdl = plot(tmp_plot_data.ax_hdl, ...
    10 .^ (fit_r_list), 10 .^ (loglog_fit_hdl.Coefficients.Estimate(1) + ...
    loglog_fit_hdl.Coefficients.Estimate(2) .* fit_r_list));
tmp_plot_data.loglog_fit_hdl.LineWidth = 3;

tmp_plot_data.legend_hdl.String{num_stack + 2} = sprintf('Exponent: %.2f \\pm %.2f\nR-squared: %.3f', ...
    loglog_fit_hdl.Coefficients.Estimate(2), loglog_fit_hdl.Coefficients.SE(2), ...
    loglog_fit_hdl.Rsquared.Adjusted);
tmp_plot_data.ax_hdl.FontSize = 14;
tmp_plot_data.legend_hdl.Location = 'southwest';

fig_fp = fullfile(im_save_folder, sprintf('%s_%s_%s_vessel_radius_w_scaling.png', dataset_name, merge_stack_name, save_folder_name));
fun_print_image_in_several_formats(tmp_plot_data.fig_hdl, fig_fp);
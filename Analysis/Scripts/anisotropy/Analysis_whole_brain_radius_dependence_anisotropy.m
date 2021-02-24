set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
num_stack = numel(stack_list);
stack_name_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
grid_version = '240_cube';
data_folder = 'Radius_dependence_anisotropy';
skel_version = '240_cube_rec';
recon_version = '240_cube_recon_sc';
whole_brain_stat_ver = 'whole_brain_stat_sc';
vis_folder_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, 'all_stack'), ...
    'Radius_dependence_anisotropy');
%% Load data 
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    stack = stack_list{iter_stack};
    wb_cube_stat = DataManager.load_data(DataManager.fp_analysis_data_file(...
        dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
        dataset_name, stack, recon_version), whole_brain_stat_ver));
    wb_cube_stat.rda_data = DataManager.load_data(DataManager.fp_analysis_data_file(...
        dataset_name, stack, sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, ...
        stack, skel_version, data_folder)));
    wb_data_cell{iter_stack} = wb_cube_stat;
end
fprintf('Finish loading whole brain data from %d stacks\n', num_stack);
%% Combined data from stacks
cap2vsl_vol_f_min = 0.0;
im_save_folder_fp = fullfile(vis_folder_fp, sprintf('internal_cube_cap2vsl_min_vf%.1f', cap2vsl_vol_f_min));
plot_info_cell = cell(2, 0);
plot_info_cell(:, end+1) = {1, stack_name_list{1}};
plot_info_cell(:, end+1) = {2, stack_name_list{2}};
plot_info_cell(:, end+1) = {3, stack_name_list{3}};
plot_info_cell(:, end+1) = {[1, 2], 'Merge ML20180815 and ML20190124'};
plot_info_cell(:, end+1) = {[1, 2, 3], 'all_stacks'};
num_plot = size(plot_info_cell, 2);
weight_method = 'length_weighted';
for iter_plot = 1 : num_plot
    vis_stack_name = plot_info_cell{2, iter_plot};
    vis_stack_idx = plot_info_cell{1, iter_plot};
    
    vis_radius_idx = 1 : 16;
    vis_radius_min = wb_data_cell{1}.rda_data.radius_min(vis_radius_idx);
    vis_radius_max = wb_data_cell{1}.rda_data.radius_max(vis_radius_idx);
    num_vis_radius = numel(vis_radius_idx);
    
    is_internal_cube_Q = cellfun(@(x) x.cube_in_brain_mask_ratio == 1, wb_data_cell(vis_stack_idx), 'UniformOutput', false);
    is_internal_cube_Q = cat(1, is_internal_cube_Q{:});
    
    is_mainly_capillary_Q = cellfun(@(x) x.cap2vsl_vol_fraction >= cap2vsl_vol_f_min , wb_data_cell(vis_stack_idx), 'UniformOutput', false);
    is_mainly_capillary_Q = cat(1, is_mainly_capillary_Q{:});
    is_internal_cube_Q = is_internal_cube_Q & is_mainly_capillary_Q;
    
    fa_merged = cellfun(@(x) x.rda_data.(weight_method).fa, wb_data_cell(vis_stack_idx), 'UniformOutput', false);
    fa_merged = cat(1, fa_merged{:});
    fa_p_merged = cellfun(@(x) x.rda_data.(weight_method).fa_p, wb_data_cell(vis_stack_idx), 'UniformOutput', false);
    fa_p_merged = cat(1, fa_p_merged{:});
    fa_z_merged = cellfun(@(x) x.rda_data.(weight_method).fa_z, wb_data_cell(vis_stack_idx), 'UniformOutput', false);
    fa_z_merged = cat(1, fa_z_merged{:});
    weight_sum_merged = cellfun(@(x) x.rda_data.(weight_method).weight_sum, wb_data_cell(vis_stack_idx), 'UniformOutput', false);
    weight_sum_merged = cat(1, weight_sum_merged{:});    
    weight_sum_merged = weight_sum_merged ./ weight_sum_merged(:, end-1);
    %% Combine all the stack data
    % Radius dependence of FA, FA_z and FA_p
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [1, 1.5];
    ax_hdl_0 = subplot(4, 1, 1);
    boxplot(ax_hdl_0, weight_sum_merged(is_internal_cube_Q, vis_radius_idx), 'Symbol', '');
    ax_hdl_0.YLim = [0, 1];
    ax_hdl_0.YLabel.String = sprintf('Fraction of %s', weight_method(1:6));
    ax_hdl_0.Box = 'on';
    grid(ax_hdl_0, 'on');
    ax_hdl_0.XAxis.Visible = 'off';
    ax_hdl_1 = subplot(4, 1, 2);
    boxplot(ax_hdl_1, fa_merged(is_internal_cube_Q, vis_radius_idx), 'Symbol', '');
    ax_hdl_1.YLim = [0, 1];
    ax_hdl_1.XAxis.Visible = 'off';
    ax_hdl_1.YLabel.String = 'FA';
    ax_hdl_1.Box = 'on';
    grid(ax_hdl_1, 'on');
    ax_hdl_2 = subplot(4, 1, 3);
    boxplot(ax_hdl_2, fa_z_merged(is_internal_cube_Q, vis_radius_idx), 'Symbol', '');
    ax_hdl_2.YLim = [-5, 10];
    ax_hdl_2.XAxis.Visible = 'off';
    ax_hdl_2.YLabel.String = 'FA_z';
    ax_hdl_2.Box = 'off';
    grid(ax_hdl_2, 'on');
    ax_hdl_3 = subplot(4, 1, 4);
    tmp_p_value = fa_p_merged(is_internal_cube_Q, vis_radius_idx);
    tmp_p_value(tmp_p_value == 0) = 5e-5;
    boxplot(ax_hdl_3, tmp_p_value, 'Symbol', '');
    ax_hdl_3.YLim = [1e-5, 1];
    ax_hdl_3.YScale = 'log';
    ax_hdl_3.YLabel.String = 'FA_p';
    ax_hdl_3.YTick = [1e-4, 1e-2, 1];
    ax_hdl_3.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), vis_radius_max, 'UniformOutput', false);
    ax_hdl_3.XLabel.String = 'Radius upper limit (\mum)';
    ax_hdl_3.Box = 'on';
    grid(ax_hdl_3, 'on');
    ax_hdl_3.YMinorGrid = 'off';
    ax_hdl_0.Title.String = vis_stack_name;
    ax_hdl_0.Title.Interpreter = 'none';
    %%
    tmp_fig_fp = fullfile(im_save_folder_fp, sprintf('%s_%s_%s_anisotropy_vs_radius_upper_lim.png', ...
        dataset_name, vis_stack_name, weight_method));
    fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
end
%% 
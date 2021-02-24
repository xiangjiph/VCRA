%% Document
% This script work for individual image stack
% 1. Plot PDF and CDF of local vessel network properties
% 2. Plot box plot of local vessel network properties distribution
% 3. This script should be run after
% Analysis_preprocess_whole_brain_local_statistics .
%%
set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
% stack = 'ML20190124';
skl_grid_name = '240_cube_auto';
reconstruction_name = '240_cube_recon_v2';
image_grid_name = '240_cube';

grid_c_version = '240_cube_combined_5_o_2';
grid_c_info = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c_info.grid_ori;
wb_stat_folder_name = 'whole_brain_stat_v2';
save_im_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), wb_stat_folder_name, 'internal_cubes');
if ~isfolder(save_im_folder)
    mkdir(save_im_folder);
end
%% Load whole brain 240-cube statistics
tic
wb_240_cube_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
    dataset_name, stack, reconstruction_name), wb_stat_folder_name);
wb_240_cube_pdf_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_pdf.mat', ...
    dataset_name, stack, reconstruction_name), wb_stat_folder_name);
fprintf('Finish loading 240 cube data\n');
toc
is_internal_240_cube_ratio = (wb_240_cube_stat_str.cube_in_brain_mask_ratio == 1);
%% Initialize structure for recording plotting information
stat_str = struct;
stat_str.info.dataset_name = dataset_name;
stat_str.info.stack = stack;
stat_str.info.process_date = datestr(now);
stat_str.info.filepath = fullfile(save_im_folder, sprintf('%s_%s_%s_%s_plot_data.mat', dataset_name, stack, wb_stat_folder_name, 'internal_cubes'));
stat_str.feature_scalar_stat = [];
stat_str.cell_data = [];
stat_str.internal_cube_ind = find(is_internal_240_cube_ratio);
stat_str.cube_volume_ratio_in_brain_mask = is_internal_240_cube_ratio;
stat_str.feature_scalar_plot_setting = [];
stat_str.feature_box_plot_setting = [];
fprintf('Finish initializing the structures\n');
%% Histograms - whole brain statistics - scalar data
fprintf('Generating plot setting cell array for scalar data\n');
plot_info_cell = cell(3, 0);
% Segments and node density
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.num_data.length ./ 0.240^3, 'Capillary_number_density', 'Capillary number density (mm^{-3})'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.num_data.length ./ 0.240^3, 'Vessel_number_density', 'Vessel number density (mm^{-3})'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.num_data.degree ./ 0.240^3, 'Node_density', 'Node density(mm^{-3})'};

% Anisotropy
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.fractional_anisotropy, 'Vessel_volume-weighted_fractional_anisotropy', 'Fractioanal Anisotropy'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.fa_z, 'Vessel_volume-weighted_FA_z-score', 'FA_z'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.fa_p, 'Vessel_volume-weighted_FA_p-Value', 'FA p-Value'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.svd_value_ratio(:, 1), 'Vessel_volume-weighted_PCV1', 'PCV1'};

plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.fractional_anisotropy, 'Capillary_volume-weighted_fractional_anisotropy', 'Fractioanal Anisotropy'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.fa_z, 'Capillary_volume-weighted_FA_z-score', 'FA_z'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.fa_p, 'Capillary_volume-weighted_FA_p-Value', 'FA p-Value'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.svd_value_ratio(:, 1), 'Capillary_volume-weighted_PCV1', 'PCV1'};
% Capillary - vessel ratio
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_length_fraction, 'Capillary-vessel_length_ratio', 'Local capillary length ratio'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_vol_fraction, 'Capillary-vessel_volume_ratio', 'Local capillary volume ratio'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_surf_area_fraction, 'Capillary-vessel_surface_area_ratio', 'Local capillary surface area ratio'};

% Global statistics
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.mask_volume_density, 'Reconstructed_vessel_volume_density', 'Vessel volume density'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.mask_surface_area_density_mm2_mm3, 'Reconstructed_vessel_surface_area', 'Vessel surface area density(mm^2/mm^3)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_length_density_m_mm3, 'Vessel_length_density', 'Vessel length density(m/mm^3)'};

plot_info_cell(:, end+1) = {wb_240_cube_stat_str.capillary_volume_density, 'Capillary_volume_density', 'Capillary volume density'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.capillary_length_density_m_mm3, 'Capillary_length_density', 'Capillary length density(m/mm^3)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.capillary_surface_area_density_mm2_mm3, 'Capillary_surface_area_density', 'Capillary surface area density(mm^2/mm^3)'};

% capillary segment statistics
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.length, 'Median_capillary_segment_length', 'Median capillary segment length(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.length, 'Average_capillary_segment_length', 'Average capillary segment length(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max, 'Median_capillary_dmax', 'Median $d_{max} (\mu m)$'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.nearest_tissue_dt_max, 'Mean_capillary_dmax', 'Average $d_{max} (\mu m)$'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.nearest_tissue_radius, 'Mean_capillary_nearest_tissue_radius', 'Mean capillary tissue radius(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_radius, 'Median_capillary_nearest_tissue_radius', 'Median capillary tissue radius(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.dt_median, 'Mean_capillary_radius', 'Local average capillary radius (\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.median.dt_median, 'Median_vessel_radius', 'Local median vessel radius (\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.dt_median, 'Mean_vessel_radius', 'Local average vessel radius (\mum)'};
% Loop
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_geodesic_length, 'Average_number_of_edges_in_the_shortest_loop_of_capillary', 'Number of edges in shortest loop'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_length, 'Average_length_of_the_shortest_loop_of_capillary', 'Length of shortest loop(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.std.shortest_loop_geodesic_length, 'STD_number_of_edges_in_the_shortest_loop_of_capillary', 'STD of number of edges in shortest loop'};

% Branching order
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.capillary_branching_order, 'Mean_capillary_branching_order', 'Capillary branching order'};
% Distance to large vessels
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.dist_to_nearest_noncapillary_mean, 'Median_distance_to_nearest_noncapillary', 'Median distance to nearest noncapillary(\mum)'};
% Geometry
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.straightness, 'Mean_capillary_straightness', 'Straightness'};
% Noncapillary
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_tissue_radius, 'Mean_noncapillary_nearest_tissue_radius' 'Average noncapillary tissue radius(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_capillary_num_vxl, 'Mean_nearest_caillary_total_num_skel_vxl', 'Average number of nearest capillary skeleton voxels'};

% Node statistics
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.mean.degree, 'Mean_node_degree', 'Node degree'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.nearest_node_dist, 'Median_distance_to_the_nearest_node', 'Median distance to the nearest node(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.mean.path_to_nearest_neighbor_geodesic_length, 'Mean_geodesic_distance_to_nearest_node', 'Number of edges to the nearest node'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_min, 'Median_length_of_the_shortest_link_of_node', 'Median of the connected shortest links length(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_max, 'Median_length_of_the_longest_link_of_node', 'Median of the connected longest link length(\mum)'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_median, 'Median_length_of_the median_link_of_node', 'Median of the connected median link length(\mum)'};

% Segmentation quality
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q, 'Fraction_of_vessel_without_endpoint', 'Local fraction of vessel w/o endpoint'};
plot_info_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.has_no_ep_Q, 'Fraction_of_capillary_without_endpoint', 'Local fraction of capillary w/o endpoint'};

stat_str.feature_scalar_plot_setting = plot_info_cell;
%% Boxplot - whole brain statistics - cell data
fprintf('Generating plot setting cella array for nonscalar data\n');
box_data_extraction_info = cell(7, 0);
box_data_extraction_info(:, end+1) = {'node_stat', 'path_to_nearest_neighbor_geodesic_length', 1 : 6, 'Number of edges to the nearest node', 'linear', 'log', []};
box_data_extraction_info(:, end+1) = {'link_all_stat', 'shortest_loop_geodesic_length', 2 : 16, 'Number of edges in the shortest loop', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'shortest_loop_geodesic_length', 2 : 16, 'Number of edges in the shortest loop', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'capillary_branching_order', 1 : 9, 'Capillary branching order', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'degree', 3 : 5, 'Node degree', 'linear', 'log', []};
num_box_data = size(box_data_extraction_info, 2);

for iter_box = 1 : num_box_data
    tmp_field_name = box_data_extraction_info{1, iter_box};
    tmp_feature_name = box_data_extraction_info{2, iter_box};
    tmp_bin_val = box_data_extraction_info{3, iter_box};
    tmp_data = wb_240_cube_pdf_str.(tmp_field_name).hist_pdf.(tmp_feature_name)(is_internal_240_cube_ratio == 1);
    tmp_y_bin_edge = wb_240_cube_pdf_str.(tmp_field_name).hist_edge.(tmp_feature_name)(is_internal_240_cube_ratio == 1);
    tmp_bin_data_str = fun_analysis_bin_network_properties_pdf(tmp_data, tmp_y_bin_edge, tmp_bin_val);
    box_data_extraction_info{end, iter_box} = tmp_bin_data_str.data;
    tmp_stat_data = rmfield(tmp_bin_data_str, 'data');
    tmp_stat_data.data_name = sprintf('%s_%s', tmp_field_name, tmp_feature_name);
    tmp_stat_data.data_label = box_data_extraction_info{4, iter_box};
    stat_str.cell_data{end+1} = tmp_stat_data;
end
stat_str.feature_box_plot_setting = box_data_extraction_info;
fprintf('Finish constructing plot information structure\n');
%% Plot
tic
scalar_pdf_cdf_plot_Q = false;
feature_scalar_plot_stat = cell(size(plot_info_cell, 2), 1);
for iter_vis = 1 : size(plot_info_cell, 2)
    vis_data = plot_info_cell{1, iter_vis};
    vis_data = vis_data(is_internal_240_cube_ratio == 1);
    vis_data = vis_data(~isnan(vis_data));
    
    tmp_stat = fun_analysis_get_basic_statistics(vis_data);
    [tmp_pdf, tmp_y_bin_edge] = fun_analysis_select_histcount_edge_by_percentile(tmp_stat.hist_pdf, ...
        tmp_stat.hist_edge, [5e-3, 1 - 1e-2], 'pdf');
    [tmp_cdf, tmp_1] = fun_analysis_select_histcount_edge_by_percentile(tmp_stat.hist_cdf, ...
        tmp_stat.hist_edge, [5e-3, 1 - 1e-2], 'cdf');
    assert(all(tmp_1 == tmp_y_bin_edge));
    
    tmp_leg_str = fun_analysis_basic_stat_str_to_string(tmp_stat);
    if scalar_pdf_cdf_plot_Q
        fig_hdl = figure('Visible', 'off');
        fig_hdl.Position(3:4) = fig_hdl.Position(3:4).* 2;
        ax_hdl = axes(fig_hdl); %#ok<LAXES>
        yyaxis(ax_hdl, 'right');
        histogram(ax_hdl, 'BinCounts', tmp_cdf, 'BinEdges', tmp_y_bin_edge);
        ax_hdl.YLabel.String = 'CDF';
        yyaxis(ax_hdl, 'left');
        histogram(ax_hdl, 'BinCounts', tmp_pdf, 'BinEdges', tmp_y_bin_edge);
        if ~isempty(strfind(plot_info_cell{3, iter_vis}, '$'))
            ax_hdl.XLabel.Interpreator = 'latex';
        else
            ax_hdl.XLabel.Interpreator = 'tex';
        end
        ax_hdl.XLabel.String = plot_info_cell{3, iter_vis};
        ax_hdl.YLabel.String = 'PDF';
        ax_hdl.FontSize = 14;
        ax_hdl.FontWeight = 'bold';
        ax_hdl.Title.String = strrep(plot_info_cell{2, iter_vis}, '_', ' ');
        grid(ax_hdl, 'on');
        leg_hdl = legend(tmp_leg_str, 'Location', 'best');
        fig_name = fullfile(save_im_folder, sprintf('%s_whole_brain_internal_cubes_%s.png', ...
            stack, plot_info_cell{2, iter_vis}));
        fun_print_image_in_several_formats(fig_hdl, fig_name);
        delete(fig_hdl);
    end
    tmp_stat.data_name = strrep(plot_info_cell{2, iter_vis}, '_', ' ');
    tmp_stat.data_label = plot_info_cell{3, iter_vis};
    feature_scalar_plot_stat{iter_vis} = tmp_stat;
end
fprintf('Finish writing saclar histogram figures\n');
toc
stat_str.feature_scalar_stat = feature_scalar_plot_stat;
%% Plot boxplot
boxplot_Q = true;
box_plot_output = cell(size(stat_str.feature_box_plot_setting, 2), 1);
tic
for iter_plot = 1 : size(stat_str.feature_box_plot_setting, 2)
    tmp_data = stat_str.feature_box_plot_setting{7, iter_plot};
    tmp_xlabel = stat_str.feature_box_plot_setting{4, iter_plot};
    tmp_tick = stat_str.feature_box_plot_setting{3, iter_plot};
    tmp_feature_name = sprintf('%s_%s', stat_str.feature_box_plot_setting{1, iter_plot}, ...
        stat_str.feature_box_plot_setting{2, iter_plot});
    tmp_xscale = stat_str.feature_box_plot_setting{5, iter_plot};
    tmp_yscale = stat_str.feature_box_plot_setting{6, iter_plot};
    if boxplot_Q
        fig_hdl = figure('Visible', 'on');
        fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
        ax_hdl = axes(fig_hdl);
        box_hdl = boxplot(ax_hdl, tmp_data);
        ax_hdl.XLabel.String = tmp_xlabel;
        ax_hdl.YLabel.String = 'PDF';
        ax_hdl.FontSize = 14;
        ax_hdl.FontWeight = 'bold';
        ax_hdl.YScale = tmp_yscale;
        ad_hdl.XScale = tmp_xscale;
        ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%d'), tmp_tick, 'UniformOutput', false);
        grid(ax_hdl, 'on');
        fig_fp = fullfile(save_im_folder, sprintf('%s_whole_brain_internal_cubes_%s_pdf_boxplot.png', ...
            stack, tmp_feature_name));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
        delete(fig_hdl);
        box_plot_output{iter_plot} = box_hdl;
    end
end
fprintf('Finish writing box plot figures\n');
toc
stat_str.feature_box_stat = box_plot_output;
%% Correlations matrix
cube_feature_stat = {wb_240_cube_stat_str.link_cap_stat.median.length, wb_240_cube_stat_str.link_cap_stat.median.surface_area, ...
    wb_240_cube_stat_str.link_cap_stat.median.shortest_loop_length, wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_geodesic_length, ...
    wb_240_cube_stat_str.link_cap_stat.median.straightness, ...
    wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_volume .^ (1/3), wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max, ...
    wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_radius, ...
    wb_240_cube_stat_str.link_cap_stat.std.length, wb_240_cube_stat_str.link_cap_stat.std.nearest_tissue_volume .^ (1/3)...
    wb_240_cube_stat_str.link_cap_stat.num_data.length .^ (1/3), wb_240_cube_stat_str.node_stat.num_data.degree .^ (1/3), ...
    wb_240_cube_stat_str.capillary_length_density_m_mm3.^ (1/2), ...
    wb_240_cube_stat_str.node_stat.mean.path_to_nearest_neighbor_geodesic_length, ...
    wb_240_cube_stat_str.node_stat.median.nearest_node_dist, ...
    wb_240_cube_stat_str.wb_ai_all_vw.fractional_anisotropy, wb_240_cube_stat_str.wb_ai_cap_vw.fractional_anisotropy, ...
    wb_240_cube_stat_str.wb_ai_cap_vw.fa_z};

corr_feature_name = {'Median capillary length', 'Median capillary surface area', ...
    'Median shortest loop length', 'Mean number of edge in shortest loop', ...
    'Median capillary straightness', ...
    '(Median nearest tissue volume)^{1/3}', 'Median nearest tissue DT max',...
    'Median nearest tissue radius', ...
    'STD Capillary length', '(STD nearest tissue volume)^{1/3}'...
    '(Number of capillaries)^{1/3}', '(Number of nodes)^{1/3}'...
    '(Capillary length density)^{1/2}', ...
    'Number of edges to nearest node', ...
    'Distance to nearest node', ...
    'Vessel fractional anisotropy', 'Capillary fractional anisotropy', ...
    'Capillary FA_z'};
num_cand = numel(cube_feature_stat);
%% Plot the correlation matrix
corr_mat = eye(num_cand);
selected_Q = is_internal_240_cube_ratio == 1;
for iter_1 = 1 : num_cand
    for iter_2 = (iter_1 + 1) : num_cand
        tmp_X = cube_feature_stat{iter_1};
        tmp_Y = cube_feature_stat{iter_2};
        tmp_is_selected = ~isnan(tmp_X) & ~isnan(tmp_Y) & selected_Q;
        tmp_X = tmp_X(tmp_is_selected);
        tmp_Y = tmp_Y(tmp_is_selected);
        tmp_corr = corrcoef(tmp_X, tmp_Y);
        corr_mat(iter_1, iter_2) = tmp_corr(1, 2);
        corr_mat(iter_2, iter_1) = tmp_corr(1, 2);
    end
end
fig_hdl = figure('Position', [1, 1, 2048, 2048]);
ax = axes(fig_hdl);
imagesc(ax, corr_mat);
cbar = colorbar;
cbar.Label.String = 'Correlation';
cbar.Label.FontSize = 14;
cbar.Limits = [-1, 1];
ax.Colormap = jet;
daspect([1,1,1]);
ax.XAxis.TickValues = 1 : num_cand;
ax.XAxis.TickLabelRotation = 90;
ax.XAxis.TickLabels = corr_feature_name;
ax.YAxis.TickValues = 1 : num_cand;
ax.YAxis.TickLabels = corr_feature_name;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
ax.Title.String = 'Whole brain internal cube local network properties correlation matrix';
ax.Title.FontSize = 18;
fun_print_image_in_several_formats(fig_hdl, fullfile(save_im_folder, 'Whole_brain_internal_cube_local_network_properties_correlation_matrix.png'));

stat_str.feature_corr_mat.value = corr_mat;
stat_str.feature_corr_mat.feature_name = corr_feature_name;
stat_str.feature_corr_mat.cube_selected_Q = selected_Q;
%% Data for linear regression
linear_fit_cell = cell(3, 0);
% Anisotropy
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.fractional_anisotropy, 'Vessel_volume-weighted_fractional_anisotropy', 'Vessel volume-weighted fractioanal anisotropy'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_all_vw.svd_value_ratio(:, 1), 'Vessel_volume-weighted_PCV1', 'Vessel volume-weighted PCV1'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.fractional_anisotropy, 'Capillary_volume-weighted_fractional_anisotropy', 'Capillary volume-weighted fractioanal anisotropy'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.wb_ai_cap_vw.svd_value_ratio(:, 1), 'Capillary_volume-weighted_PCV1', 'Capillary volume-weighted PCV1'};
% Capillary - vessel ratio
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_length_fraction, 'Capillary-vessel_length_ratio', 'Length ratio'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_vol_fraction, 'Capillary-vessel_volume_ratio', 'Volume ratio'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.cap2vsl_surf_area_fraction, 'Capillary-vessel_surface_area_ratio', 'Surface area ratio'};

% Global statistics
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.mask_volume_density, 'Reconstructed_vessel_volume_density', 'Vessel volume density'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.mask_surface_area_density_mm2_mm3, 'Reconstructed_vessel_surface_area', 'Vessel surface area density(mm^2/mm^3)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_length_density_m_mm3, 'Vessel_length_density', 'Vessel length density(m/mm^3)'};
linear_fit_cell(:, end+1) = {(wb_240_cube_stat_str.link_length_density_m_mm3).^(1/2), 'Vessel_length_density_squared_root', '(Vessel length density)^{1/2}'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.capillary_volume_density, 'Capillary_volume_density', 'Capillary volume density'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.capillary_length_density_m_mm3, 'Capillary_length_density', 'Capillary length density(m/mm^3)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.capillary_surface_area_density_mm2_mm3, 'Capillary_surface_area_density', 'Capillary surface area density(mm^2/mm^3)'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.capillary_length_density_m_mm3 .^ (1/2), 'Capillary_length_density_squared_root', '(Capillary length density)^{1/2}'};
% capillary segment statistics
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.length, 'Median_capillary_segment_length', 'Median capillary segment length(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.length, 'Average_capillary_segment_length', 'Average capillary segment length(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max, 'Median_maximum_distance_between_tissue_and_the_nearest_capillary', 'Median capillary-tissue distance(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.nearest_tissue_radius, 'Mean_capillary_nearest_tissue_radius', 'Mean capillary tissue radius(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_radius, 'Median_capillary_nearest_tissue_radius', 'Median capillary tissue radius(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.num_data.length, 'Number_of_capillaries', 'Number of capillary'};
linear_fit_cell(:, end+1) = {(wb_240_cube_stat_str.link_cap_stat.num_data.length).^(1/3), 'Number_of_capillaries_cubic_root', '(Number of capillary)^{1/3}'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_volume, 'Median_capillary_nearest_tissue_volume', 'Median capillary nearest tissue volume(\mum^3)'};
linear_fit_cell(:, end+1) = {(wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_volume).^(1/3), 'Median_capillary_nearest_tissue_volume_cubic_root', '(Median capillary nearest tissue volume)^{1/3}(\mum)'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.std.nearest_tissue_dt_max, 'STD_capillary_nearest_tissue_DT_max', 'STD of maximum tissue-capillary distance(\mum)'};
linear_fit_cell(:, end+1) = {(wb_240_cube_stat_str.link_cap_stat.std.nearest_tissue_volume) .^ (1/3), 'STD_capillary_nearest_tissue_volume_cubic_root', '(STD of capillary nearest tissue volume)^{1/3}(\mum)'};

% Loop
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_geodesic_length, 'Average_number_of_edges_in_the_shortest_loop_of_capillary', 'Number of edges in shortest loop'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_length, 'Average_length_of_the_shortest_loop_of_capillary', 'Length of shortest loop(\mum)'};
% Distance to large vessels
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.dist_to_nearest_noncapillary_mean, 'Median_distance_to_nearest_noncapillary', 'Median distance to nearest noncapillary(\mum)'};
% Geometry
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.straightness, 'Mean_capillary_straightness', 'Straightness'};
% Noncapillary
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_tissue_radius, 'Mean_noncapillary_nearest_tissue_radius' 'Average noncapillary tissue radius(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_all_stat.mean.noncapillary_nearest_capillary_num_vxl, 'Mean_nearest_caillary_total_num_skel_vxl', 'Average number of nearest capillary skeleton voxels'};

% Node statistics
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.nearest_node_dist, 'Median_distance_to_the_nearest_node', 'Median distance to the nearest node(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.mean.path_to_nearest_neighbor_geodesic_length, 'Mean_geodesic_distance_to_nearest_node', 'Number of edges to the nearest node'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.num_data.degree ./ 0.240^3, 'Node_density', 'Node density(mm^{-3})'};
linear_fit_cell(:, end+1) = {(wb_240_cube_stat_str.node_stat.num_data.degree ./ 0.240^3) .^(1/3), 'Node_density_cubic_root', '(Node density)^{1/3}(mm^{-1})'};

linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_min, 'Median_length_of_the_shortest_link_of_node', 'Median of the connected shortest links length(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_max, 'Median_length_of_the_longest_link_of_node', 'Median of the connected longest link length(\mum)'};
linear_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_median, 'Median_length_of_the median_link_of_node', 'Median of the connected median link length(\mum)'};
%% Plot linear regression between pairs of local network features
% Remove about 4% of the data for improving the linear regression
lr_selected_Q = wb_240_cube_stat_str.link_cap_stat.num_data.length > 50 & wb_240_cube_stat_str.link_cap_stat.num_data.length < 1000 &...
    is_internal_240_cube_ratio == 1 & wb_240_cube_stat_str.link_cap_stat.median.length > 35 & wb_240_cube_stat_str.mask_volume_density < 0.1 & ...
    wb_240_cube_stat_str.link_cap_stat.median.length < 120;
num_linear_reg_features = size(linear_fit_cell, 2);

[fit_slope_list, fit_intercept_list, fit_rsquared_adjust_list] = deal(nan(num_linear_reg_features));
linear_fit_plot_Q = false;
for iter_fix_ind = 1 : num_linear_reg_features
    tmp_fix_feature_data = linear_fit_cell(:, iter_fix_ind);
    for iter_feature = 1 : num_linear_reg_features
        if iter_feature == iter_fix_ind
            continue;
        end
        tmp_fit_feature_data = linear_fit_cell(:, iter_feature);
        
        tmp_X = tmp_fix_feature_data{1};
        tmp_Y = tmp_fit_feature_data{1};
        
        tmp_is_selected = ~isnan(tmp_X) & ~isnan(tmp_Y) & lr_selected_Q;
        tmp_X = tmp_X(tmp_is_selected);
        tmp_Y = tmp_Y(tmp_is_selected);
        linear_fit_hdl = fitlm(tmp_X, tmp_Y);
        fit_slope_list(iter_fix_ind, iter_feature) = linear_fit_hdl.Coefficients.Estimate(2);
        fit_intercept_list(iter_fix_ind, iter_feature) = linear_fit_hdl.Coefficients.Estimate(1);
        fit_rsquared_adjust_list(iter_fix_ind, iter_feature) = linear_fit_hdl.Rsquared.Adjusted;
        if linear_fit_hdl.Rsquared.Adjusted > 0.5
            % Plot the fiture
            if linear_fit_plot_Q
                fig_hdl = figure('Visible', 'on');
                fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
                ax_hdl = axes(fig_hdl); %#ok<LAXES>
                histogram2(ax_hdl, tmp_X, tmp_Y, 'DisplayStyle', 'tile');
                ax_hdl.XLabel.String = tmp_fix_feature_data{3};
                ax_hdl.YLabel.String = tmp_fit_feature_data{3};
                ax_hdl.FontSize = 14;
                ax_hdl.FontWeight = 'bold';
                ax_hdl.ColorScale = 'log';
                hold(ax_hdl, 'on');
                fit_data_x = linspace(min(tmp_X, [], 'all'), max(tmp_X, [], 'all'), 30);
                plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
                    'LineWidth', 2, 'Color', 'k');
                cbar_hdl = colorbar;
                cbar_hdl.Label.String = 'Number of data points';
                leg_hdl = legend(plt_hdl, sprintf('Slope: %.2e \\pm %.1e\nIntercept: %.2e \\pm %.1e\nR-squared: %.2f\nData size: %d', ...
                    linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
                    linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
                    linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations), 'Location', 'best');
                
                
                
                fun_print_image_in_several_formats(fig_hdl, fullfile(save_im_folder, 'cube_feature_linear_fit', ...
                    tmp_fix_feature_data{2}, sprintf('%s_%s_int_cube_feature_fit_%s_to_%s.png', ...
                    dataset_name, stack, tmp_fix_feature_data{2}, tmp_fit_feature_data{2})));
                delete(fig_hdl);
            end
        end
    end
end
stat_str.int_cap_feat_fit.slope = fit_slope_list;
stat_str.int_cap_feat_fit.intercept = fit_intercept_list;
stat_str.int_cap_feat_fit.rsquared_adjust = fit_rsquared_adjust_list;
stat_str.int_cap_feat_fit.feature_name = linear_fit_cell(2, :);
stat_str.int_cap_feat_fit.cube_selected_Q = lr_selected_Q;
%% Log - log fit - network length properties
log_log_fit_cell = cell(3, 0);
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_length_density_m_mm3, 'Vessel_length_density', 'Vessel length density(m/mm^3)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.capillary_length_density_m_mm3, 'Capillary_length_density', 'Capillary length density(m/mm^3)'};
% capillary segment statistics
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.length, 'Median_capillary_segment_length', 'Median capillary segment length(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.length, 'Average_capillary_segment_length', 'Average capillary segment length(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max, 'Median_maximum_distance_between_tissue_and_the_nearest_capillary', 'Median of maximum capillary-tissue distance(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.nearest_tissue_radius, 'Mean_capillary_nearest_tissue_radius', 'Mean capillary tissue radius(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_radius, 'Median_capillary_nearest_tissue_radius', 'Median capillary tissue radius(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_volume, 'Median_capillary_nearest_tissue_volume', 'Median capillary nearest tissue volume(\mum^3)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.std.nearest_tissue_dt_max, 'STD_capillary_nearest_tissue_DT_max', 'STD of maximum tissue-capillary distance(\mum)'};

log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.num_data.length, 'Number_of_capillaries', 'Number of capillary'};

log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.link_cap_stat.mean.shortest_loop_length, 'Average_length_of_the_shortest_loop_of_capillary', 'Length of shortest loop(\mum)'};

% Node statistics
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.nearest_node_dist, 'Median_distance_to_the_nearest_node', 'Median distance to the nearest node(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.num_data.degree ./ 0.240^3, 'Node_density', 'Node density(mm^{-3})'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_min, 'Median_length_of_the_shortest_link_of_node', 'Median of the connected shortest links length(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_max, 'Median_length_of_the_longest_link_of_node', 'Median of the connected longest link length(\mum)'};
log_log_fit_cell(:, end+1) = {wb_240_cube_stat_str.node_stat.median.link_length_median, 'Median_length_of_the median_link_of_node', 'Median of the connected median link length(\mum)'};

num_log_fit = size(log_log_fit_cell, 2);

lglg_selected_Q = wb_240_cube_stat_str.link_cap_stat.num_data.length > 50 & wb_240_cube_stat_str.link_cap_stat.num_data.length < 1000 &...
    wb_240_cube_stat_str.link_cap_stat.median.length > 35 & wb_240_cube_stat_str.mask_volume_density < 0.2 & ...
    wb_240_cube_stat_str.link_cap_stat.median.length < 120 & wb_240_cube_stat_str.capillary_length_density_m_mm3 > 0.3;
%% Plot linear regression in log-log plot
[loglog_fit_slope, loglog_fit_R] = deal(nan(num_log_fit));
loglog_plot_Q = false;
for iter_fix_feature = 1 : num_log_fit
    tmp_fix_data = log_log_fit_cell(:, iter_fix_feature);
    for iter_fit_feature = 1 : num_log_fit
        tmp_fit_data = log_log_fit_cell(:, iter_fit_feature);
        if iter_fix_feature == iter_fit_feature
            continue;
        end
        tmp_X = tmp_fix_data{1};
        tmp_Y = tmp_fit_data{1};
        tmp_is_valid_Q = ~isnan(tmp_X) & ~isnan(tmp_Y) & is_internal_240_cube_ratio == 1;
        num_bbox_before_selection = nnz(tmp_is_valid_Q);
        tmp_X = log10(tmp_X(lglg_selected_Q & tmp_is_valid_Q));
        tmp_Y = log10(tmp_Y(lglg_selected_Q & tmp_is_valid_Q));
        linear_fit_hdl = fitlm(tmp_X, tmp_Y);
        loglog_fit_slope(iter_fix_feature, iter_fit_feature) = linear_fit_hdl.Coefficients.Estimate(2);
        loglog_fit_R(iter_fix_feature, iter_fit_feature) = linear_fit_hdl.Rsquared.Ordinary;
        % Plot the fiture
        if loglog_plot_Q || (iter_fix_feature == 5 || iter_fit_feature == 5)
            fig_hdl = figure('Visible', 'on');
            fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
            ax_hdl = axes(fig_hdl); %#ok<LAXES>
            hist_hdl = histogram2(ax_hdl, tmp_X, tmp_Y, 'DisplayStyle', 'tile');
            hold(ax_hdl, 'on');
            fit_data_x = linspace(min(tmp_X, [], 'all'), max(tmp_X, [], 'all'), 30);
            fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
                'LineWidth', 2, 'Color', 'k');
            cbar_hdl = colorbar;
            cbar_hdl.Label.String = 'Number of data points';
            ax_hdl.XLabel.String = tmp_fix_data{3};
            ax_hdl.YLabel.String = tmp_fit_data{3};
            ax_hdl.FontSize = 14;
            ax_hdl.DataAspectRatio = [1,1,1];
            ax_hdl.FontWeight = 'bold';
            ax_hdl.ColorScale = 'log';
            ax_hdl.XScale = 'linear';
            ax_hdl.YScale = 'linear';
            ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%.1e'), 10 .^ (ax_hdl.XTick), 'UniformOutput', false);
            ax_hdl.YTickLabel = arrayfun(@(x) num2str(x, '%.1e'), 10 .^ (ax_hdl.YTick), 'UniformOutput', false);
            leg_hdl = legend(fit_plt_hdl, {sprintf('Slope: %.2e \\pm %.1e\nIntercept: %.2e \\pm %.1e\nR-squared: %.2f\nData size: %d (%d%%)', ...
                linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
                linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
                linear_fit_hdl.Rsquared.Ordinary, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/num_bbox_before_selection))});
            if ax_hdl.Position(3) < ax_hdl.Position(4)
                leg_hdl.Location = 'northoutside';
            else
                leg_hdl.Location = 'eastoutside';
            end
            
            fun_print_image_in_several_formats(fig_hdl, fullfile(save_im_folder, 'cube_length_feature_scaling', ...
                tmp_fix_data{2}, sprintf('%s_%s_int_cube_length_feature_scaling_%s_to_%s.png', ...
                dataset_name, stack, tmp_fix_data{2}, tmp_fit_data{2})));
            delete(fig_hdl);
        end
    end
end
stat_str.int_cap_length_scaling_fit.slope = loglog_fit_slope;
stat_str.int_cap_length_scaling_fit.rsquared = loglog_fit_R;
stat_str.int_cap_length_scaling_fit.feature_name = log_log_fit_cell(2, :);
stat_str.int_cap_feat_fit.cube_selected_Q = lglg_selected_Q;
%% Special case: fit median maximum tissue-capillary denstance vs (capillary length density)^{-1/2}
lglg_selected_Q = wb_240_cube_stat_str.link_cap_stat.num_data.length > 50 & wb_240_cube_stat_str.link_cap_stat.num_data.length < 1000 &...
    wb_240_cube_stat_str.link_cap_stat.median.length > 35 & wb_240_cube_stat_str.mask_volume_density < 0.05 & ...
    wb_240_cube_stat_str.link_cap_stat.median.length < 120 & ...
    wb_240_cube_stat_str.capillary_length_density_m_mm3 > 0.3 & wb_240_cube_stat_str.capillary_length_density_m_mm3 < 3.5 & ...
    wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max < 60;

tmp_plot_y = wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max;
tmp_plot_x = wb_240_cube_stat_str.capillary_length_density_m_mm3 .* 1e-3;
tmp_is_valid_Q = ~isnan(tmp_plot_x) & ~isnan(tmp_plot_y) & (wb_240_cube_stat_str.cube_in_brain_mask_ratio == 1);
tmp_selected_Q =  tmp_is_valid_Q & lglg_selected_Q;

tmp_plot_y = tmp_plot_y(tmp_selected_Q);
tmp_plot_x = tmp_plot_x(tmp_selected_Q);

tmp_plot_x = (tmp_plot_x) .^ (-1/2);
linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y);
% linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y, 'Intercept', false);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.4;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, tmp_plot_x, tmp_plot_y, 'DisplayStyle', 'tile');
ax_hdl.FontSize = 14;
ax_hdl.XLim(1) = 0;
ax_hdl.YLim(1) = 0;
hold(ax_hdl, 'on');
% fit_data_x = linspace(min(tmp_plot_x, [], 'all'), max(tmp_plot_x, [], 'all'), 30);
fit_data_x = linspace(0, max(tmp_plot_x, [], 'all'), 30);
lat_cubic_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 1.225, 'LineWidth', 3, 'LineStyle', '-.');
lat_8_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 1.138, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_b_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 1.015, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.951, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_c_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.868, 'LineWidth', 3, 'LineStyle', '-.');

fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
    'LineWidth', 3.5, 'Color', 'k');

leg_hdl_array = [fit_plt_hdl, lat_cubic_hdl, lat_8_3_a_hdl, lat_10_3_b_hdl, lat_10_3_a_hdl, lat_10_3_c_hdl];

leg_str_array = {sprintf('Slope: %.3f \\pm %.3f\nIntercept: %.2f \\pm %.2f\nR^2-Adjusted: %.2f\nData size: %d (%d%%)', ...
    linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
    linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
    linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/nnz(tmp_is_valid_Q))), ...
    'Cubic lattice: 1.225', '(8, 3)-a lattice: 1.138',  '(10, 3)-b lattice: 1.015', '(10, 3)-a lattice: 0.951', '(10, 3)-c lattice: 0.868'};
leg_hdl = legend(leg_hdl_array, leg_str_array, 'Location', 'northwest', 'Box', false);
leg_hdl.FontSize = 12;
cbar_hdl = colorbar;
cbar_hdl.Label.String = 'Number of data points';
cbar_hdl.FontSize = 12;
ax_hdl.XLabel.String = '$\rho_l^{-1/2} (\mu m)$';
ax_hdl.XLabel.FontSize = 18;
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$l_{max} (\mu m)$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.FontSize = 18;
ax_hdl.ColorScale = 'log';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_median_nearest_tissue_dt_max_vs_sqrt_cap_length_density.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Special case: fit median maximum tissue-capillary denstance vs (capillary length density)^{-1/2}
lglg_selected_Q = wb_240_cube_stat_str.link_cap_stat.num_data.length > 50 & wb_240_cube_stat_str.link_cap_stat.num_data.length < 1000 &...
    wb_240_cube_stat_str.link_cap_stat.median.length > 35 & wb_240_cube_stat_str.mask_volume_density < 0.05 & ...
    wb_240_cube_stat_str.link_cap_stat.median.length < 120 & ...
    wb_240_cube_stat_str.capillary_length_density_m_mm3 > 0.3 & wb_240_cube_stat_str.capillary_length_density_m_mm3 < 3.5 & ...
    wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max < 60 & wb_240_cube_stat_str.cap2vsl_vol_fraction > 0.75;

tmp_plot_y = wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_max;
tmp_plot_x = wb_240_cube_stat_str.capillary_length_density_m_mm3 .* 1e-3;
tmp_is_valid_Q = ~isnan(tmp_plot_x) & ~isnan(tmp_plot_y) & (wb_240_cube_stat_str.cube_in_brain_mask_ratio == 1);
tmp_selected_Q =  tmp_is_valid_Q & lglg_selected_Q;

tmp_plot_y = tmp_plot_y(tmp_selected_Q);
tmp_plot_x = tmp_plot_x(tmp_selected_Q);

tmp_plot_x = (tmp_plot_x) .^ (-1/2);
linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y);
% linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y, 'Intercept', false);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.4;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, tmp_plot_x, tmp_plot_y, 'DisplayStyle', 'tile');
ax_hdl.FontSize = 14;
ax_hdl.XLim(1) = 0;
% ax_hdl.YLim(1) = 0;
hold(ax_hdl, 'on');
% fit_data_x = linspace(min(tmp_plot_x, [], 'all'), max(tmp_plot_x, [], 'all'), 30);
fit_data_x = linspace(0, max(tmp_plot_x, [], 'all'), 30);
lat_cubic_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 1.225 + linear_fit_hdl.Coefficients.Estimate(1), 'LineWidth', 3, 'LineStyle', '-.');
lat_8_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 1.138 + linear_fit_hdl.Coefficients.Estimate(1), 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_b_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 1.015 + linear_fit_hdl.Coefficients.Estimate(1), 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.951 + linear_fit_hdl.Coefficients.Estimate(1), 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_c_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.868 + linear_fit_hdl.Coefficients.Estimate(1), 'LineWidth', 3, 'LineStyle', '-.');

fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
    'LineWidth', 3.5, 'Color', 'k');

leg_hdl_array = [fit_plt_hdl, lat_cubic_hdl, lat_8_3_a_hdl, lat_10_3_b_hdl, lat_10_3_a_hdl, lat_10_3_c_hdl];

leg_str_array = {sprintf('Slope: %.3f \\pm %.3f\nIntercept: %.2f \\pm %.2f\nR^2-Adjusted: %.2f\nData size: %d (%d%%)', ...
    linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
    linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
    linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/nnz(tmp_is_valid_Q))), ...
    'Cubic lattice: 1.225', '(8, 3)-a lattice: 1.138',  '(10, 3)-b lattice: 1.015', '(10, 3)-a lattice: 0.951', '(10, 3)-c lattice: 0.868'};

leg_hdl = legend(leg_hdl_array, leg_str_array, 'Location', 'northwest', 'Box', false);
leg_hdl.FontSize = 12;
cbar_hdl = colorbar;
cbar_hdl.Label.String = 'Number of data points';
cbar_hdl.FontSize = 12;
ax_hdl.XLabel.String = '$\rho_l^{-1/2} (\mu m)$';
ax_hdl.XLabel.FontSize = 18;
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$(l_{max} - r) (\mu m)$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.YLabel.FontSize = 18;
ax_hdl.ColorScale = 'log';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_median_nearest_tissue_dt_max_vs_sqrt_cap_length_density_cap2vsl_vol_gt_075.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Special case: Mean DT mean vs (capillary length density)^{-1/2}
tmp_plot_y = wb_240_cube_stat_str.link_cap_stat.mean.nearest_tissue_dt_mean;
tmp_plot_x = wb_240_cube_stat_str.capillary_length_density_m_mm3 .* 1e-3;
tmp_is_valid_Q = ~isnan(tmp_plot_x) & ~isnan(tmp_plot_y) & (wb_240_cube_stat_str.cube_in_brain_mask_ratio == 1);
tmp_selected_Q =  tmp_is_valid_Q & lglg_selected_Q;

tmp_plot_y = tmp_plot_y(tmp_selected_Q);
tmp_plot_x = tmp_plot_x(tmp_selected_Q);

tmp_plot_x = (tmp_plot_x) .^ (-1/2);
linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y);
% linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y, 'Intercept', false);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.4;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, tmp_plot_x, tmp_plot_y, 'DisplayStyle', 'tile');
ax_hdl.XLim(1) = 0;
ax_hdl.YLim(1) = 0;
hold(ax_hdl, 'on');
% % fit_data_x = linspace(min(tmp_plot_x, [], 'all'), max(tmp_plot_x, [], 'all'), 30);
fit_data_x = linspace(0, max(tmp_plot_x, [], 'all'), 30);
lat_cubic_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.482, 'LineWidth', 3, 'LineStyle', '-.');
lat_8_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.516, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_b_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.465, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.486, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_c_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.460, 'LineWidth', 3, 'LineStyle', '-.');
% 
fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
    'LineWidth', 3.5, 'Color', 'k');

leg_hdl_array = [fit_plt_hdl, lat_cubic_hdl, lat_8_3_a_hdl, lat_10_3_b_hdl, lat_10_3_a_hdl, lat_10_3_c_hdl];

leg_str_array = {sprintf('Slope: %.3f \\pm %.3f\nIntercept: %.2f \\pm %.2f\nR^2-Adjusted: %.2f\nData size: %d (%d%%)', ...
    linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
    linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
    linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/nnz(tmp_is_valid_Q))), ...
    'Cubic lattice: 0.482', '(8, 3)-a lattice: 0.516',  '(10, 3)-b lattice: 0.465', '(10, 3)-a lattice: 0.486', '(10, 3)-c lattice: 0.460'};

leg_hdl = legend(leg_hdl_array, leg_str_array, 'Location', 'northwest', 'Box', false);

cbar_hdl = colorbar;
cbar_hdl.Label.String = 'Number of data points';

ax_hdl.XLabel.String = '$\rho_l^{-1/2} (\mu m)$';
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$d_{avg} (\mu m)$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.XLabel.FontSize = 18;
ax_hdl.YLabel.FontSize = 18;
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
ax_hdl.ColorScale = 'log';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_mean_nearest_tissue_dt_mean_vs_isqrt_cap_length_density_cap2vsl_vol_gt_075.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%% Special case: Median DT mean vs (capillary length density)^{-1/2}
tmp_plot_y = wb_240_cube_stat_str.link_cap_stat.median.nearest_tissue_dt_mean;
tmp_plot_x = wb_240_cube_stat_str.capillary_length_density_m_mm3 .* 1e-3;
tmp_is_valid_Q = ~isnan(tmp_plot_x) & ~isnan(tmp_plot_y) & (wb_240_cube_stat_str.cube_in_brain_mask_ratio == 1);
tmp_selected_Q =  tmp_is_valid_Q & lglg_selected_Q;

tmp_plot_y = tmp_plot_y(tmp_selected_Q);
tmp_plot_x = tmp_plot_x(tmp_selected_Q);

tmp_plot_x = (tmp_plot_x) .^ (-1/2);
linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y);
% linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y, 'Intercept', false);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.4;
ax_hdl = axes(fig_hdl);
histogram2(ax_hdl, tmp_plot_x, tmp_plot_y, 'DisplayStyle', 'tile');
ax_hdl.XLim(1) = 0;
ax_hdl.YLim(1) = 0;
hold(ax_hdl, 'on');
% % fit_data_x = linspace(min(tmp_plot_x, [], 'all'), max(tmp_plot_x, [], 'all'), 30);
fit_data_x = linspace(0, max(tmp_plot_x, [], 'all'), 30);
lat_cubic_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.482, 'LineWidth', 3, 'LineStyle', '-.');
lat_8_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.497, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_b_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.464, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_a_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.484, 'LineWidth', 3, 'LineStyle', '-.');
lat_10_3_c_hdl = plot(ax_hdl, fit_data_x, fit_data_x .* 0.460, 'LineWidth', 3, 'LineStyle', '-.');
% 
fit_plt_hdl = plot(ax_hdl, fit_data_x, fit_data_x * linear_fit_hdl.Coefficients.Estimate(2) + linear_fit_hdl.Coefficients.Estimate(1), ...
    'LineWidth', 3.5, 'Color', 'k');

leg_hdl_array = [fit_plt_hdl, lat_cubic_hdl, lat_8_3_a_hdl, lat_10_3_b_hdl, lat_10_3_a_hdl, lat_10_3_c_hdl];

leg_str_array = {sprintf('Slope: %.3f \\pm %.3f\nIntercept: %.2f \\pm %.2f\nR^2-Adjusted: %.2f\nData size: %d (%d%%)', ...
    linear_fit_hdl.Coefficients.Estimate(2), linear_fit_hdl.Coefficients.SE(2),...
    linear_fit_hdl.Coefficients.Estimate(1), linear_fit_hdl.Coefficients.SE(1),...
    linear_fit_hdl.Rsquared.Adjusted, linear_fit_hdl.NumObservations, round(100 * linear_fit_hdl.NumObservations/nnz(tmp_is_valid_Q))), ...
    'Cubic lattice: 0.482', '(8, 3)-a lattice: 0.484',  '(10, 3)-b lattice: 0.464', '(10, 3)-a lattice: 0.484', '(10, 3)-c lattice: 0.460'};

leg_hdl = legend(leg_hdl_array, leg_str_array, 'Location', 'northwest', 'Box', false);

cbar_hdl = colorbar;
cbar_hdl.Label.String = 'Number of data points';

ax_hdl.XLabel.String = '$\rho_l^{-1/2} (\mu m)$';
ax_hdl.XLabel.Interpreter = 'latex';
ax_hdl.YLabel.String = '$d_{avg} (\mu m)$';
ax_hdl.YLabel.Interpreter = 'latex';
ax_hdl.XLabel.FontSize = 18;
ax_hdl.YLabel.FontSize = 18;
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
ax_hdl.ColorScale = 'log';
fig_fp = fullfile(save_im_folder, sprintf('%s_%s_median_nearest_tissue_dt_mean_vs_isqrt_cap_length_density_cap2vsl_vol_gt_075.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
%%
tmp_plot_y = wb_240_cube_stat_str.link_all_stat.num_data.length;
tmp_plot_x = wb_240_cube_stat_str.node_stat.num_data.degree;
tmp_is_valid_Q = ~isnan(tmp_plot_x) & ~isnan(tmp_plot_y) & (wb_240_cube_stat_str.cube_in_brain_mask_ratio == 1);
tmp_selected_Q =  tmp_is_valid_Q & lglg_selected_Q;

tmp_plot_y = tmp_plot_y(tmp_selected_Q);
tmp_plot_x = tmp_plot_x(tmp_selected_Q);
linear_fit_hdl = fitlm(tmp_plot_x, tmp_plot_y, 'Intercept', false);
%% 2D Histogram
%% Visualize the outliers
% plot_info_cell = cell(3, 0);
% plot_info_cell(:, end+1) = {double(wb_240_cube_stat_str.link_cap_stat.median.length < 40 & is_valid_Q), 'Capillary_median_length_outliers', 'is_outlier'};
% plot_info_cell{1}(isnan(wb_240_cube_stat_str.link_cap_stat.median.length)) = nan;
% fun_vis_generate_local_stat_video(plot_info_cell, grid_info, proj_im_str, false, 0, 100);
%
% implay(plot_info_cell{1}(isnan(wb_240_cube_stat_str.link_cap_stat.median.length)))
%% Visualize the image and reconstruction
% check_cube_grid_ind = find(plot_info_cell{1} == 1);
% check_list_ind = 1200;
% check_grid_sub = fun_ind2sub(grid_info.grid_size, check_cube_grid_ind(check_list_ind));
% check_im = DataManager.load_block_data(dataset_name, stack, grid_info.version, ...
%     check_grid_sub(1), check_grid_sub(2), check_grid_sub(3));
% check_mask = DataManager.load_block_mask(dataset_name, stack, mask_name, ...
%     check_grid_sub(1), check_grid_sub(2), check_grid_sub(3));
% check_mask = fun_reconstruct_block_mask(check_mask);
% DataManager.visualize_itksnap(check_im, check_mask);
%% Fraction of link without endpoints - internal 240-cube
is_not_nan_Q = ~isnan(wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q) & is_internal_240_cube_ratio == 1;
stat_str.int_cube_link_wo_ep_fraction.num_internal_cube = nnz(is_not_nan_Q);
stat_str.int_cube_link_wo_ep_fraction.ratio_not_ep = wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q(is_not_nan_Q);

stat_str.int_cube_link_wo_ep_fraction.num_cube_wo_ep = nnz(stat_str.int_cube_link_wo_ep_fraction.ratio_not_ep == 1);
stat_str.int_cube_link_wo_ep_fraction.fraction_cube_wo_ep = stat_str.int_cube_link_wo_ep_fraction.num_cube_wo_ep/ ...
    stat_str.int_cube_link_wo_ep_fraction.num_internal_cube;

tmp_fig = figure;
tmp_fig.Position(3:4) = tmp_fig.Position(3:4) .* 2;
ax_hdl = axes(tmp_fig);
hist_hdl = histogram(ax_hdl, stat_str.int_cube_link_wo_ep_fraction.ratio_not_ep, 'Normalization', 'cdf');
ax_hdl.YScale = 'log';
ax_hdl.XLabel.String = 'Fraction of vessels without endpoint';
ax_hdl.YLabel.String = 'CDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
grid(ax_hdl, 'on');
ax_hdl.Title.String = sprintf('Fraction of internal 240-cubes without endpoint: %.2f%%', ...
    100 * stat_str.int_cube_link_wo_ep_fraction.fraction_cube_wo_ep);
stat_str.int_cube_link_wo_ep_fraction.hist_internal_cube = fun_get_stat_from_histogram_hdl(hist_hdl);

fun_print_image_in_several_formats(tmp_fig, fullfile(save_im_folder, 'Fraction_of_links_wo_ep_in_internal_240_cube.png'));
%% Fraction of capillary without endpoints - internal 240-cube
is_not_nan_Q = ~isnan(wb_240_cube_stat_str.link_cap_stat.mean.has_no_ep_Q) & is_internal_240_cube_ratio == 1;
stat_str.int_cube_cap_wo_ep_fraction.num_internal_cube = nnz(is_not_nan_Q);
stat_str.int_cube_cap_wo_ep_fraction.ratio_not_ep = wb_240_cube_stat_str.link_cap_stat.mean.has_no_ep_Q(is_not_nan_Q);

stat_str.int_cube_cap_wo_ep_fraction.num_cube_wo_ep = nnz(stat_str.int_cube_cap_wo_ep_fraction.ratio_not_ep == 1);
stat_str.int_cube_cap_wo_ep_fraction.fraction_cube_wo_ep = stat_str.int_cube_cap_wo_ep_fraction.num_cube_wo_ep/ ...
    stat_str.int_cube_cap_wo_ep_fraction.num_internal_cube;
tmp_fig = figure;
tmp_fig.Position(3:4) = tmp_fig.Position(3:4) .* 2;
ax_hdl = axes(tmp_fig);
hist_hdl = histogram(ax_hdl, stat_str.int_cube_cap_wo_ep_fraction.ratio_not_ep, 'Normalization', 'cdf');
ax_hdl.YScale = 'log';
ax_hdl.XLabel.String = 'Fraction of vessels without endpoint';
ax_hdl.YLabel.String = 'CDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
grid(ax_hdl, 'on');
ax_hdl.Title.String = sprintf('Fraction of internal 240-cubes without endpoint: %.2f%%', ...
    100 * stat_str.int_cube_cap_wo_ep_fraction.fraction_cube_wo_ep);
stat_str.int_cube_cap_wo_ep_fraction.hist_internal_cube = fun_get_stat_from_histogram_hdl(hist_hdl);

fun_print_image_in_several_formats(tmp_fig, fullfile(save_im_folder, 'Fraction_of_cap_wo_ep_in_internal_240_cube.png'));
%% Fraction of capillaries without endpoints - all valid
is_not_nan_Q = ~isnan(wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q);
stat_str.cube_cap_wo_ep_fraction.num_internal_cube = nnz(is_not_nan_Q);
stat_str.cube_cap_wo_ep_fraction.ratio_not_ep = wb_240_cube_stat_str.link_all_stat.mean.has_no_ep_Q(is_not_nan_Q);
stat_str.cube_cap_wo_ep_fraction.num_cube_wo_ep = nnz(stat_str.cube_cap_wo_ep_fraction.ratio_not_ep == 1);
stat_str.cube_cap_wo_ep_fraction.fraction_cube_wo_ep = stat_str.cube_cap_wo_ep_fraction.num_cube_wo_ep/ ...
    stat_str.cube_cap_wo_ep_fraction.num_internal_cube;

tmp_fig = figure;
tmp_fig.Position(3:4) = tmp_fig.Position(3:4) .* 2;
ax_hdl = axes(tmp_fig);
hist_hdl = histogram(ax_hdl, stat_str.cube_cap_wo_ep_fraction.ratio_not_ep, 'Normalization', 'cdf');
stat_str.cube_cap_wo_ep_fraction.hist_cube = fun_get_stat_from_histogram_hdl(hist_hdl);

ax_hdl.YScale = 'log';
ax_hdl.XLabel.String = 'Fraction of capillaries without endpoint';
ax_hdl.YLabel.String = 'CDF';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
grid(ax_hdl, 'on');
ax_hdl.Title.String = sprintf('Fraction of 240-cubes without endpoint: %.2f%%', ...
    100 * stat_str.cube_cap_wo_ep_fraction.fraction_cube_wo_ep);

fun_print_image_in_several_formats(tmp_fig, fullfile(save_im_folder, 'Fraction_of_link_wo_ep_in_240_cube.png'));
%% Anisotropy FA, PCV1 and p-Value
% All vessels - volume weighted
ai_vessel_stat = struct;
ai_vessel_stat.fa_p = wb_240_cube_stat_str.wb_ai_all_vw.fa_p(is_internal_240_cube_Q);
ai_vessel_stat.fa = wb_240_cube_stat_str.wb_ai_all_vw.fractional_anisotropy(is_internal_240_cube_Q);
ai_vessel_stat.svd_rmax = wb_240_cube_stat_str.wb_ai_all_vw.svd_value_ratio(is_internal_240_cube_Q, 1);
ai_vessel_stat.fa_z = wb_240_cube_stat_str.wb_ai_all_vw.fa_z(is_internal_240_cube_Q);
ai_vessel_stat.svd_1_p = wb_240_cube_stat_str.wb_ai_all_vw.svd_1_p(is_internal_240_cube_Q);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1] .* 2;
ax_hdl_1 = subplot(1,3,1);
tmp_x = ai_vessel_stat.fa;
tmp_x = tmp_x(~isnan(tmp_x));
histogram(ax_hdl_1, tmp_x, 'Normalization', 'pdf');
ax_hdl_1.XLabel.String = 'FA';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
box(ax_hdl_1, 'off');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)), 'Location', 'best');

ax_hdl_2 = subplot(1,3,2);
tmp_x = ai_vessel_stat.svd_rmax;
tmp_x = tmp_x(~isnan(tmp_x));
histogram(ax_hdl_2, tmp_x, 'Normalization', 'pdf');
ax_hdl_2.XLabel.String = 'Normalized PCV';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)), 'Location', 'best');

ax_hdl_3 = subplot(1,3,3);
tmp_x = ai_vessel_stat.fa_p;
tmp_x = tmp_x(~isnan(tmp_x));
tmp_x_bin = [0, 10.^(-4 : 0.5 : 0)];
histogram(ax_hdl_3, tmp_x, tmp_x_bin, 'Normalization', 'cdf');
ax_hdl_3.XScale = 'log';
ax_hdl_3.YLabel.String = 'CDF';
ax_hdl_3.XLabel.String = 'FA p-value';
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.FontSize = 14;
box(ax_hdl_3, 'off');
legend(ax_hdl_3, sprintf('# cubes: %d\np < 0.05: %d%%\np < 0.01: %d%%\np < 10^{-3}: %d%%\np < 10^{-4}: %d%%', ...
    numel(ai_vessel_stat.fa_p), round(100*nnz(ai_vessel_stat.fa_p < 0.05)/numel(ai_vessel_stat.fa_p)), ...
    round(100*nnz(ai_vessel_stat.fa_p < 0.01)/numel(ai_vessel_stat.fa_p)), ...
    round(100*nnz(ai_vessel_stat.fa_p < 0.001)/numel(ai_vessel_stat.fa_p)), ...
    round(100*nnz(ai_vessel_stat.fa_p < 0.0001) / numel(ai_vessel_stat.fa_p))), 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'on');

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_internal_cubes_volume_weighted_vessel_anisotropy_stat.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

stat_str.ai_all_vw = ai_vessel_stat;
%% Capillaries - volume weighted
ai_cap_vw_stat = struct;
ai_cap_vw_stat.fa_p = wb_240_cube_stat_str.wb_ai_cap_vw.fa_p(is_internal_240_cube_Q);
ai_cap_vw_stat.fa = wb_240_cube_stat_str.wb_ai_cap_vw.fractional_anisotropy(is_internal_240_cube_Q);
ai_cap_vw_stat.svd_rmax = wb_240_cube_stat_str.wb_ai_cap_vw.svd_value_ratio(is_internal_240_cube_Q, 1);
ai_cap_vw_stat.fa_z = wb_240_cube_stat_str.wb_ai_cap_vw.fa_z(is_internal_240_cube_Q);
ai_cap_vw_stat.svd_1_p = wb_240_cube_stat_str.wb_ai_cap_vw.svd_1_p(is_internal_240_cube_Q);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1] .* 2;
ax_hdl_1 = subplot(1,3,1);
tmp_x = ai_cap_vw_stat.fa;
tmp_x = tmp_x(~isnan(tmp_x));
histogram(ax_hdl_1, tmp_x, 'Normalization', 'pdf');
ax_hdl_1.XLabel.String = 'FA';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
box(ax_hdl_1, 'off');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)), 'Location', 'best');

ax_hdl_2 = subplot(1,3,2);
tmp_x = ai_cap_vw_stat.svd_rmax;
tmp_x = tmp_x(~isnan(tmp_x));
histogram(ax_hdl_2, tmp_x, 'Normalization', 'pdf');
ax_hdl_2.XLabel.String = 'Normalized PCV';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)), 'Location', 'best');

ax_hdl_3 = subplot(1,3,3);
tmp_x = ai_cap_vw_stat.fa_p;
tmp_x = tmp_x(~isnan(tmp_x));
tmp_x_bin = [0, 10.^(-4 : 0.5 : 0)];
histogram(ax_hdl_3, tmp_x, tmp_x_bin, 'Normalization', 'cdf');
ax_hdl_3.XScale = 'log';
ax_hdl_3.YLabel.String = 'CDF';
ax_hdl_3.XLabel.String = 'FA p-value';
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.FontSize = 14;
box(ax_hdl_3, 'off');
legend(ax_hdl_3, sprintf('# cubes: %d\np < 0.05: %d%%\np < 0.01: %d%%\np < 10^{-3}: %d%%\np < 10^{-4}: %d%%', ...
    numel(ai_cap_vw_stat.fa_p), round(100*nnz(ai_cap_vw_stat.fa_p < 0.05)/numel(ai_cap_vw_stat.fa_p)), ...
    round(100*nnz(ai_cap_vw_stat.fa_p < 0.01)/numel(ai_cap_vw_stat.fa_p)), ...
    round(100*nnz(ai_cap_vw_stat.fa_p < 0.001)/numel(ai_cap_vw_stat.fa_p)), ...
    round(100*nnz(ai_cap_vw_stat.fa_p < 0.0001) / numel(ai_cap_vw_stat.fa_p))), 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'on');

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_internal_cubes_volume_weighted_capillary_anisotropy_stat.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

stat_str.ai_cap_vw = ai_cap_vw_stat;
%% Capillaries - length weight
ai_cap_lw_stat = struct;
ai_cap_lw_stat.fa_p = wb_240_cube_stat_str.wb_ai_cap_lw.fa_p(is_internal_240_cube_Q);
ai_cap_lw_stat.fa = wb_240_cube_stat_str.wb_ai_cap_lw.fractional_anisotropy(is_internal_240_cube_Q);
ai_cap_lw_stat.svd_rmax = wb_240_cube_stat_str.wb_ai_cap_lw.svd_value_ratio(is_internal_240_cube_Q, 1);
ai_cap_lw_stat.fa_z = wb_240_cube_stat_str.wb_ai_cap_lw.fa_z(is_internal_240_cube_Q);
ai_cap_lw_stat.svd_1_p = wb_240_cube_stat_str.wb_ai_cap_lw.svd_1_p(is_internal_240_cube_Q);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3,1] .* 2;
ax_hdl_1 = subplot(1,3,1);
tmp_x = ai_cap_lw_stat.fa;
tmp_x = tmp_x(~isnan(tmp_x));
histogram(ax_hdl_1, tmp_x, 'Normalization', 'pdf');
ax_hdl_1.XLabel.String = 'FA';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_1.FontSize = 14;
box(ax_hdl_1, 'off');
legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)), 'Location', 'best');

ax_hdl_2 = subplot(1,3,2);
tmp_x = ai_cap_lw_stat.svd_rmax;
tmp_x = tmp_x(~isnan(tmp_x));
histogram(ax_hdl_2, tmp_x, 'Normalization', 'pdf');
ax_hdl_2.XLabel.String = 'Normalized PCV';
ax_hdl_2.YLabel.String = 'PDF';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_x)), 'Location', 'best');

ax_hdl_3 = subplot(1,3,3);
tmp_x = ai_cap_lw_stat.fa_p;
tmp_x = tmp_x(~isnan(tmp_x));
tmp_x_bin = [0, 10.^(-4 : 0.5 : 0)];
histogram(ax_hdl_3, tmp_x, tmp_x_bin, 'Normalization', 'cdf');
ax_hdl_3.XScale = 'log';
ax_hdl_3.YLabel.String = 'CDF';
ax_hdl_3.XLabel.String = 'FA p-value';
ax_hdl_3.YLim = [0, 1];
ax_hdl_3.FontSize = 14;
box(ax_hdl_3, 'off');
legend(ax_hdl_3, sprintf('# cubes: %d\np < 0.05: %d%%\np < 0.01: %d%%\np < 10^{-3}: %d%%\np < 10^{-4}: %d%%', ...
    numel(ai_cap_lw_stat.fa_p), round(100*nnz(ai_cap_lw_stat.fa_p < 0.05)/numel(ai_cap_lw_stat.fa_p)), ...
    round(100*nnz(ai_cap_lw_stat.fa_p < 0.01)/numel(ai_cap_lw_stat.fa_p)), ...
    round(100*nnz(ai_cap_lw_stat.fa_p < 0.001)/numel(ai_cap_lw_stat.fa_p)), ...
    round(100*nnz(ai_cap_lw_stat.fa_p < 0.0001) / numel(ai_cap_lw_stat.fa_p))), 'Location', 'northwest', 'Interpreter', 'tex', 'Box', 'on');

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_internal_cubes_length_weighted_capillary_anisotropy_stat.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

stat_str.ai_cap_lw = ai_cap_lw_stat;
%% Anisotropy capillary volume fraction vs p-Value joint histogram
% Bin p-value by capillary volume fraction
tmp_x_bin_edge = 0 : 0.25 : 1;
tmp_y_bin_edge = [0, 10 .^ (-4 : 1 : 0)];

[fig_hdl, ax_hdl, cap_ai_vr_stat] = fun_vis_stack_histogram(wb_240_cube_stat_str.cap2vsl_vol_fraction(stat_str.cube_volume_ratio_in_brain_mask == 1), ...
    wb_240_cube_stat_str.wb_ai_cap_vw.fa_p(stat_str.cube_volume_ratio_in_brain_mask == 1), tmp_x_bin_edge, tmp_y_bin_edge);
leg_hdl = legend(ax_hdl, '[0, 10^{-4})', '[10^{-4}, 10^{-3})', '[10^{-3}, 10^{-2})', '[10^{-2}, 10^{-1})', '[10^{-1}, 1]', 'Location', 'northwest');
leg_hdl.Title.String = 'p-Value';
ax_hdl.XLabel.String = 'Capillary - vessel volume ratio';
ax_hdl.YLabel.String = 'Probability';
ax_hdl.Title.String = 'Volume-weighted capillary FA test p-Value';

fig_fp = fullfile(save_im_folder, sprintf('%s_%s_internal_cubes_volume_weighted_cap_FA_test_pValue_vs_capvsl_vol_ratio_hist.png', dataset_name, stack));
fun_print_image_in_several_formats(fig_hdl, fig_fp);

stat_str.ai_cap_vw_fap_vs_cap2vsl_vol_r = cap_ai_vr_stat;
%% Save data
if ~isfile(stat_str.info.filepath)
    save(stat_str.info.filepath, '-struct', 'stat_str');
else
    warning('The file already exist. Overwrite?')
end
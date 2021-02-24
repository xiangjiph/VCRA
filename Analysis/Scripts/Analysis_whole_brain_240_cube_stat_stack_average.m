%% Document
% This script plot PDF of local network features for each stack, as well as
% the averaged PDF.
% This script need to be run after
% Analysis_preprocess_whole_brain_local_statistics
%%
clc;clear;
DataManager = FileManager;
dataset_name = 'WholeBrain';
image_grid_version = '240_cube';
reconstruction_version = '240_cube_recon_sc';
skel_version = '240_cube_re';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
num_stack = numel(stack_list);
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_nonrigid.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');
merge_stack_name = 'all_stack';
saved_stat_folder_name = 'whole_brain_stat_sc';
wb_stat_folder_name = 'wb_internal';
stat_data_name = 'internal_cubes';
vis_ptl_min = 5e-3;
vis_ptl_max = 1 - vis_ptl_min;
%%  Load registration and brain mask
[wb_data_cube_stat_cell, wb_data_cube_pdf_str] = deal(cell(num_stack, 1));
tic_load_data = tic;
for iter_stack = 1 : num_stack
    wb_data_cube_stat_cell{iter_stack} = DataManager.load_analysis_data(dataset_name,...
        stack_list{iter_stack}, sprintf('%s_%s_%s_240_cube_stat_data.mat',...
        dataset_name, stack_list{iter_stack}, reconstruction_version), ...
        saved_stat_folder_name);
    wb_data_cube_pdf_str{iter_stack} = DataManager.load_analysis_data(dataset_name, ...
        stack_list{iter_stack}, sprintf('%s_%s_%s_240_cube_stat_pdf.mat', ...
        dataset_name, stack_list{iter_stack}, reconstruction_version), ...
        saved_stat_folder_name);
end
fprintf('Finish loading data. Elapsed time is %f seconds.\n', toc(tic_load_data));
%% Initialize structure for recording plotting information
stat_str = struct;
stat_str.internal_cube_min_in_mask_vol_f = 1;
stat_str.minimal_cap2vsl_vol_ratio = 0;
stat_str.info.folder = fullfile(DataManager.fp_visualization_folder(dataset_name, merge_stack_name), sprintf('%s_min_capVr%d', ...
    wb_stat_folder_name, stat_str.minimal_cap2vsl_vol_ratio * 100));

stat_str.info.dataset_name = dataset_name;
stat_str.info.stack_list = stack_list;
stat_str.info.process_date = datestr(now);
stat_str.info.filepath = fullfile(stat_str.info.folder,...
    sprintf('%s_%s_%s_plot_data.mat', dataset_name, merge_stack_name, wb_stat_folder_name));
%% Add features
for iter_stack = 1 : num_stack
    wb_data_cube_stat_cell{iter_stack}.node_density_mm3 = wb_data_cube_stat_cell{iter_stack}.node_stat.num_data.degree ./ (0.24 ^ 3);
    wb_data_cube_stat_cell{iter_stack}.link_all_stat.mean.has_ep_Q = 1 - wb_data_cube_stat_cell{iter_stack}.link_all_stat.mean.has_no_ep_Q;
    wb_data_cube_stat_cell{iter_stack}.link_cap_stat.mean.has_ep_Q = 1 - wb_data_cube_stat_cell{iter_stack}.link_cap_stat.mean.has_no_ep_Q;
    wb_data_cube_stat_cell{iter_stack}.noncap2vsl_vol_fraction = 1 - wb_data_cube_stat_cell{iter_stack}.cap2vsl_vol_fraction;
end
%% Histograms - whole brain statistics - scalar data
fprintf('Generating plot setting cell array for scalar data\n');
plot_info_cell = cell(3, 0);
% Anisotropy
plot_info_cell(:, end+1) = {'wb_ai_all_vw.fractional_anisotropy', 'Vessel_volume-weighted_fractional_anisotropy', 'Vessel fractioanal anisotropy'};
plot_info_cell(:, end+1) = {'wb_ai_all_vw.fa_z', 'Vessel_volume-weighted_FA_z-score', 'Vessel FA_z'};
plot_info_cell(:, end+1) = {'wb_ai_all_vw.fa_p', 'Vessel_volume-weighted_FA_p-Value', 'Vessel FA p-Value'};
plot_info_cell(:, end+1) = {'wb_ai_all_vw.svd_value_ratio', 'Vessel_volume-weighted_PCV1', 'Vessel PCV1'};

plot_info_cell(:, end+1) = {'wb_ai_cap_vw.fractional_anisotropy', 'Capillary_volume-weighted_fractional_anisotropy', 'Capillary fractioanal anisotropy'};
plot_info_cell(:, end+1) = {'wb_ai_cap_vw.fa_z', 'Capillary_volume-weighted_FA_z-score', 'Capillary FA_z'};
plot_info_cell(:, end+1) = {'wb_ai_cap_vw.fa_p', 'Capillary_volume-weighted_FA_p-Value', 'Capillary FA p-Value'};
plot_info_cell(:, end+1) = {'wb_ai_cap_vw.svd_value_ratio', 'Capillary_volume-weighted_PCV1', 'Capillary PCV1'};
% Capillary - vessel ratio
plot_info_cell(:, end+1) = {'cap2vsl_length_fraction', 'Capillary-vessel_length_ratio', 'Capillary-vessel length ratio'};
plot_info_cell(:, end+1) = {'cap2vsl_vol_fraction', 'Capillary-vessel_volume_ratio', 'Capillary-vessel volume ratio'};
plot_info_cell(:, end+1) = {'noncap2vsl_vol_fraction', 'Noncapillary-vessel_volume_ratio', 'Noncapillary-vessel volume ratio'};
plot_info_cell(:, end+1) = {'cap2vsl_surf_area_fraction', 'Capillary-vessel_surface_area_ratio', 'Capillary-vessel surface area ratio'};

% Global statistics
plot_info_cell(:, end+1) = {'mask_volume_density', 'Reconstructed_vessel_volume_density', 'Vessel volume density'};
plot_info_cell(:, end+1) = {'mask_surface_area_density_mm2_mm3', 'Reconstructed_vessel_surface_area', 'Vessel surface area density (mm^2/mm^3)'};
plot_info_cell(:, end+1) = {'link_length_density_m_mm3', 'Vessel_length_density', 'Vessel length density (m/mm^3)'};

plot_info_cell(:, end+1) = {'capillary_volume_density', 'Capillary_volume_density', 'Capillary volume density'};
plot_info_cell(:, end+1) = {'capillary_length_density_m_mm3', 'Capillary_length_density', 'Capillary length density (m/mm^3)'};
plot_info_cell(:, end+1) = {'capillary_surface_area_density_mm2_mm3', 'Capillary_surface_area_density', 'Capillary surface area density (mm^2/mm^3)'};

% capillary segment statistics
plot_info_cell(:, end+1) = {'link_cap_stat.median.length', 'Median_capillary_segment_length', 'Median capillary segment length (\mum)'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.length', 'Average_capillary_segment_length', 'Average capillary segment length (\mum)'};
plot_info_cell(:, end+1) = {'link_cap_stat.median.nearest_tissue_dt_max', 'Median_maximum_distance_between_tissue_and_the_nearest_capillary', 'Maximum capillary-tissue distance (\mum)'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.nearest_tissue_radius', 'Mean_capillary_nearest_tissue_radius', 'Mean capillary tissue radius (\mum)'};
plot_info_cell(:, end+1) = {'link_cap_stat.median.nearest_tissue_radius', 'Median_capillary_nearest_tissue_radius', 'Median capillary tissue radius (\mum)'};
plot_info_cell(:, end+1) = {'link_cap_stat.num_data.length', 'Number_of_capillaries', 'Number of capillary'};
plot_info_cell(:, end+1) = {'link_cap_stat.median.nearest_tissue_volume', 'Median_capillary_nearest_tissue_volume', 'Median capillary nearest tissue volume (\mum^3)'};
% Loop
plot_info_cell(:, end+1) = {'link_cap_stat.mean.shortest_loop_geodesic_length', 'Average_number_of_edges_in_the_shortest_loop_of_capillary', 'Average number of edges in shortest loop'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.shortest_loop_length', 'Average_length_of_the_shortest_loop_of_capillary', 'Average length of shortest loop (\mum)'};
plot_info_cell(:, end+1) = {'link_cap_stat.std.shortest_loop_geodesic_length', 'STD_number_of_edges_in_the_shortest_loop_of_capillary', 'STD of the number of edges in shortest loop'};
% Branching order
plot_info_cell(:, end+1) = {'link_cap_stat.mean.capillary_branching_order', 'Mean_capillary_branch_order', 'Capillary branch order'};
% Distance to large vessels
plot_info_cell(:, end+1) = {'link_cap_stat.median.dist_to_nearest_noncapillary_mean', 'Median_distance_to_nearest_noncapillary', 'Median distance to nearest noncapillary(\mum)'};
% Geometry
plot_info_cell(:, end+1) = {'link_cap_stat.mean.straightness', 'Mean_capillary_straightness', 'Average capillary segment straightness'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.tortuosity', 'Mean_capillary_tortuosity', 'Average capillary segment tortuosity'};
plot_info_cell(:, end+1) = {'link_all_stat.mean.tortuosity', 'Mean_vessel_tortuosity', 'Average vessel segment tortuosity'};
% Noncapillary
plot_info_cell(:, end+1) = {'link_all_stat.mean.noncapillary_nearest_tissue_radius', 'Mean_noncapillary_nearest_tissue_radius' 'Average noncapillary tissue radius(\mum)'};
plot_info_cell(:, end+1) = {'link_all_stat.mean.noncapillary_nearest_capillary_num_vxl', 'Mean_nearest_caillary_total_num_skel_vxl', 'Average number of nearest capillary skeleton voxels'};

% Node statistics
plot_info_cell(:, end+1) = {'node_stat.mean.degree', 'Average_node_degree', 'Average node degree'};
plot_info_cell(:, end+1) = {'node_stat.median.nearest_node_dist', 'Median_distance_to_the_nearest_node', 'Median distance to the nearest node(\mum)'};
plot_info_cell(:, end+1) = {'node_stat.mean.path_to_nearest_neighbor_geodesic_length', 'Mean_geodesic_distance_to_nearest_node', 'Average number of edges to the nearest node'};
plot_info_cell(:, end+1) = {'node_density_mm3', 'Node_density', 'Node density(mm^{-3})'};
plot_info_cell(:, end+1) = {'node_stat.median.link_length_min', 'Median_length_of_the_shortest_link_of_node', 'Median of the connected shortest links length (\mum)'};
plot_info_cell(:, end+1) = {'node_stat.median.link_length_max', 'Median_length_of_the_longest_link_of_node', 'Median of the connected longest link length (\mum)'};
plot_info_cell(:, end+1) = {'node_stat.median.link_length_median', 'Median_length_of_the median_link_of_node', 'Median of the connected median link length (\mum)'};

plot_info_cell(:, end+1) = {'node_stat.mean.link_ori_ep_sin_elevation', 'Mean_vsl_ori_ep_sin_ele', 'Mean link tangent sin(\theta_{ele})'};
plot_info_cell(:, end+1) = {'node_stat.mean.link_ep2ep_sin_elevation', 'Mean_vessel_ep2ep_vec_sin_ele', 'Mean link ep2ep vector sin(\theta_{ele})'};

plot_info_cell(:, end+1) = {'node_stat.mean.link_ori_ep_cos_elevation', 'Mean_vsl_ori_ep_cos_ele', 'Mean link tangent cos(\theta_{ele})'};
plot_info_cell(:, end+1) = {'node_stat.mean.link_ep2ep_cos_elevation', 'Mean_vessel_ep2ep_vec_cos_ele', 'Mean link ep2ep vector cos(\theta_{ele})'};
% Connected node distance rank
plot_info_cell(:, end+1) = {'node_stat.mean.conn_node_drk_med', 'Average_median_conn_node_distance_rank', 'Average connected node median distance rank'};
plot_info_cell(:, end+1) = {'node_stat.mean.conn_node_drk_min', 'Average_of_min_conn_node_distance_rank' 'Average connected node minimum distance rank'};
plot_info_cell(:, end+1) = {'node_stat.mean.conn_node_drk_max', 'Average_of_max_conn_node_distance_rank', 'Average connected node maximum distance rank'};
% Link branching angle - endpoint tengent
plot_info_cell(:, end+1) = {'node_stat.median.link_ori_agl_min', 'Median_of_min_link_ori_agl', 'Median of minimum branch tangent angle (^{\circ})'};
plot_info_cell(:, end+1) = {'node_stat.median.link_ori_agl_med', 'Median_of_med_link_ori_agl', 'Median of median branch tangent angle (^{\circ})'};
plot_info_cell(:, end+1) = {'node_stat.median.link_ori_agl_max', 'Median_of_max_link_ori_agl', 'Median of maximum branch tangent angle (^{\circ})'};
% Link branching angle - endpoint to endpoint
plot_info_cell(:, end+1) = {'node_stat.median.link_ep2ep_agl_min', 'Median_of_mim_link_ep2ep_agl', 'Median of minimum branch angle (^{\circ})'};
plot_info_cell(:, end+1) = {'node_stat.median.link_ep2ep_agl_med', 'Median_of_med_link_ep2ep_agl', 'Median of median branch angle (^{\circ})'};
plot_info_cell(:, end+1) = {'node_stat.median.link_ep2ep_agl_max', 'Median_of_max_link_ep2ep_agl', 'Median of maximum branch angle (^{\circ})'};

% Segmentation quality
plot_info_cell(:, end+1) = {'link_all_stat.mean.has_no_ep_Q', 'Fraction_of_vessel_without_endpoint', 'Local fraction of links without endpoint'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.has_no_ep_Q', 'Fraction_of_capillary_without_endpoint', 'Local fraction of links without endpoint'};

plot_info_cell(:, end+1) = {'link_all_stat.mean.has_ep_Q', 'Fraction_of_links_with_ep', 'Local fraction of links with endpoint'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.has_ep_Q', 'Fraction_of_cap_with_endpoint', 'Local fraction of capillaries with endpoint'};

stat_str.feature_scalar_plot_setting = plot_info_cell;
num_plot_scalar_features = size(plot_info_cell, 2);
%% Generate histograms of 240-cube properties
% 1. For each brain
% 2. Average of brains.
for iter_feature = 1 : num_plot_scalar_features
    tmp_plot_info = stat_str.feature_scalar_plot_setting(:, iter_feature);
    tmp_vis_plt_max = vis_ptl_max;
    tmp_vis_plt_min = vis_ptl_min;
    
    stack_pdf_fit = struct;
    stack_pdf_fit.stack_list = stack_list;
    [stack_pdf_fit.hist_edge, stack_pdf_fit.interpolate_fun] = deal(cell(num_stack, 1));
    [stack_pdf_fit.mean, stack_pdf_fit.std] = deal(nan(num_stack, 1));
    
    stack_pdf_fit.data_cell = cell(num_stack, 1);
    stack_pdf_fit.feature_string_cell = cell(num_stack + 1, 1);
    
    for iter_stack = 1 : num_stack
        tmp_wb_data = wb_data_cube_stat_cell{iter_stack};
        tmp_is_internal_cube_Q = (tmp_wb_data.cube_in_brain_mask_ratio == stat_str.internal_cube_min_in_mask_vol_f) & ...
            (tmp_wb_data.cap2vsl_vol_fraction >= stat_str.minimal_cap2vsl_vol_ratio);
        tmp_hist_data = fun_getfield(tmp_wb_data, tmp_plot_info{1});
        if size(tmp_hist_data, 2) > 1
            tmp_hist_data = tmp_hist_data(:, 1);
        end
        tmp_hist_data = tmp_hist_data(tmp_is_internal_cube_Q & ~isnan(tmp_hist_data));
        if all(tmp_hist_data <= 1)
            tmp_vis_plt_max = 1;
            if all(tmp_hist_data >= 0)
                tmp_vis_plt_min = 0;
            end
        end
        
        tmp_stat_str = fun_analysis_get_basic_statistics(tmp_hist_data);
        
        stack_pdf_fit.data_cell{iter_stack} = tmp_stat_str;
        stack_pdf_fit.mean(iter_stack) = tmp_stat_str.mean;
        stack_pdf_fit.std(iter_stack) = tmp_stat_str.std;
        
        tmp_stack_string = sprintf('%s\n%.2e \\pm %.1e', strrep(stack_list{iter_stack}, ...
            '_', ''), tmp_stat_str.mean, tmp_stat_str.std);
        
        stack_pdf_fit.feature_string_cell{iter_stack} = tmp_stack_string;
        stack_pdf_fit.hist_edge{iter_stack} = tmp_stat_str.hist_edge;
        stack_pdf_fit.interpolate_fun{iter_stack} = griddedInterpolant(tmp_stat_str.hist_edge, ...
            [0, tmp_stat_str.hist_cdf], 'linear', 'nearest');
    end
    stack_pdf_fit.feature_name = tmp_plot_info{2};
    stack_pdf_fit.filepath = fullfile(stat_str.info.folder, 'feature_stat_hist', ...
        sprintf('%s_%s_%s_%s_pdf_fit_data.mat', dataset_name, merge_stack_name, stat_data_name, ...
        strrep(stack_pdf_fit.feature_name, ' ' , '_')));
    
    % Determine interval size
    num_bin = cellfun(@numel, stack_pdf_fit.hist_edge) - 1;
    edge_min = cellfun(@min, stack_pdf_fit.hist_edge);
    edge_max = cellfun(@max, stack_pdf_fit.hist_edge);
    stack_pdf_fit.int_x = linspace(min(edge_min), max(edge_max), min(num_bin));
    
    stack_pdf_fit.cdf_int_val = cellfun(@(x) x(stack_pdf_fit.int_x), ...
        stack_pdf_fit.interpolate_fun, 'UniformOutput', false);
    
    % Find the common plot bin range
    stack_pdf_fit.selected_int_x_Q = cellfun(@(x) (x >= tmp_vis_plt_min) & (x <= tmp_vis_plt_max), ...
        stack_pdf_fit.cdf_int_val, 'UniformOutput', false);
    stack_pdf_fit.selected_int_x_Q = all(cat(1, stack_pdf_fit.selected_int_x_Q{:}), 1);
    stack_pdf_fit.selected_int_x_ind = find(stack_pdf_fit.selected_int_x_Q) - 1;
    stack_pdf_fit.selected_int_x_ind = stack_pdf_fit.selected_int_x_ind(stack_pdf_fit.selected_int_x_ind > 0);
    
    stack_pdf_fit.pdf_for_avg = cell(num_stack, 1);
    fig_hdl = figure('Visible', 'off');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
    ax_hdl = axes(fig_hdl);
    hold(ax_hdl, 'on');
    for iter_stack = 1 : num_stack
        tmp_bin_edge = stack_pdf_fit.int_x;
        tmp_bin_width = diff(tmp_bin_edge);
        tmp_bin_val = stack_pdf_fit.cdf_int_val{iter_stack};
        tmp_bin_val = diff(tmp_bin_val) ./ tmp_bin_width;
        
        stack_pdf_fit.pdf_for_avg{iter_stack} = tmp_bin_val(stack_pdf_fit.selected_int_x_ind);
        
        [tmp_bin_val, tmp_bin_edge] = fun_analysis_select_histcount_edge_by_percentile(...
            tmp_bin_val, tmp_bin_edge, [tmp_vis_plt_min, tmp_vis_plt_max]);
        scatter(ax_hdl, movmean(tmp_bin_edge, 2, 'Endpoints', 'discard'), ...
            tmp_bin_val, 'o', 'filled');
        %         histogram(ax_hdl, 'BinCounts', tmp_bin_val, 'BinEdges', tmp_bin_edge);
    end
    ax_hdl.XLabel.String = tmp_plot_info{3};
    ax_hdl.YLabel.String = 'PDF';
    ax_hdl.FontSize = 14;
    ax_hdl.FontWeight = 'bold';
    grid(ax_hdl, 'off');
    box(ax_hdl, 'off');
    
    % Add average curve
    stack_pdf_fit.pdf_avg = cat(1, stack_pdf_fit.pdf_for_avg{:});
    stack_pdf_fit.pdf_avg = mean(stack_pdf_fit.pdf_avg, 1);
    stack_pdf_fit.pdf_std = cat(1, stack_pdf_fit.pdf_for_avg{:});
    stack_pdf_fit.pdf_std = std(stack_pdf_fit.pdf_std, 1);
    
    stack_pdf_fit.pdf_avg_bin_val = movmean(stack_pdf_fit.int_x, 2, 'Endpoints', 'discard');
    stack_pdf_fit.pdf_avg_bin_val = stack_pdf_fit.pdf_avg_bin_val(stack_pdf_fit.selected_int_x_ind);
    
    % Plot average PDF with shaded area PDF
    [ax_hdl, line_hdl, error_area_hdl]= fun_vis_errorbar_shaded(stack_pdf_fit.pdf_avg_bin_val, ...
        stack_pdf_fit.pdf_avg, stack_pdf_fit.pdf_std, ax_hdl);
    stack_pdf_fit.feature_string_cell{end} = sprintf('Average\n%.2e \\pm %.1e', ...
        mean(stack_pdf_fit.mean), mean(stack_pdf_fit.std));
    legend(ax_hdl, stack_pdf_fit.feature_string_cell{:}, 'Location', 'best');
    
    tmp_fp = strrep(stack_pdf_fit.filepath, '_fit_data.mat', '.png');
    
    fun_print_image_in_several_formats(fig_hdl, tmp_fp);
    save(stack_pdf_fit.filepath, 'stack_pdf_fit');
    delete(fig_hdl);
end
%% CDF for quantifying graph refinement quality
plot_type = 'CDF';
plot_info_cell = cell(5, 0);
% Segmentation quality
plot_info_cell(:, end+1) = {'link_all_stat.mean.has_no_ep_Q', 'Fraction_vsl_wo_ep',...
    'Fraction of links without endpoint', 0 : 0.01 : 1, 'linear'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.has_no_ep_Q', 'Fraction_cap_wo_ep',...
    'Fraction of capillary without endpoint', 0 : 0.01 : 1, 'linear'};

plot_info_cell(:, end+1) = {'link_all_stat.mean.has_ep_Q', 'Fraction_vsl_w_ep', ...
    'Fraction of links with endpoint', 10 .^ linspace(-2.8, 0, 50), 'log'};
plot_info_cell(:, end+1) = {'link_cap_stat.mean.has_ep_Q', 'Fraction_cap_w_ep', ...
    'Fraction of capillaries with endpoint', 10 .^ linspace(-2.8, 0, 50), 'log'};
%% Generate figures
for iter_feature = 1 : size(plot_info_cell, 2)
    %%
    tmp_plot_info = plot_info_cell(:, iter_feature);
    
    vis_field_name = tmp_plot_info{1};
    output_file_name = tmp_plot_info{2};
    vis_x_label_name = tmp_plot_info{3};
    stack_vis_x_val = tmp_plot_info{4};
    vis_x_scale = tmp_plot_info{5};
    
    stack_pdf_fit = struct;
    stack_pdf_fit.stack_list = stack_list;
    [stack_pdf_fit.hist_edge, stack_pdf_fit.interpolate_fun] = deal(cell(num_stack, 1));
    [stack_pdf_fit.mean, stack_pdf_fit.std, stack_pdf_fit.fraction_of_0] = deal(nan(num_stack, 1));
    
    stack_pdf_fit.data_cell = cell(num_stack, 1);
    stack_pdf_fit.feature_string_cell = cell(num_stack + 1, 1);
    
    for iter_stack = 1 : num_stack
        tmp_wb_data = wb_data_cube_stat_cell{iter_stack};
        tmp_is_internal_cube_Q = (tmp_wb_data.cube_in_brain_mask_ratio == stat_str.internal_cube_min_in_mask_vol_f) & ...
            (tmp_wb_data.cap2vsl_vol_fraction >= stat_str.minimal_cap2vsl_vol_ratio);
        tmp_hist_data = fun_getfield(tmp_wb_data, vis_field_name);
        if size(tmp_hist_data, 2) > 1
            tmp_hist_data = tmp_hist_data(:, 1);
        end
        tmp_hist_data = tmp_hist_data(tmp_is_internal_cube_Q & ~isnan(tmp_hist_data));
        if all(tmp_hist_data <= 1)
            tmp_vis_plt_max = 1;
            if all(tmp_hist_data >= 0)
                tmp_vis_plt_min = 0;
            end
        end
        
        tmp_stat_str = fun_analysis_get_basic_statistics(tmp_hist_data);
        
        stack_pdf_fit.data_cell{iter_stack} = tmp_stat_str;
        stack_pdf_fit.fraction_of_0(iter_stack) = nnz(tmp_hist_data == 0) / numel(tmp_hist_data);
        stack_pdf_fit.mean(iter_stack) = tmp_stat_str.mean;
        stack_pdf_fit.std(iter_stack) = tmp_stat_str.std;
        
        tmp_stack_string = sprintf('%s Mean: %.1f%%\nNo endpoint: %.1f%%', strrep(stack_list{iter_stack}, ...
            '_', ''), tmp_stat_str.mean * 100, stack_pdf_fit.fraction_of_0(iter_stack) * 100);
        
        stack_pdf_fit.feature_string_cell{iter_stack} = tmp_stack_string;
        stack_pdf_fit.hist_edge{iter_stack} = tmp_stat_str.hist_edge;
        stack_pdf_fit.interpolate_fun{iter_stack} = griddedInterpolant(tmp_stat_str.hist_edge, ...
            [0, tmp_stat_str.hist_cdf], 'linear', 'nearest');
    end
    stack_pdf_fit.feature_name = output_file_name;
    stack_pdf_fit.filepath = fullfile(stat_str.info.folder, 'feature_stat_hist', ...
        sprintf('%s_%s_%s_%s_CDF_fit_data.mat', dataset_name, merge_stack_name, stat_data_name, ...
        strrep(stack_pdf_fit.feature_name, ' ' , '_')));
    
    stack_vis_x_cell = cellfun(@(x) x.hist_bin_val, stack_pdf_fit.data_cell, 'UniformOutput', false);
    stack_vis_y_cell = cellfun(@(x) x.hist_cdf, stack_pdf_fit.data_cell, 'UniformOutput', false);
    
    stack_vis_itp_str = fun_analysis_get_xy_curve_avgNstd_by_interpolation(...
        stack_vis_x_cell, stack_vis_y_cell, stack_vis_x_val);
    
    fig_hdl = figure('Visible', 'on');
    % fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
    ax_hdl = axes(fig_hdl);
    for iter_stack = 1 : num_stack
        tmp_y_itp = stack_vis_itp_str.y_interpolation{iter_stack};
        tmp_y = tmp_y_itp(stack_vis_x_val);
        scatter(ax_hdl, stack_vis_x_val, tmp_y, 'o', 'filled');
        hold(ax_hdl, 'on');
    end
    [ax_hdl, line_hdl, error_area_hdl] = fun_vis_errorbar_shaded(stack_vis_x_val, ...
        stack_vis_itp_str.y_avg_interpolation(stack_vis_x_val),...
        stack_vis_itp_str.y_std_interpolation(stack_vis_x_val), ax_hdl);
    stack_pdf_fit.feature_string_cell{end} = sprintf('Average: %.1f%%\nNo endpoint: %.1f%%', ...
        mean(stack_pdf_fit.mean) * 100, mean(stack_pdf_fit.fraction_of_0) * 100);
    legend(ax_hdl, stack_pdf_fit.feature_string_cell{:}, 'Location', 'best');
    
    ax_hdl.XScale = vis_x_scale;
    ax_hdl.YLim = [0, 1];
    ax_hdl.XLabel.String = vis_x_label_name;
    ax_hdl.YLabel.String = plot_type;
    ax_hdl.FontSize = 14;
    grid(ax_hdl, 'off');
    box(ax_hdl, 'off');
    
    fun_print_image_in_several_formats(fig_hdl, strrep(stack_pdf_fit.filepath, '_fit_data.mat', '.png'));
    delete(fig_hdl);
end
%% Capillary - vessel volume fraction
plot_type = 'CDF';
plot_info_cell = cell(5, 0);
plot_info_cell(:, end+1) = {'cap2vsl_length_fraction', 'Cap2vsl_length_ratio', ...
    'Capillary-vessel length ratio', 0 : 0.05 : 1, 'linear'};
plot_info_cell(:, end+1) = {'cap2vsl_vol_fraction', 'Cap2vsl_volume_ratio', ...
    'Capillary-vessel volume ratio', 0 : 0.05 : 1, 'linear'};
plot_info_cell(:, end+1) = {'noncap2vsl_vol_fraction', 'Noncap2vsl_volume_ratio', ...
    'Noncapillary-vessel volume ratio', 0 : 0.05 : 1, 'linear'};
plot_info_cell(:, end+1) = {'cap2vsl_surf_area_fraction', 'Cap2vsl_surfArea_ratio', ...
    'Capillary-vessel surface area ratio', 0 : 0.05 : 1, 'linear'};
%% Generate figures
for iter_feature = 1 : size(plot_info_cell, 2)
    tmp_plot_info = plot_info_cell(:, iter_feature);
    tmp_plot_feature_name = tmp_plot_info{1};
    tmp_output_fn = tmp_plot_info{2};
    tmp_xlabel = tmp_plot_info{3};
    stack_vis_x_val = tmp_plot_info{4};
    tmp_x_scale = tmp_plot_info{5};
    
    stack_pdf_fit = struct;
    stack_pdf_fit.stack_list = stack_list;
    [stack_pdf_fit.hist_edge, stack_pdf_fit.interpolate_fun] = deal(cell(num_stack, 1));
    [stack_pdf_fit.mean, stack_pdf_fit.std, stack_pdf_fit.fraction_of_0] = deal(nan(num_stack, 1));
    
    stack_pdf_fit.data_cell = cell(num_stack, 1);
    stack_pdf_fit.feature_string_cell = cell(num_stack + 1, 1);
    
    for iter_stack = 1 : num_stack
        tmp_wb_data = wb_data_cube_stat_cell{iter_stack};
        tmp_is_internal_cube_Q = (tmp_wb_data.cube_in_brain_mask_ratio == stat_str.internal_cube_min_in_mask_vol_f) & ...
            (tmp_wb_data.cap2vsl_vol_fraction >= stat_str.minimal_cap2vsl_vol_ratio);
        tmp_hist_data = fun_getfield(tmp_wb_data, tmp_plot_feature_name);
        if size(tmp_hist_data, 2) > 1
            tmp_hist_data = tmp_hist_data(:, 1);
        end
        tmp_hist_data = tmp_hist_data(tmp_is_internal_cube_Q & ~isnan(tmp_hist_data));
        if all(tmp_hist_data <= 1)
            tmp_vis_plt_max = 1;
            if all(tmp_hist_data >= 0)
                tmp_vis_plt_min = 0;
            end
        end
        
        tmp_stat_str = fun_analysis_get_basic_statistics(tmp_hist_data);
        
        stack_pdf_fit.data_cell{iter_stack} = tmp_stat_str;
        stack_pdf_fit.fraction_of_0(iter_stack) = nnz(tmp_hist_data == 0) / numel(tmp_hist_data);
        stack_pdf_fit.mean(iter_stack) = tmp_stat_str.mean;
        stack_pdf_fit.std(iter_stack) = tmp_stat_str.std;
        
        tmp_stack_string = sprintf('%s Mean: %.1f%%', strrep(stack_list{iter_stack}, ...
            '_', ''), tmp_stat_str.mean * 100);
        
        stack_pdf_fit.feature_string_cell{iter_stack} = tmp_stack_string;
        stack_pdf_fit.hist_edge{iter_stack} = tmp_stat_str.hist_edge;
    end
    stack_vis_x_cell = cellfun(@(x) x.hist_edge, stack_pdf_fit.data_cell, 'UniformOutput', false);
    stack_vis_y_cell = cellfun(@(x) [0, x.hist_cdf], stack_pdf_fit.data_cell, 'UniformOutput', false);
    stack_vis_itp_str = fun_analysis_get_xy_curve_avgNstd_by_interpolation(...
        stack_vis_x_cell, stack_vis_y_cell, stack_vis_x_val);
    
    fig_hdl = figure('Visible', 'on');
    % fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.5;
    ax_hdl = axes(fig_hdl);
    for iter_stack = 1 : num_stack
        tmp_y_itp = stack_vis_itp_str.y_interpolation{iter_stack};
        tmp_y = tmp_y_itp(stack_vis_x_val);
        scatter(ax_hdl, stack_vis_x_val, tmp_y, 'o', 'filled');
        hold(ax_hdl, 'on');
    end
    [ax_hdl, line_hdl, error_area_hdl] = fun_vis_errorbar_shaded(stack_vis_x_val, ...
        stack_vis_itp_str.y_avg_interpolation(stack_vis_x_val),...
        stack_vis_itp_str.y_std_interpolation(stack_vis_x_val), ax_hdl);
    
    ax_hdl.XLabel.String = tmp_xlabel;
    ax_hdl.YLabel.String = plot_type;
    ax_hdl.FontSize = 14;
    grid(ax_hdl, 'off');
    box(ax_hdl, 'off');
    
    stack_pdf_fit.feature_name = tmp_output_fn;
    stack_pdf_fit.filepath = fullfile(stat_str.info.folder, 'feature_stat_hist', ...
        sprintf('%s_%s_%s_%s_CDF_fit_data.mat', dataset_name, merge_stack_name, stat_data_name, ...
        strrep(stack_pdf_fit.feature_name, ' ' , '_')));
    
    stack_pdf_fit.feature_string_cell{end} = sprintf('Average: %.1f%%\n', ...
        mean(stack_pdf_fit.mean) * 100);
    legend(ax_hdl, stack_pdf_fit.feature_string_cell{:}, 'Location', 'best');
    
    fun_print_image_in_several_formats(fig_hdl, strrep(stack_pdf_fit.filepath, '_fit_data.mat', '.png'));
    delete(fig_hdl);
end
%% Merge all PDF from all stack for visualization
%% Boxplot - whole brain statistics - cell data
stat_str.cell_data = [];
stat_str.feature_scalar_plot_setting = [];
fprintf('Generating plot setting cella array for nonscalar data\n');
box_data_extraction_info = cell(7, 0);
box_data_extraction_info(:, end+1) = {'node_stat', 'path_to_nearest_neighbor_geodesic_length', 1 : 6, 'Number of edges to the nearest node', 'linear', 'log', []};
box_data_extraction_info(:, end+1) = {'link_all_stat', 'shortest_loop_geodesic_length', 2 : 16, 'Number of edges in the shortest loop', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'shortest_loop_geodesic_length', 2 : 16, 'Number of edges in the shortest loop', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'capillary_branching_order', 1 : 9, 'Capillary branch order', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'degree', 3 : 5, 'Node degree', 'linear', 'log', []};

box_data_extraction_info(:, end+1) = {'node_stat', 'link_ori_ep_sin_elevation', 0.05 : 0.1 : 0.95, '$|\hat{n}_{t12}\cdot\vec{v}_{t3}|$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'link_ep2ep_sin_elevation', 0.05 : 0.1 : 0.95, '$|\hat{n}_{e12}\cdot\vec{v}_{e3}|$', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'node_stat', 'link_ori_ep_cos_elevation', 0.05 : 0.1 : 0.95, '$|sin(\theta_{t})|$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'link_ep2ep_cos_elevation', 0.05 : 0.1 : 0.95, '$|sin(\theta_{ep}|$', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'node_stat', 'link_ep2ep_agl_min', 10 : 20 : 170, '$\phi_{min}\;(^{\circ})$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'link_ep2ep_agl_med', 10 : 20 : 170, '$\phi_{median}\;(^{\circ})$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'link_ep2ep_agl_max', 10 : 20 : 170, '$\phi_{max}\;(^{\circ})$', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'node_stat', 'link_ori_agl_min', 10 : 20 : 170, '$\phi_{t,min}\;(^{\circ})$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'link_ori_agl_med', 10 : 20 : 170, '$\phi_{t,median}\;(^{\circ})$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'link_ori_agl_max', 10 : 20 : 170, '$\phi_{t,max}\;(^{\circ})$', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'node_stat', 'conn_node_drk_min', 1 : 1 : 10, 'Minimum node neighbor distance rank', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'conn_node_drk_med', 2.5 : 5 : 52.5, 'Median node neighbor distance rank', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'node_stat', 'conn_node_drk_max', 5 : 10 : 95, 'Maximum node neighbor distance rank', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'node_stat', 'nearest_node_dist', 5 : 10 : 55, 'Distance to the nearest node (\mum)', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'straightness', 0.05 : 0.1 : 0.95, 'Capilalry segment straightness', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_all_stat', 'straightness', 0.05 : 0.1 : 0.95, 'Vessel segment straightness', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'tortuosity', 1.05 : 0.1 : 1.95, 'Capilalry segment tortuosity', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_all_stat', 'tortuosity', 1.05 : 0.1 : 1.95, 'Vessel segment tortuosity', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'link_cap_stat', 'nearest_tissue_dt_max', 10 : 5 : 50, '$d_{max}\;(\mu m)$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'nearest_tissue_radius', 5 : 5 : 30, '$d_{r}\;(\mu m)$', 'linear', 'linear', []};
box_data_extraction_info(:, end+1) = {'link_cap_stat', 'nearest_tissue_dt_mean', 5 : 5 : 25, '$d_{avg}\;(\mu m)$', 'linear', 'linear', []};

box_data_extraction_info(:, end+1) = {'link_cap_stat', 'dt_median', 1 : 0.5 : 3.5, 'Capillary radius (\mum)', 'linear', 'linear', []};
num_box_data = size(box_data_extraction_info, 2);

for iter_box = 1 : num_box_data
    tmp_field_name = box_data_extraction_info{1, iter_box};
    tmp_feature_name = box_data_extraction_info{2, iter_box};
    tmp_bin_val = box_data_extraction_info{3, iter_box};
    [tmp_pdf_cell, tmp_y_bin_edge_cell, tmp_valid_Q_cell, tmp_mean_cell] = deal(cell(num_stack, 1));
    for iter_stack = 1 : num_stack
        tmp_pdf_cell{iter_stack} = fun_getfield(wb_data_cube_pdf_str{iter_stack}, ...
            sprintf('%s.hist_pdf.%s', tmp_field_name, tmp_feature_name));
        tmp_y_bin_edge_cell{iter_stack} = fun_getfield(wb_data_cube_pdf_str{iter_stack}, ...
            sprintf('%s.hist_edge.%s', tmp_field_name, tmp_feature_name));
        tmp_valid_Q_cell{iter_stack} = (wb_data_cube_stat_cell{iter_stack}.cube_in_brain_mask_ratio == ...
            stat_str.internal_cube_min_in_mask_vol_f) & ...
            (wb_data_cube_stat_cell{iter_stack}.cap2vsl_vol_fraction >= stat_str.minimal_cap2vsl_vol_ratio);
        tmp_mean_cell{iter_stack} = fun_getfield(wb_data_cube_stat_cell{iter_stack}, ...
            sprintf('%s.mean.%s', tmp_field_name, tmp_feature_name));
    end
    tmp_mean_data = cat(1, tmp_mean_cell{:});
    tmp_pdf_data = cat(1, tmp_pdf_cell{:});
    tmp_y_bin_edge = cat(1, tmp_y_bin_edge_cell{:});
    tmp_valid_Q = cat(1, tmp_valid_Q_cell{:});
    tmp_pdf_data = tmp_pdf_data(tmp_valid_Q);
    tmp_y_bin_edge = tmp_y_bin_edge(tmp_valid_Q);
    tmp_mean_data = tmp_mean_data(tmp_valid_Q);
    
    tmp_bin_data_str = fun_analysis_bin_network_properties_pdf(tmp_pdf_data, tmp_y_bin_edge, tmp_bin_val);
    box_data_extraction_info{end, iter_box} = tmp_bin_data_str.data;
    tmp_stat_data = rmfield(tmp_bin_data_str, 'data');
    tmp_stat_data.data_name = sprintf('%s_%s', tmp_field_name, tmp_feature_name);
    tmp_stat_data.data_label = box_data_extraction_info{4, iter_box};
    tmp_stat_data.mean = mean(tmp_mean_data, 'omitnan');
    tmp_stat_data.std = std(tmp_mean_data, 'omitnan');
    stat_str.cell_data{end+1} = tmp_stat_data;
end
stat_str.feature_box_plot_setting = box_data_extraction_info;
fprintf('Finish constructing plot information structure\n');
%%
boxplot_Q = true;
box_stat = cell(size(stat_str.feature_box_plot_setting, 2), 1);
tic
for iter_plot = 1 : size(stat_str.feature_box_plot_setting, 2)
    % for iter_plot = 8 : 13
    tmp_data = stat_str.feature_box_plot_setting{7, iter_plot};
    tmp_xlabel = stat_str.feature_box_plot_setting{4, iter_plot};
    tmp_tick = stat_str.feature_box_plot_setting{3, iter_plot};
    tmp_feature_name = sprintf('%s_%s', stat_str.feature_box_plot_setting{1, iter_plot}, ...
        stat_str.feature_box_plot_setting{2, iter_plot});
    tmp_xscale = stat_str.feature_box_plot_setting{5, iter_plot};
    tmp_yscale = stat_str.feature_box_plot_setting{6, iter_plot};
    tmp_stat_data = stat_str.cell_data{iter_plot};
    
    fig_hdl = figure('Visible', 'on');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1;
    ax_hdl = axes(fig_hdl);
    if boxplot_Q
        boxplot(ax_hdl, tmp_data, 'Symbol', 'r.', 'Widths', 1, 'Whisker', 1.5);
    end
    %         [tmp_violine_str, ax_hdl] = violinplot(tmp_data, tmp_tick, 'ShowData', false, 'ViolinAlpha', 0.6, 'ShowNotches', true);
    %         patch_hdl = patch(ax_hdl, [tmp_stat_data.mean - tmp_stat_data.std, tmp_stat_data.mean + tmp_stat_data.std, ...
    %             tmp_stat_data.mean + tmp_stat_data.std, tmp_stat_data.mean - tmp_stat_data.std], ...
    %             repelem([0, 1], 1, 2), 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
    %          hold(ax_hdl, 'on');
    %         line_hdl = line(ax_hdl, [tmp_stat_data.mean, tmp_stat_data.mean], ...
    %             [0, 1], 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-.');
    
    ax_hdl.YLabel.String = 'PDF';
    ax_hdl.FontSize = 14;
    if any(strfind(tmp_xlabel, '$'))
        ax_hdl.XLabel.Interpreter = 'latex';
        ax_hdl.XLabel.FontSize = 18;
    else
        ax_hdl.XLabel.Interpreter = 'tex';
    end
    ax_hdl.XLabel.String = tmp_xlabel;
    ax_hdl.YScale = tmp_yscale;
    ad_hdl.XScale = tmp_xscale;
    
    tmp_tick_value_diff = unique(round(diff(tmp_tick), 6, 'significant'));
    assert(isscalar(tmp_tick_value_diff), 'Uneven spacing of bin value');
    if tmp_tick_value_diff == 1
        tmp_tick_disp_value = tmp_tick;
    else
        tmp_tick_pos_value = [ax_hdl.XTick - 0.5, ax_hdl.XTick(end) + 0.5];
        tmp_tick_disp_value = [tmp_tick - 0.5 * tmp_tick_value_diff, tmp_tick(end) + 0.5 * tmp_tick_value_diff];
        
        ax_hdl.XTick = tmp_tick_pos_value;   
    end
    if all(mod(tmp_tick_disp_value, 1) == 0)
        ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%d'), tmp_tick_disp_value, 'UniformOutput', false);
    else
        ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), tmp_tick_disp_value, 'UniformOutput', false);
    end
    ax_hdl.XTickLabelRotation = 90;
    box(ax_hdl, 'off');
    tmp_stat_str_cell = num2cell(tmp_data, 1);
    tmp_stat_str_cell = cellfun(@fun_analysis_get_basic_statistics, tmp_stat_str_cell, 'UniformOutput', false);
    
    tmp_not_outlier_max = nan(size(tmp_data, 2), 1);
    tmp_not_outlier_min = nan(size(tmp_data, 2), 1);
    for iter_bin = 1 : numel(tmp_not_outlier_max)
        [~, tmp_not_outlier_min(iter_bin), tmp_not_outlier_max(iter_bin)]  = fun_analysis_is_outlier_by_percentile(tmp_data(:, iter_bin), 1.5);
    end
    tmp_not_outlier_max = min(max(tmp_data(:)), max(tmp_not_outlier_max));
    %         tmp_not_outlier_min = max(0, min(tmp_not_outlier_min) * 0.95);
    tmp_min_idx = find(ax_hdl.YTick >= tmp_not_outlier_max, 1, 'first');
    if ~isempty(tmp_min_idx)
        ax_hdl.YLim(2) = ax_hdl.YTick(tmp_min_idx);
    end
    ax_hdl.YLim(1) = 0;
    % Add whole brain Average + std of local average
    %         leg_hdl = legend([line_hdl, patch_hdl], sprintf('Mean: %.2e', tmp_stat_data.mean), ...
    %             sprintf('STD: %.2e (%.1f%%)', tmp_stat_data.std, tmp_stat_data.std / tmp_stat_data.mean * 100));
    ax_hdl.Title.String = sprintf('Mean: %.2e \\pm %.1e (%.1f%%)', tmp_stat_data.mean, ...
        tmp_stat_data.std, tmp_stat_data.std / tmp_stat_data.mean * 100);
    set(ax_hdl.Children.Children, 'MarkerEdgeColor', 'b');
    set(ax_hdl.Children.Children, 'MarkerSize', 5);
    set(ax_hdl.Children.Children, 'LineWidth', 2);
    ax_hdl.LineWidth = 2;
    
    fig_fp = fullfile(stat_str.info.folder, 'feature_pdf_boxplot', sprintf('%s_%s_%s_%s_PDF_boxplot.png', ...
        dataset_name, merge_stack_name, stat_data_name, tmp_feature_name));
    %%
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    delete(fig_hdl);
end
fprintf('Finish writing box plot figures\n');
toc
save(stat_str.info.filepath, '-struct', 'stat_str');
%% For presentation
boxplot_Q = true;
box_stat = cell(size(stat_str.feature_box_plot_setting, 2), 1);
tic
for iter_plot = 1 : size(stat_str.feature_box_plot_setting, 2)
    tmp_data = stat_str.feature_box_plot_setting{7, iter_plot};
    tmp_xlabel = stat_str.feature_box_plot_setting{4, iter_plot};
    tmp_tick = stat_str.feature_box_plot_setting{3, iter_plot};
    tmp_feature_name = sprintf('%s_%s', stat_str.feature_box_plot_setting{1, iter_plot}, ...
        stat_str.feature_box_plot_setting{2, iter_plot});
    tmp_xscale = stat_str.feature_box_plot_setting{5, iter_plot};
    tmp_yscale = stat_str.feature_box_plot_setting{6, iter_plot};
    tmp_stat_data = stat_str.cell_data{iter_plot};
    if boxplot_Q
        fig_hdl = figure('Visible', 'on');
        fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 1.2;
        ax_hdl = axes(fig_hdl);
        %         patch_hdl = patch(ax_hdl, [tmp_stat_data.mean - tmp_stat_data.std, tmp_stat_data.mean + tmp_stat_data.std, ...
        %             tmp_stat_data.mean + tmp_stat_data.std, tmp_stat_data.mean - tmp_stat_data.std], ...
        %             repelem([0, 1], 1, 2), 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
        %          hold(ax_hdl, 'on');
        %         line_hdl = line(ax_hdl, [tmp_stat_data.mean, tmp_stat_data.mean], ...
        %             [0, 1], 'LineWidth', 2, 'Color', 'k', 'LineStyle', '-.');
        boxplot(ax_hdl, tmp_data, 'Symbol', 'r.', 'Widths', 1);
        %         [tmp_violine_str, ax_hdl] = violinplot(tmp_data, tmp_tick, 'ShowData', false, 'ViolinAlpha', 0.6, 'ShowNotches', true);
        ax_hdl.YLabel.String = 'PDF';
        ax_hdl.FontSize = 14;
        if any(strfind(tmp_xlabel, '$'))
            ax_hdl.XLabel.Interpreter = 'latex';
            ax_hdl.XLabel.FontSize = 18;
        else
            ax_hdl.XLabel.Interpreter = 'tex';
        end
        ax_hdl.XLabel.String = tmp_xlabel;
        ax_hdl.YScale = tmp_yscale;
        ad_hdl.XScale = tmp_xscale;
        if all(mod(tmp_tick, 1) == 0)
            ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%d'), tmp_tick, 'UniformOutput', false);
        else
            ax_hdl.XTickLabel = arrayfun(@(x) num2str(x, '%.2f'), tmp_tick, 'UniformOutput', false);
        end
        box(ax_hdl, 'off');
        tmp_stat_str_cell = num2cell(tmp_data, 1);
        tmp_stat_str_cell = cellfun(@fun_analysis_get_basic_statistics, tmp_stat_str_cell, 'UniformOutput', false);
        [~, tmp_min_idx] = min(abs(ax_hdl.YTick - max(prctile(tmp_data, 99.5, 1))));
        ax_hdl.YLim(2) = ax_hdl.YTick(tmp_min_idx);
        ax_hdl.YLim(1) = 0;
        % Add whole brain Average + std of local average
        %         leg_hdl = legend([line_hdl, patch_hdl], sprintf('Mean: %.2e', tmp_stat_data.mean), ...
        %             sprintf('STD: %.2e (%.1f%%)', tmp_stat_data.std, tmp_stat_data.std / tmp_stat_data.mean * 100));
        %         ax_hdl.Title.String = sprintf('Mean: %.2e \\pm %.1e (%.1f%%)', tmp_stat_data.mean, ...
        %             tmp_stat_data.std, tmp_stat_data.std / tmp_stat_data.mean * 100);
        set(ax_hdl.Children.Children, 'MarkerEdgeColor', 'b');
        set(ax_hdl.Children.Children, 'MarkerSize', 5);
        set(ax_hdl.Children.Children, 'LineWidth', 2);
        ax_hdl.LineWidth = 2;
        
        fig_fp = fullfile(stat_str.info.folder, 'feature_pdf_boxplot_wo_title', sprintf('%s_%s_%s_%s_PDF_boxplot.png', ...
            dataset_name, merge_stack_name, stat_data_name, tmp_feature_name));
        fun_print_image_in_several_formats(fig_hdl, fig_fp);
        delete(fig_hdl);
    end
end
fprintf('Finish writing box plot figures\n');
toc
save(stat_str.info.filepath, '-struct', 'stat_str');
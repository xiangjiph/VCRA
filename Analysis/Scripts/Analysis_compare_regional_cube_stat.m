clc;clear;close all;
%%
DataManager = FileManager;
dataset_name = 'WholeBrain';
image_grid_version = '240_cube';
stack_list = {'ML_2018_08_15', 'ML20190124', 'ML20200201'};
skel_version = '240_cube_rec';
pO2_folder_name = 'pO2_min_r_2um_sc';
reconstruction_version = '240_cube_recon_sc';
wb_stat_folder_name = 'whole_brain_stat_sc';

plot_stack_list = cellfun(@(x) strrep(x, '_', ''), stack_list, 'UniformOutput', false);
linear_scaling_factor = 1.0521;
num_stack = numel(stack_list);
allen_atlas = load('Allen_atlas.mat');
registration_version = 'Allen_2017_25um_landmark.mat';
Allen_atlas_id = load('Allen_atlas_id.mat');
tic_load_data = tic;
wb_data_cell = cell(num_stack, 1);
for iter_stack = 1 : num_stack
    tmp_tic = tic;
    tmp_stack = stack_list{iter_stack};
    wb_region_stat_data = fun_analysis_load_whole_brain_data_for_regional_analysis(...
        dataset_name, tmp_stack, image_grid_version, reconstruction_version, ...
        skel_version, registration_version, wb_stat_folder_name, false);
    wb_region_stat_data.cube_stat.num_link = wb_region_stat_data.cube_stat.link_all_stat.num_data.length;
    wb_region_stat_data.cube_stat.num_cap = wb_region_stat_data.cube_stat.link_cap_stat.num_data.length;
    
    tmp_pO2_stat_filepath = fullfile(DataManager.fp_analysis_data_folder(dataset_name, tmp_stack), ...
        sprintf('%s_%s_%s_%s_stat_data.mat', dataset_name, tmp_stack, skel_version, pO2_folder_name));
    wb_region_stat_data.cube_stat.pO2_data = DataManager.load_data(tmp_pO2_stat_filepath);
    %%
    wb_data_cell{iter_stack} = wb_region_stat_data;
    clearvars wb_region_stat_data
    fprintf('Finish loading grid, registration and brain mask distance transform for stack %s. Elapsed time is %f seconds\n', tmp_stack, toc(tmp_tic));
end
fprintf('Finish loading data. Elapsed time is %f seconds.\n', toc(tic_load_data));
%% Feature extraction options
% Compare capillary length density and length distribution in major brain
% regions 
de_opt_str = struct;
de_opt_str.cube_stat_Q = true;
de_opt_str.node_feature_Q = false;
de_opt_str.link_feature_Q = false;
de_opt_str.vessel_graph_Q = false;
de_opt_str.depth_dependent_density_Q = false;
de_opt_str.merge_cc_Q = false;
de_opt_str.save_Q = false;
de_opt_str.save_fp = [];
merge_stack_name = 'all_stack';
allen_atlas_id = load('Allen_atlas_id.mat');
scatter_marker_type_list = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
num_maker_type = numel(scatter_marker_type_list);
if de_opt_str.merge_cc_Q
    num_cc = 1;
else
    num_cc = 2;
end
%%
allen_atlas_map_old_id_to_new_id = wb_data_cell{1}.registration.map_oldID_to_newID;
parent_str_name = 'Major_brain_regions';
analysis_region_id_list = {453, 500, 247, 669, 895, 961 1080 477 549 [302 294] 4 771 354 528 776, 1097};
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
analysis_region_name{10} = 'Superior colliculus';

visualization_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, ...
    merge_stack_name), vis_subfolder_name, parent_str_name);
analysis_region_name_abbrv = analysis_region_name;
%% Collect all the region statistics in all the stack
region_data = cell(num_region, num_stack);
collect_data_tic = tic;
for iter_stack = 1 : num_stack
    for iter_region = 1 : num_region
        tmp_tic = tic;
        tmp_region_id = analysis_region_id_list{iter_region};
        tmp_region_name = strrep(analysis_region_name{iter_region}, ' ', '_');
        tmp_wb_data = wb_data_cell{iter_stack};
        tmp_region_data = fun_analysis_get_atlas_regional_data(tmp_wb_data, ...
            tmp_region_id, tmp_region_name, de_opt_str);
        %% Process region pO2 data
        tmp_pO2_data_field_name = setdiff(fieldnames(tmp_wb_data.cube_stat.pO2_data), ...
            {'dataset_name', 'stack', 'filepath', 'stack', 'skel_version', 'pO2_data_folder'});
        if all(isfield(tmp_region_data, {'cc_1', 'cc_2'}))
            for iter_cc = 1 : tmp_region_data.num_cc
                tmp_cc_data = tmp_region_data.(sprintf('cc_%d', iter_cc));
                tmp_cc_data.local_cube_stat.num_cap = tmp_cc_data.local_cube_stat.link_cap_stat.num_data.length;
                tmp_cc_data.local_cube_stat.num_link = tmp_cc_data.local_cube_stat.link_all_stat.num_data.length;
                tmp_cc_data.local_cube_stat.node_density_mm3 = tmp_cc_data.local_cube_stat.node_stat.num_data.degree ./ (0.24 * linear_scaling_factor)^3;
                tmp_cc_data.local_cube_stat = rmfield(tmp_cc_data.local_cube_stat, {'mask_volume', 'num_surface_voxel'});
                
                for iter_field_name = 1 : numel(tmp_pO2_data_field_name)
                    tmp_field_name = tmp_pO2_data_field_name{iter_field_name};
                    tmp_cc_data.pO2_data.(tmp_field_name) = fun_structure_field_indexing(tmp_wb_data.cube_stat.pO2_data.(tmp_field_name), ...
                        tmp_cc_data.is_in_cc_cube_Q);
                end
                %% Compute local statistics
                tmp_cc_data.local_cube_stat.Krogh_corr_coeff = tmp_cc_data.pO2_data.fit_krogh.coefficient;
                
                tmp_cc_data.local_cube_stat.local_dt_mean = tmp_cc_data.pO2_data.local_dt_stat.mean;
                tmp_cc_data.local_cube_stat.local_dt_median = tmp_cc_data.pO2_data.local_dt_stat.median;
                
                tmp_cc_data.local_cube_stat.local_pO2_mean = tmp_cc_data.pO2_data.local_pO2_stat.mean;
                tmp_cc_data.local_cube_stat.local_pO2_median = tmp_cc_data.pO2_data.local_pO2_stat.median;
                
                tmp_cc_data.local_cube_stat.dt_lm_dt_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.dt_lm.v);
                tmp_cc_data.local_cube_stat.dt_lm_dt_v_median = cellfun(@(x) double(median(x)), tmp_cc_data.pO2_data.dt_lm.v);
                tmp_cc_data.local_cube_stat.dt_lm_pO2_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.dt_lm.pO2_v);
                tmp_cc_data.local_cube_stat.dt_lm_pO2_v_median = cellfun(@(x) double(median(x)), tmp_cc_data.pO2_data.dt_lm.pO2_v);
                
                tmp_cc_data.local_cube_stat.pO2_lm_dt_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.pO2_lm.dt_v);
                tmp_cc_data.local_cube_stat.pO2_lm_dt_v_median = cellfun(@(x) double(median(x)), tmp_cc_data.pO2_data.pO2_lm.dt_v);
                tmp_cc_data.local_cube_stat.pO2_lm_pO2_v_mean = cellfun(@(x) double(mean(x)), tmp_cc_data.pO2_data.pO2_lm.v);
                tmp_cc_data.local_cube_stat.pO2_lm_pO2_v_median = cellfun(@(x) double(median(x)), tmp_cc_data.pO2_data.pO2_lm.v);
                
                tmp_region_data.(sprintf('cc_%d', iter_cc)) = tmp_cc_data;
            end
        end
        region_data{iter_region, iter_stack} = tmp_region_data;
        toc(tmp_tic);
    end
end
fprintf('Finish collecting all the region data. Elapsed time is %f seconds.\n', toc(collect_data_tic));
%% Collect data
% 1. Capillary length density in 240-cube
% 2. Capillary segment length distribution PDF
% 3. Capillary-tissue maximum distance 
min_in_cc_ratio = 0.75;
min_in_mask_ratio = 1;
min_cap2vsl_vol_ratio = 0;
selection_folder_name = sprintf('min_in_cc_f_%.2f_min_cap2vsl_vol_r_%.2f', ...
    min_in_cc_ratio, min_cap2vsl_vol_ratio);
plot_data = struct;
[ plot_data.vessel_length_density, plot_data.capillary_length_density, ...
    plot_data.avg_dt_lm, plot_data.avg_d, ...
    plot_data.avg_u, plot_data.avg_u_lm, ...
    plot_data.fractional_anisotropy, plot_data.FA_z] = deal(cell(num_cc, num_stack, num_region));
for iter_region = 1 : num_region
    for iter_stack = 1 : num_stack
        for iter_cc = 1 : num_cc
            tmp_cc_data = region_data{iter_region, iter_stack}.(sprintf('cc_%d', iter_cc));
            tmp_240_cube_selected_Q = (tmp_cc_data.local_cube_stat.cube_in_brain_mask_ratio >= min_in_mask_ratio & ...
                tmp_cc_data.local_cube_stat.in_cc_vol_f >= min_in_cc_ratio) & ...
                tmp_cc_data.local_cube_stat.cap2vsl_vol_fraction >= min_cap2vsl_vol_ratio;
            plot_data.vessel_length_density{iter_cc, iter_stack, iter_region} = ...
                tmp_cc_data.local_cube_stat.link_length_density_m_mm3(tmp_240_cube_selected_Q);
            plot_data.capillary_length_density{iter_cc, iter_stack, iter_region} = ...
                tmp_cc_data.local_cube_stat.capillary_length_density_m_mm3(tmp_240_cube_selected_Q);
            
            plot_data.fractional_anisotropy{iter_cc, iter_stack, iter_region} = ...
                tmp_cc_data.local_cube_stat.wb_ai_cap_vw.fractional_anisotropy(tmp_240_cube_selected_Q);
            plot_data.FA_z{iter_cc, iter_stack, iter_region} = ...
                tmp_cc_data.local_cube_stat.wb_ai_cap_vw.fa_z(tmp_240_cube_selected_Q);
            
            plot_data.avg_d{iter_cc, iter_stack, iter_region} = ...
                tmp_cc_data.pO2_data.local_dt_stat.mean(tmp_240_cube_selected_Q);
            
            plot_data.avg_dt_lm{iter_cc, iter_stack, iter_region} = ...
                cellfun(@(x) double(mean(x, 'omitnan')), tmp_cc_data.pO2_data.dt_lm.v(tmp_240_cube_selected_Q));
            
            plot_data.avg_u{iter_cc, iter_stack, iter_region} = ...                
                tmp_cc_data.pO2_data.local_pO2_stat.mean(tmp_240_cube_selected_Q);
            plot_data.avg_u_lm{iter_cc, iter_stack, iter_region} = ...
                cellfun(@(x) double(mean(x, 'omitnan')), tmp_cc_data.pO2_data.pO2_lm.v(tmp_240_cube_selected_Q));
        end
    end
end
fprintf('Finish collecting data for visualization\n');
%% Compute the histogram 
capillary_hist = struct;
compute_hist_field_name = fieldnames(plot_data);
avg_info_cell = cell(0, 7);
% density_vis_region_list_ind = [1, 2, 7, 8, 11, 14];
density_vis_region_list_ind = [1, 15, 7, 8, 11, 14];
avg_info_cell(end+1, :) = {'vessel_length_density', 0 : 0.1 : 2.4, 'vessel_length_density', 'Vessel length density (m/mm^3)', 'linear', 'linear', density_vis_region_list_ind};
avg_info_cell(end+1, :) = {'capillary_length_density', 0 : 0.1 : 2.4, 'capillary_length_density', 'Capillary length density (m/mm^3)', 'linear', 'linear', density_vis_region_list_ind};
avg_info_cell(end+1, :) = {'avg_dt_lm', 0 : 2.5 : 60, 'avg_dt_lm', 'Average DT local maximum (\mum)', 'linear', 'linear', density_vis_region_list_ind};
avg_info_cell(end+1, :) = {'avg_d', 0 : 2.5 : 30, 'avg_dt', 'Average vessel-tissue distance (\mum)', 'linear', 'linear', density_vis_region_list_ind};

avg_info_cell(end+1, :) = {'avg_u_lm', -1500 : 50 : 0, 'avg_u_lm', 'Average normalized pO_2 minimum', 'linear', 'linear', density_vis_region_list_ind};
%%
for iter_field = 1 : size(avg_info_cell, 1)
    tmp_field_name = avg_info_cell{iter_field, 1};
    tmp_bin_edge = avg_info_cell{iter_field, 2};
    tmp_save_file_name = avg_info_cell{iter_field, 3};
    tmp_xlabel = avg_info_cell{iter_field, 4};
    tmp_xscale = avg_info_cell{iter_field, 5};
    tmp_yscale = avg_info_cell{iter_field, 6};
    tmp_overlay_structure_list_idx = avg_info_cell{iter_field, 7};
    
    tmp_hist_cell =  cellfun(@(x) fun_analysis_get_basic_statistics_with_edges(x, tmp_bin_edge, true), ...
        plot_data.(tmp_field_name), 'UniformOutput', false);
    capillary_hist.(tmp_field_name) = tmp_hist_cell;
    
    tmp_hist_itp_cell = cell(num_region, 1);
    tmp_mean = cellfun(@(x) x.mean, tmp_hist_cell);
    tmp_median = cellfun(@(x) x.median, tmp_hist_cell);
    tmp_std = cellfun(@(x) x.std, tmp_hist_cell);
    for iter_region = 1 : num_region
        tmp_feature_cell = tmp_hist_cell(:, :, iter_region);
        tmp_pdf_vec_cell = cellfun(@(x) x.hist_pdf, tmp_feature_cell, 'UniformOutput', false);
        tmp_bin_val_vec_cell = cellfun(@(x) x.hist_bin_val, tmp_feature_cell, 'UniformOutput', false);
        tmp_bin_edge_vec_cell = cellfun(@(x) x.hist_edge, tmp_feature_cell, 'UniformOutput', false);
        tmp_hist_itp_cell{iter_region} = fun_analysis_get_xy_curve_avgNstd_by_interpolation(tmp_bin_val_vec_cell, ...
            tmp_pdf_vec_cell, tmp_bin_val_vec_cell{1});
    end
    tmp_region_mean = mean(reshape(tmp_mean, 2 * num_stack, num_region), 1);
    [~, tmp_sorted_ind] = sort(tmp_region_mean, 'ascend');
%% Average histogram for each region - one structure one subplot
    num_row = ceil(sqrt(num_region));
    num_col = ceil(num_region / num_row);
    fig_hdl = figure('Visible', 'off');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [num_col, num_row];
    for iter_region = 1 : num_region
        ax_hdl = subplot(num_col, num_row, iter_region);
        tmp_region_ind = tmp_sorted_ind(iter_region);
        tmp_hist_str_cell = tmp_hist_cell(:, :, iter_region);
        tmp_hist_str_cell = tmp_hist_str_cell(:);
        hist_hdl = histogram(ax_hdl, 'BinCounts',  tmp_hist_itp_cell{tmp_region_ind}.y_avg, 'BinEdges', tmp_bin_edge);
        hold(ax_hdl, 'on');
        for iter_cell = 1 : numel(tmp_hist_str_cell)
            tmp_y = tmp_hist_itp_cell{tmp_region_ind}.y_interpolation{iter_cell};
            tmp_x = tmp_hist_itp_cell{tmp_region_ind}.interpolate_x;  
            scatter(ax_hdl, tmp_x, tmp_y.Values, 'ko');
        end        
        ebar_hdl = errorbar(ax_hdl, tmp_hist_itp_cell{tmp_region_ind}.interpolate_x, ...
           tmp_hist_itp_cell{tmp_region_ind}.y_avg, tmp_hist_itp_cell{tmp_region_ind}.y_std, 'LineStyle', 'none', 'Color', 'b', 'LineWidth', 2);
       
       tmp_legend_mean = tmp_mean(:, :, tmp_region_ind);
       tmp_legend_median = tmp_median(:, :, tmp_region_ind);
       tmp_legend_std = tmp_std(:, :, tmp_region_ind);
       
       tmp_legent_str = sprintf('%s\nMean: %.2f \\pm %.2f\nSD: %.2f \\pm %.2f', analysis_region_name{tmp_region_ind}, ...
           mean(tmp_legend_mean(:), 'all'), std(tmp_legend_mean(:)), ...
           mean(tmp_legend_std(:), 'all'), std(tmp_legend_std(:)));
       legend(ebar_hdl, tmp_legent_str);
       ax_hdl.XLabel.String = tmp_xlabel;
       ax_hdl.YLabel.String = 'PDF';      
       ax_hdl.YLim(1) = 0;
       ax_hdl.FontSize = 14;
       ax_hdl.XScale = tmp_xscale;
       ax_hdl.YScale = tmp_yscale;
    end
    tmp_fig_fp = fullfile(visualization_folder, selection_folder_name, sprintf('%s_%s_%s_%s_%s.png', ...
        dataset_name, merge_stack_name, parent_str_name, tmp_save_file_name, ...
        'pdf_ebar_per_str'));
    fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);    
    capillary_hist.(sprintf('avg_itp_%s', tmp_field_name)) = tmp_hist_itp_cell;
    delete(fig_hdl);
%% Average histogram for each region - overlay histogram 
    fig_hdl = figure('Visible', 'on');
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4);
    ax_hdl = axes(fig_hdl);
    leg_hdl_array = [];
%     tmp_plot_mean = tmp_region_mean(tmp_overlay_structure_list_idx);
%     [~, tmp_sorted_ind] = sort(tmp_plot_mean, 'ascend');
%     tmp_plot_list_idx = tmp_overlay_structure_list_idx(tmp_sorted_ind);
    tmp_plot_list_idx = tmp_overlay_structure_list_idx;
    for iter_region = 1 : numel(tmp_overlay_structure_list_idx)
        tmp_region_idx = tmp_plot_list_idx(iter_region);
        tmp_x = tmp_hist_itp_cell{tmp_region_idx}.interpolate_x;
        tmp_y = tmp_hist_itp_cell{tmp_region_idx}.y_avg;
        tmp_y_error = tmp_hist_itp_cell{tmp_region_idx}.y_std;
        [ax_hdl, plt_hdl, patch_hdl] = fun_vis_errorbar_shaded(tmp_x, tmp_y, tmp_y_error, ax_hdl);
        leg_hdl_array(end+1) = plt_hdl;
        hold(ax_hdl, 'on');
    end
    ax_hdl.XLabel.String = tmp_xlabel;
    ax_hdl.YLabel.String = 'PDF';
    ax_hdl.YLim(1) = 0;
    ax_hdl.XScale = tmp_xscale;
    ax_hdl.YScale = tmp_yscale;
    legend(ax_hdl, leg_hdl_array, analysis_region_name{tmp_plot_list_idx}, 'Location', 'northeast');
    tmp_fig_fp = fullfile(visualization_folder, selection_folder_name, sprintf('%s_%s_%s_%s_%s.png', ...
        dataset_name, merge_stack_name, parent_str_name, tmp_save_file_name, ...
        'selected_subregions_pdf_ebar'));
    fun_print_image_in_several_formats(fig_hdl, tmp_fig_fp);
    delete(fig_hdl);
end
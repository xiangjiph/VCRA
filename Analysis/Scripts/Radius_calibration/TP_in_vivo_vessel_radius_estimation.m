set_env
DataManager = FileManager;
%% Parameters
medfilt_tile_image_Q = true;
radius_estimation_method = 'DT';
%%
simulation_file_fp = fullfile(DataManager.fp_metadata_file('Vessel_radius_calibration', 'In_vivo_image_simulation', 'In_vivo_image_simulation_v2'));
invivo_im_prof_str = DataManager.load_data(simulation_file_fp);
%% Preprocess simulation result
num_parameter_combination = numel(invivo_im_prof_str.edge_int_str);
for iter_p = 1 : num_parameter_combination
    tmp_psf = invivo_im_prof_str.edge_int_str{iter_p};
    tmp_num_dist = numel(tmp_psf.normalized_radial_int_dist);
    tmp_psf.normalized_radial_int_dist_itp = cell(size(tmp_psf.normalized_radial_int_dist));
    tmp_psf.normalized_vertical_radial_int_dist_itp = cell(size(tmp_psf.normalized_vertical_radial_int_dist));
    
    for iter_dist = 1 : tmp_num_dist
        tmp_int_xy_dist = tmp_psf.max_normalized_radial_int_dist{iter_dist};
        if ~isempty(tmp_int_xy_dist)
            tmp_psf.normalized_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_xy_dist(:, 1), ...
                tmp_int_xy_dist(:, 2), 'linear', 'nearest');
            
            tmp_int_z_dist = tmp_psf.max_normalized_vertical_radial_int_dist{iter_dist};
            tmp_psf.normalized_vertical_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_z_dist(:, 1), ...
                tmp_int_z_dist(:, 2), 'linear', 'nearest');
        end
    end
    invivo_im_prof_str.edge_int_str{iter_p} = tmp_psf;
end
fprintf('Finish pre-computing the intensity interpolation\n');
%% Local two photon image peprocessing
% Convert raw iamge to frame-averaged image. Correct pixel shift due to scanner mechanical recovery error.
dataset_name = 'Vessel_radius_calibration';
% average_method = 'Rank_3_w_gmc';
average_method = 'Average';
image_group_list = {'In_vivo'};
stack_list = {'DK20200504_WT2', 'DK20200504_WT1', 'DK20200427_WT3', 'DK20200427_WT4', 'DK20200511_WT5', 'DK20200517_WT7'};
num_image_group = numel(image_group_list);
num_data_stack = numel(stack_list); 
%%
% Parallel computing:
num_processes = 16;
p_hdl = gcp('nocreate');
delete(p_hdl);
p_hdl = parpool(num_processes);
num_thread_per_prop = ceil(48  / num_processes);

for iter_data_stack = 1 : num_data_stack
    stack = stack_list{iter_data_stack};
    for iter_IG = 1 : num_image_group
        image_group = image_group_list{iter_IG};
        averaged_image_folder = fullfile(DataManager.fp_processed_data(dataset_name, ...
            stack), sprintf('%s_image', average_method), image_group);
        image_group_info = DataManager.load_data(fullfile(DataManager.fp_raw_data_folder(dataset_name, stack), ...
            image_group, 'image_group_info.mat'));
        averaged_image_file_str = dir(fullfile(averaged_image_folder, '*.tiff'));
        num_roi = numel(averaged_image_file_str);
        assert(num_roi == image_group_info.num_stack);
        %%
        parfor iter_roi = 1 : num_roi
            maxNumCompThreads(num_thread_per_prop);
            %% Load data
            tmp_tic = tic;
            tmp_filename = fullfile(averaged_image_folder,...
                sprintf('%s_%s_%s_%s_%d.tiff', dataset_name, stack, image_group, average_method, iter_roi));
            avg_data = DataManager.load_single_tiff(tmp_filename);
            
            tmp_filename = fullfile(averaged_image_folder,...
                sprintf('%s_%s_%s_%s_%d_info.mat', dataset_name, stack, image_group, average_method, iter_roi));
            avg_data_info = DataManager.load_data(tmp_filename);
            avg_im_size = size(avg_data);
            voxel_size_um = avg_data_info.voxel_size_um;
            tmp_grid_version = sprintf('240_cube_%s_%d', image_group, iter_roi);  
            
            tmp_grid = DataManager.load_grid(dataset_name, stack, tmp_grid_version);
            vessel_image = DataManager.load_blocks_files('image', dataset_name, stack, ...
                tmp_grid.version, 1 : tmp_grid.grid_size(1), ...
                1 : tmp_grid.grid_size(2), 1 : tmp_grid.grid_size(3));
            
            vessel_graph = DataManager.load_graph_in_block(dataset_name, stack, ...
                tmp_grid_version, 0, 0, 0);
            fprintf('Finish loading data. Elapsed time is %f seconds\n', toc(tmp_tic));
            %% Radius estimation with fixed PSF parameters    
            if medfilt_tile_image_Q
                avg_data = medfilt3(avg_data);
            end
            skeleton_cc_sub_um = cellfun(@(x) fun_ind2sub(vessel_graph.num.mask_size, x) .* ...
                vessel_graph.info.voxel_size_um, vessel_graph.link.cc_ind, 'UniformOutput', false);
            skeleton_cc_r_um = cellfun(@(x) double(full(vessel_graph.radius(x))), ...
                vessel_graph.link.cc_ind, 'UniformOutput', false);
            % Use vessel with RBC simulation
            try
                est_radius_str = fun_in_vivo_radius_estimation_in_image_stack(invivo_im_prof_str, ...
                    avg_data, voxel_size_um, skeleton_cc_sub_um, skeleton_cc_r_um, radius_estimation_method);
            catch ME
               fprintf('Unable to estimate radius in ROI %d\n', iter_roi);
               rethrow(ME);
%                fprintf('Error message: %s\n', getReport(ME,'extended','hyperlinks','off'));
%                fprintf('Skip saving the graph\n');
%                continue;
            end
            vessel_graph.info.radius_estimation.medfilt3_Q = medfilt_tile_image_Q;
            vessel_graph.info.radius_estimation.method = radius_estimation_method;
            
            vessel_graph.link.features.radius_median = est_radius_str.cc_stat.median;
            vessel_graph.link.features.radius_mean = est_radius_str.cc_stat.mean;
            vessel_graph.link.features.radius_std = est_radius_str.cc_stat.std;
            
            vessel_graph.link.features.cc_com_sub = cellfun(@(x) mean(x, 1), ...
                skeleton_cc_sub_um, 'UniformOutput', false);
            vessel_graph.link.cc_r = est_radius_str.r;
            vessel_graph.link.cc_int = cellfun(@(x) vessel_image(x), vessel_graph.link.cc_ind, 'UniformOutput', false);
            vessel_graph.link.cc_local_ori_vec_z = est_radius_str.local_ori_vec_z;
            vessel_graph.link.cc_est_parms = est_radius_str.best_fit_parameters;
            vessel_graph.link.cc_est_corr = est_radius_str.best_fit_corr;
            assert(all(cellfun(@numel, est_radius_str.local_ori_vec_z) == vessel_graph.link.num_voxel_per_cc)); 
            assert(all(cellfun(@(x) size(x, 1), est_radius_str.best_fit_parameters) == vessel_graph.link.num_voxel_per_cc)); 
            assert(all(cellfun(@numel, est_radius_str.best_fit_corr) == vessel_graph.link.num_voxel_per_cc)); 
            %% Visualization
%             recon_mask_0 = uint8(fun_graph_reconstruction(vessel_graph));
%             recon_mask_0(vessel_graph.link.pos_ind) = 2;
%             recon_mask_0(vessel_graph.node.pos_ind) = 3;            
%             DataManager.visualize_itksnap(vessel_image, recon_mask_0);
            %%
%             vis_link_sub = [95 102 429];
%             vis_link_ind = sub2ind(vessel_graph.num.mask_size, vis_link_sub(1), ...
%                 vis_link_sub(2), vis_link_sub(3));
%             vis_link_label = full(vessel_graph.link.map_ind_2_label(vis_link_ind));
%             fprintf('Link radius: %.3f, %.3f, %.3f, %.3f\n', est_radius_str.cc_stat.median(vis_link_label, :));
            %% Save vessel graph
            output_graph_version = sprintf('%s_re', tmp_grid_version);
            vessel_graph.info.graph_version = output_graph_version;
            DataManager.write_graph_in_block(vessel_graph, dataset_name, stack, ...
                output_graph_version, 0, 0, 0);
            fprintf('Finish radius estimation. Elapsed time is %f seconds\n', toc(tmp_tic));
        end
    end
    fprintf('Finish processing %s %s\n', dataset_name, stack);
end
%% Comparing Estimation accoridng to simulation with vs without RBC
% tmp_bin = 0 : 0.25 : 8;
% tmp_x = r_est_dt.cc_stat.median;
% tmp_y = est_radius_str.cc_stat.median;
% tmp_selected_Q = tmp_x < tmp_bin(end) & tmp_y < tmp_bin(end);
% tmp_x = tmp_x(tmp_selected_Q);
% tmp_y = tmp_y(tmp_selected_Q);
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% histogram2(ax_hdl, tmp_x, tmp_y, tmp_bin, tmp_bin, 'DisplayStyle', 'tile');
% ax_hdl.XLabel.String = 'Profile w RBC (\mum)';
% ax_hdl.YLabel.String = 'Profile w/o RBC (\mum)';
% hold(ax_hdl, 'on');
% plot(ax_hdl, tmp_bin, tmp_bin, 'k-.', 'LineWidth', 1);
% ax_hdl.FontSize = 14;
% ax_hdl.DataAspectRatio = [1,1,1];
% cbar_hdl = colorbar(ax_hdl);
% cbar_hdl.Label.String = 'Number of vessels';
% ax_hdl.ColorScale = 'log';
% 
% tmp_fit = fitlm(tmp_x, tmp_y, 'Intercept', false);
% leg_string = sprintf('Slope: %.2f\nR^2-Adj: %.2f', tmp_fit.Coefficients.Estimate(1), ...
%     tmp_fit.Rsquared.Adjusted);
% leg_hdl = legend(ax_hdl, leg_string, 'Location', 'northwest');
% fun_print_image_in_several_formats(fig_hdl, fullfile(DataManager.fp_visualization_folder(...
%     dataset_name, stack), 'Radius_estimation_comparison', sprintf('%s_%s_%s_%d_re_dt_wRBC_vs_woRBC_medfilt3.png', ...
%     dataset_name, stack, image_group, iter_roi)));
%% Debug
% check_cc_idx = find(abs(r_est_dt.cc_stat.median - est_radius_str.cc_stat.median) > 2);
%% 
% tmp_bin = 0 : 0.25 : 6;
% tmp_x = r_est_dt.cc_stat.median;
% tmp_y = r_est_fit.cc_stat.median;
% tmp_num_vxl_per_cc = cellfun(@numel, r_est_dt.r);
% tmp_selected_Q = tmp_x < tmp_bin(end) & tmp_y < tmp_bin(end) & tmp_num_vxl_per_cc > 8;
% tmp_x = tmp_x(tmp_selected_Q);
% tmp_y = tmp_y(tmp_selected_Q);
% 
% fig_hdl = figure;
% ax_hdl = axes(fig_hdl);
% histogram2(ax_hdl, tmp_x, tmp_y, tmp_bin, tmp_bin, 'DisplayStyle', 'tile');
% ax_hdl.XLabel.String = 'Radius estimated by DT (\mum)';
% ax_hdl.YLabel.String = 'Radius estiamted by profile (\mum)';
% hold(ax_hdl, 'on');
% plot(ax_hdl, tmp_bin, tmp_bin, 'k-.', 'LineWidth', 1);
% ax_hdl.FontSize = 14;
% ax_hdl.DataAspectRatio = [1,1,1];
% cbar_hdl = colorbar(ax_hdl);
% cbar_hdl.Label.String = 'Number of vessels';
% tmp_fit = fitlm(tmp_x, tmp_y, 'Intercept', false);
% ax_hdl.ColorScale = 'log';
% leg_string = sprintf('Slope: %.2f\nR^2-Adj: %.2f', tmp_fit.Coefficients.Estimate(1), ...
%     tmp_fit.Rsquared.Adjusted);
% leg_hdl = legend(ax_hdl, leg_string, 'Location', 'northwest');
% fun_print_image_in_several_formats(fig_hdl, fullfile(DataManager.fp_visualization_folder(...
%     dataset_name, stack), 'Radius_estimation_comparison', sprintf('%s_%s_%s_%d_re_dt_vs_profile_medfilt3.png', ...
%     dataset_name, stack, image_group, iter_roi))); 
% est_diff = r_est_dt.cc_stat.median - r_est_fit.cc_stat.median;
% check_cc_idx = find(abs(est_diff) > 2);
%% Visualize PSF estimation result
%             tmp_vis_data = psf_stat;
%             vis_link_ori_vec = cat(1, tmp_vis_data.link_ori_vec{:});
%             vis_best_fit_FWHM = cat(1, tmp_vis_data.best_fit_PSF_FWHM{:});
%             vis_r = cat(1, tmp_vis_data.refine_r_est{:});
%             vis_best_fit_corr = cat(1, tmp_vis_data.best_fit_corr{:});
%
%             tmp_vis_Q = (vis_r <= 4 & vis_best_fit_corr > 0.99 & ...
%                 abs(vis_link_ori_vec(:, 3)) < sqrt(1/2));
%             % tmp_vis_Q = ~isnan(best_fit_corr);
%             fig_hdl = figure;
%             fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1] .* 1.5;
%             ax_hdl_1 = subplot(1,3,1);
%             tmp_data = vis_best_fit_FWHM(tmp_vis_Q, 1);
%             histogram(ax_hdl_1, tmp_data, psf_est_int.psf_FWHM_xy_list - 0.05, 'Normalization', 'pdf');
%             legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)));
%             ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
%             ax_hdl_1.YLabel.String = 'PDF';
%             ax_hdl_1.FontSize = 14;
%             ax_hdl_2 = subplot(1,3,2);
%             tmp_data = vis_best_fit_FWHM(tmp_vis_Q, 3);
%             histogram(ax_hdl_2, tmp_data, [0, psf_est_int.psf_FWHM_z_list], 'Normalization', 'pdf');
%             legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)), 'Location', 'northwest');
%             ax_hdl_2.XLabel.String = 'FWHM_{z} (\mum)';
%             ax_hdl_2.YLabel.String = 'PDF';
%             ax_hdl_2.FontSize = 14;
%
%             ax_hdl_3 = subplot(1,3,3);
%             histogram(ax_hdl_3, vis_r(tmp_vis_Q), 'Normalization', 'pdf');
%             legend(ax_hdl_3, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(vis_r(tmp_vis_Q))), 'Location', 'northeast');
%             ax_hdl_3.XLabel.String = 'Estimated radius (\mum)';
%             ax_hdl_3.YLabel.String = 'PDF';
%             ax_hdl_3.FontSize = 14;
%
%             if medfilt_tile_image_Q
%                 fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation', image_group, sprintf('%s_%s_%s_ROI_%d_PSF_FWHM_estimation_hists.png', ...
%                     dataset_name, stack, image_group, iter_roi));
%             else
%                 fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'PSF_estimation_medfilt', image_group, sprintf('%s_%s_%s_ROI_%d_PSF_FWHM_estimation_hists.png', ...
%                     dataset_name, stack, image_group, iter_roi));
%             end
%             fun_print_image_in_several_formats(fig_hdl, fig_fp);
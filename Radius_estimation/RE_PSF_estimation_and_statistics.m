set_env;
DataManager = FileManager;
%% Load psf simulation data
psf_est_int = DataManager.load_data(DataManager.fp_metadata_file('DKLab', 'PSF_simulation_20200103', 'psf_fitting_data'));
% psf_est_int = load('psf_fitting_data_v2.mat');
psf_FWHM_array = cat(1, psf_est_int.psf_FWHM{:});
initial_est_PSF_FWHM_xy = 0.8;
initial_est_PSF_FWHM_z = 8;
[~, initial_est_PSF_FWHM_xy_ind] = min(abs(psf_est_int.psf_FWHM_xy_list - initial_est_PSF_FWHM_xy));
[~, initial_est_PSF_FWHM_z_ind] = min(abs(psf_est_int.psf_FWHM_z_list - initial_est_PSF_FWHM_z));
initial_est_PSF_list_ind = sub2ind(size(psf_est_int.psf_FWHM), initial_est_PSF_FWHM_xy_ind, initial_est_PSF_FWHM_z_ind);
num_psf = numel(psf_est_int.psf_FWHM);
for iter_psf = 1 : num_psf
    tmp_psf = psf_est_int.edge_int_str{iter_psf};
    tmp_num_dist = numel(tmp_psf.normalized_radial_int_dist);
    tmp_psf.normalized_radial_int_dist_itp = cell(size(tmp_psf.normalized_radial_int_dist));
    tmp_psf.normalized_vertical_radial_int_dist_itp = cell(size(tmp_psf.normalized_vertical_radial_int_dist));
    for iter_dist = 1 : tmp_num_dist
        tmp_int_xy_dist = tmp_psf.normalized_radial_int_dist{iter_dist};
        tmp_psf.normalized_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_xy_dist(:, 1), ...
            tmp_int_xy_dist(:, 2), 'linear', 'nearest');
        
        tmp_int_z_dist = tmp_psf.normalized_vertical_radial_int_dist{iter_dist};
        tmp_psf.normalized_vertical_radial_int_dist_itp{iter_dist} = griddedInterpolant(tmp_int_z_dist(:, 1), ...
            tmp_int_z_dist(:, 2), 'linear', 'nearest');
    end
    psf_est_int.edge_int_str{iter_psf} = tmp_psf;
end
fprintf('Finish pre-computing the intensity interpolation\n');
%% Set stack information
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
% stack = 'ML20200201';
grid_version = '240_cube';
skel_version = sprintf('%s_auto', grid_version);
mask_version = sprintf('%s_recon', grid_version);
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
data_octree = grid_info.octree;
tile_size = data_octree.block_size;
grid_voxel_size = grid_info.voxel_size_um;

im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
    'PSF_estimation_wo_medfilt');
if ~isfolder(im_save_folder)
    mkdir(im_save_folder);
end
%% load 240-cube data
% fun_vis_grid(dataset_name, stack, grid_version);
psf_est_info = cell(2, 0);
switch stack
    case 'ML_2018_08_15'
        psf_est_info(:, end+1) = {'Neocortex', [13266 24188 24289 24664 24843]};
        psf_est_info(:, end+1) = {'Hippocampus', [24181 24642]};
        psf_est_info(:, end+1) = {'Superior colliculus', [36303 35398]};
        psf_est_info(:, end+1) = {'Brainstem', [36466, 40314]};
        psf_est_info(:, end+1) = {'Cerebellum', [43779 45024]};
        medfilt_raw_Q = false;
    case 'ML20190124'
        psf_est_info(:, end+1) = {'Neocortex', [18284 19506 12801 17268 12885 22635 14099 16319 23742]};
        psf_est_info(:, end+1) = {'Hippocampus', [27203 27315 32198 33236]};
        psf_est_info(:, end+1) = {'Superior colliculus', [33731 30761 35590]};
        psf_est_info(:, end+1) = {'Brainstem', [33872 37584 38396 35709]};
        psf_est_info(:, end+1) = {'Cerebellum', [42508 44778 42559 43368 44898]};
        medfilt_raw_Q = false;
    case 'ML20200201'
        psf_est_info(:, end+1) = {'Neocortex', [15135 10933 15434 10250 13355]};
        psf_est_info(:, end+1) = {'Hippocampus', [18522 21269 23831 25727 32693]};
        psf_est_info(:, end+1) = {'Brainstem', [35191 36193 37148 38053 45187 46548]};
        psf_est_info(:, end+1) = {'Thalamus', [15208 15058 15093 15305]};
        psf_est_info(:, end+1) = {'Piriform cortex', [15677 14518 13601 13625]};
        medfilt_raw_Q = true;
end
num_region = size(psf_est_info, 2);
compute_region_info = cell(num_region, 1);
for iter_region = 1 : num_region
    tmp_cell = psf_est_info(:, iter_region);
    tmp_cell_1 = arrayfun(@(x) {tmp_cell{1}, x}, tmp_cell{2}', 'UniformOutput', false);
    compute_region_info{iter_region} = cat(1, tmp_cell_1{:});
end
compute_region_info = cat(1, compute_region_info{:});
compute_region_info = compute_region_info';
num_cube = size(compute_region_info, 2);
psf_est_result = cell(num_cube, 1);
%% Radius - PSF joint estimation
parpool_hdl = parpool(12);
parfor iter_cube = 1 : num_cube
    maxNumCompThreads(5);
    tmp_im_save_folder = fullfile(im_save_folder, compute_region_info{1, iter_cube});
    test_grid_ind = compute_region_info{2, iter_cube};
    psf_stat = fun_radius_estimation_get_PSF_stat_in_240_cube(psf_est_int, grid_info, skel_version, test_grid_ind, medfilt_raw_Q);
    psf_stat.region_name = compute_region_info{1, iter_cube};
    psf_stat.cube_label = test_grid_ind;
    %% PSF estimation statistics
    % Histogram of FWHM for:
    % 1. Small vessel
    % 2. Small z-component
    if isfield(grid_info.octree, 'clearing_linear_correction_factor')
        fprintf('Use the acquisition voxel size for PSF size estimation\n');
        psf_stat.best_fit_PSF_FWHM = psf_stat.best_fit_PSF_FWHM ./ grid_info.octree.clearing_linear_correction_factor;
        psf_stat.clearing_linear_correction_factor = grid_info.octree.clearing_linear_correction_factor;
    end
    tmp_vis_Q = (psf_stat.refine_r_est <= 4 & psf_stat.best_fit_corr > 0.99 & ...
        abs(psf_stat.link_ori_vec(:, 3)) < sqrt(1/2));
    % tmp_vis_Q = ~isnan(best_fit_corr);
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1] .* 1.5;
    ax_hdl_1 = subplot(1,2,1);
    tmp_data = psf_stat.best_fit_PSF_FWHM(tmp_vis_Q, 1);
    histogram(ax_hdl_1, tmp_data, [0, psf_est_int.psf_FWHM_xy_list], 'Normalization', 'pdf');
    legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)));
    ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
    ax_hdl_1.YLabel.String = 'PDF';
    ax_hdl_1.FontSize = 14;
    ax_hdl_1.XLim = [0, 2];
    ax_hdl_2 = subplot(1,2,2);
    tmp_data = psf_stat.best_fit_PSF_FWHM(tmp_vis_Q, 3);
    histogram(ax_hdl_2, tmp_data, [0, psf_est_int.psf_FWHM_z_list], 'Normalization', 'pdf');
    legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)), 'Location', 'northwest');
    ax_hdl_2.XLabel.String = 'FWHM_{z} (\mum)';
    ax_hdl_2.YLabel.String = 'PDF';
    ax_hdl_2.FontSize = 14;
    ax_hdl_2.XLim = [0, 12];
    fig_fp = fullfile(tmp_im_save_folder, sprintf('%s_%s_%s_%d_PSF_FWHM_estimation_hists.png', ...
        dataset_name, stack, grid_version, test_grid_ind));
    fun_print_image_in_several_formats(fig_hdl, fig_fp);
    data_fp = fullfile(tmp_im_save_folder, sprintf('%s_%s_%s_%d_PSF_estimation_data.mat', ...
        dataset_name, stack, grid_version, test_grid_ind));
    DataManager.write_data(data_fp, psf_stat);
    psf_est_result{iter_cube} = psf_stat;
end
delete(parpool_hdl);
%% Reploting
% for iter_region = 1 : num_region
%     tmp_im_save_folder = fullfile(im_save_folder, psf_est_info{1, iter_region});
%     tmp_cube_ind_list = psf_est_info{2, iter_region};
%     tmp_cube_est = psf_est_result{iter_region};
%     for iter_cube = 1 : numel(tmp_cube_ind_list)
%         test_grid_ind = tmp_cube_ind_list(iter_cube);
%         % test_grid_ind = 28187;
%         psf_stat = tmp_cube_est{iter_cube};
%         %% PSF estimation statistics
%         % Histogram of FWHM for:
%         % 1. Small vessel
%         % 2. Small z-component
%         tmp_vis_Q = (psf_stat.refine_r_est <= 4 & psf_stat.best_fit_corr > 0.99 & ...
%             abs(psf_stat.link_ori_vec(:, 3)) < sqrt(1/2));
%         % tmp_vis_Q = ~isnan(best_fit_corr);
%         fig_hdl = figure;
%         fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1] .* 1.5;
%         ax_hdl_1 = subplot(1,2,1);
%         tmp_data = psf_stat.best_fit_PSF_FWHM(tmp_vis_Q, 1);
%         histogram(ax_hdl_1, tmp_data, psf_est_int.psf_FWHM_xy_list - 0.05, 'Normalization', 'pdf');
%         legend(ax_hdl_1, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)));
%         ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
%         ax_hdl_1.YLabel.String = 'PDF';
%         ax_hdl_1.FontSize = 14;
%         ax_hdl_2 = subplot(1,2,2);
%         tmp_data = psf_stat.best_fit_PSF_FWHM(tmp_vis_Q, 3);
%         histogram(ax_hdl_2, tmp_data, [0, psf_est_int.psf_FWHM_z_list], 'Normalization', 'pdf');
%         legend(ax_hdl_2, fun_analysis_basic_stat_str_to_string(fun_analysis_get_basic_statistics(tmp_data)), 'Location', 'northwest');
%         ax_hdl_2.XLabel.String = 'FWHM_{z} (\mum)';
%         ax_hdl_2.YLabel.String = 'PDF';
%         ax_hdl_2.FontSize = 14;
%         fig_fp = fullfile(tmp_im_save_folder, sprintf('%s_%s_%s_%d_PSF_FWHM_estimation_hist_pdf.png', ...
%             dataset_name, stack, grid_version, test_grid_ind));
%         fun_print_image_in_several_formats(fig_hdl, fig_fp);
%         delete(fig_hdl);
%     end
% end
%% Reload the data from hard drive
% psf_est_result = cell(num_cube, 1);
for iter_cube = 1 : num_cube
    tmp_im_save_folder = fullfile(im_save_folder, compute_region_info{1, iter_cube});
    test_grid_ind = compute_region_info{2, iter_cube};
    data_fp = fullfile(tmp_im_save_folder, sprintf('%s_%s_%s_%d_PSF_estimation_data.mat', ...
            dataset_name, stack, grid_version, test_grid_ind));
    psf_stat = DataManager.load_data(data_fp);
    psf_est_result{iter_cube} = psf_stat;
end
%%
wb_psf_est_result = psf_est_result;
[wb_psf_stat_xy.mean, wb_psf_stat_xy.median, wb_psf_stat_xy.std] = deal(nan(size(wb_psf_est_result)));
[wb_psf_stat_z.mean, wb_psf_stat_z.median, wb_psf_stat_z.std] = deal(nan(size(wb_psf_est_result)));
for iter_region = 1 : numel(wb_psf_est_result)
    tmp_data = wb_psf_est_result{iter_region};
    tmp_vis_Q = (tmp_data.refine_r_est <= 3 & tmp_data.best_fit_corr > 0.99 & ...
        abs(tmp_data.link_ori_vec(:, 3)) < sqrt(2)/2);
    tmp_xy_stat = fun_analysis_get_basic_statistics(tmp_data.best_fit_PSF_FWHM(tmp_vis_Q, 1));
    wb_psf_stat_xy.mean(iter_region) = tmp_xy_stat.mean;
    wb_psf_stat_xy.median(iter_region) = tmp_xy_stat.median;
    wb_psf_stat_xy.std(iter_region) = tmp_xy_stat.std;
    
    tmp_z_stat = fun_analysis_get_basic_statistics(tmp_data.best_fit_PSF_FWHM(tmp_vis_Q, 3));
    wb_psf_stat_z.mean(iter_region) = tmp_z_stat.mean;
    wb_psf_stat_z.median(iter_region) = tmp_z_stat.median;
    wb_psf_stat_z.std(iter_region) = tmp_z_stat.std;
end
fprintf('FWHM xy: %.2f +/- %.2f\n', mean(wb_psf_stat_xy.mean), std(wb_psf_stat_xy.mean));
fprintf('FWHM z: %.2f +/- %.2f\n', mean(wb_psf_stat_z.mean), std(wb_psf_stat_z.mean));
%% Save
psf_estimation_str = struct;
psf_estimation_str.psf_est_setting = psf_est_info;
psf_estimation_str.result = psf_est_result;
psf_estimation_str.PSF_FWHM_stat_xy = wb_psf_stat_xy;
psf_estimation_str.PSF_FWHM_stat_z = wb_psf_stat_z;

psf_estimation_str.all_avg_PSF_FWHM_xy.mean = mean(wb_psf_stat_xy.mean);
psf_estimation_str.all_avg_PSF_FWHM_xy.std = std(wb_psf_stat_xy.mean);
psf_estimation_str.all_avg_PSF_FWHM_z.mean = mean(wb_psf_stat_z.mean);
psf_estimation_str.all_avg_PSF_FWHM_z.std = std(wb_psf_stat_z.mean);
% This dataset is wirted. Select the valid PSF for simulation: 
valid_psf_Q = wb_psf_stat_xy.mean < 1;
psf_estimation_str.avg_PSF_FWHM_xy.mean = mean(wb_psf_stat_xy.mean(valid_psf_Q));
psf_estimation_str.avg_PSF_FWHM_xy.std = std(wb_psf_stat_xy.mean(valid_psf_Q));
psf_estimation_str.avg_PSF_FWHM_z.mean = mean(wb_psf_stat_z.mean(valid_psf_Q));
psf_estimation_str.avg_PSF_FWHM_z.std = std(wb_psf_stat_z.mean(valid_psf_Q));

if isfield(grid_info.octree, 'clearing_linear_correction_factor')
    psf_estimation_str.clearing_linear_correction_factor = grid_info.octree.clearing_linear_correction_factor;
    psf_estimation_str.avg_PSF_FWHM_xy.mean_c = psf_estimation_str.avg_PSF_FWHM_xy.mean * psf_estimation_str.clearing_linear_correction_factor;
    psf_estimation_str.avg_PSF_FWHM_z.mean_c = psf_estimation_str.avg_PSF_FWHM_z.mean * psf_estimation_str.clearing_linear_correction_factor;
end

psf_estimation_str.dataset_name = dataset_name;
psf_estimation_str.stack = stack;
psf_estimation_str.filepath = DataManager.fp_metadata_file(dataset_name, ...
    stack, 'PSF_estimation_v2_w_medfilt3');
save(psf_estimation_str.filepath, '-struct', 'psf_estimation_str');
%% Generate PSF Look up table at higher resolution.
% Only care about the xy threshold
% gpuDevice(2);
% vessel_radius_list = [0.3 : 0.1 : 15];
% z_component_list = 0 : 0.05 : 1;
% theta_list = asin(z_component_list);
% if isfield(psf_estimation_str, 'clearing_linear_correction_factor')
%     psf_sigma = [ones(1,2) .* psf_estimation_str.avg_PSF_FWHM_xy.mean_c, ...
%         psf_estimation_str.avg_PSF_FWHM_z.mean_c] ./ (2 * sqrt(2 * log(2)));    
% else
%     psf_sigma = [ones(1,2) .* psf_estimation_str.avg_PSF_FWHM_xy.mean, ...
%         psf_estimation_str.avg_PSF_FWHM_z.mean] ./ (2 * sqrt(2 * log(2)));
% end
% dx_min_scale_ratio = 0.1;
% tmp_tic = tic;
% edge_intensity = fun_radius_estimation_compute_vessel_edge_int_given_psf(psf_sigma, ...
%     vessel_radius_list, theta_list, dx_min_scale_ratio, false);
% fprintf('Finish computation. Elapsed time is %f seconds.\n', toc(tmp_tic));
% edge_intensity.vsl_ori_vec_z_list = z_component_list;
% edge_intensity.dataset_name = dataset_name;
% edge_intensity.stack = stack;
% edge_intensity.filepath = DataManager.fp_metadata_file(dataset_name, ...
%     stack, 'Estimated_PSF_edge_intensity');
% save(edge_intensity.filepath, '-struct', 'edge_intensity');
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
% stack = 'ML_2019_01_24';
mask_version = '240_cube_recon';
skel_version = '240_cube_re';
grid_version = '240_cube';
grid_info = DataManager.load_grid(dataset_name, stack, grid_version);
im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
    'whole_brain_stat', 'internal_cubes', 'anisotropy');
registration_name = 'Allen_2017_25um_nonrigid.mat';
registration_str = DataManager.load_registration_data(dataset_name, stack, registration_name);
%% Find the 240-cubes that are completely inside the brain mask
wb_mask_ds_ratio = 16;
wb_mask = DataManager.load_data(sprintf('%s_mask.nii.gz', fullfile(DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration'),...
    sprintf('%s_%s_d16x_registration', dataset_name, stack)))) > 0;
is_internal_240_cube_ratio = fun_analysis_get_bbox_in_mask_vol_ratio(wb_mask, grid_info.bbox_xyz_mmxx_list ./ wb_mask_ds_ratio);
is_internal_240_cube_Q = is_internal_240_cube_ratio == 1;
internal_cube_ind = grid_info.bbox_grid_ind_list(is_internal_240_cube_ratio == 1 );
%% Load whole brain 240-cube statistics
tic
wb_240_cube_stat_str = DataManager.load_analysis_data(dataset_name, stack, sprintf('%s_%s_%s_240_cube_stat_data.mat', ...
    dataset_name, stack, reconstruction_name), wb_stat_folder_name);
fprintf('Finish loading 240 cube data\n');
toc

is_valid_Q = (is_internal_240_cube_ratio == 1);
%% Correlation 
% end to end anisotropy and volume-weighted anisotropy, for capillary-only
% and all the vessels
candidate_stat = {wholebrain_stat_str.anisotropy_all_vw_fa, ...
    wholebrain_stat_str.anisotropy_all_vw_min2max_z, local_recon_stat.all_link_isotropy, ...
    wholebrain_stat_str.anisotropy_all_vw_svd1_z, local_recon_stat.all_link_ori_svd_1, ...
    wholebrain_stat_str.anisotropy_cap_vw_fa, ...
    wholebrain_stat_str.anisotropy_cap_vw_min2max_z, local_recon_stat.capillary_isotropy, ...
    wholebrain_stat_str.anisotropy_cap_vw_svd1_z, local_recon_stat.capillary_ori_svd_1};
candidate_name = {'FA_VW_All', ...
    'MXZ_VW_All', 'MXZ_All', ...
    'V1Z_VW_All',  'V1Z_All', ...
    'FA_VW_Cap', ...
    'MXZ_VW_Cap', 'MXZ_Cap', ...
    'V1Z_VW_Cap', 'V1Z_Cap'};
num_cand = numel(candidate_stat);
corr_mat = eye(num_cand);
for iter_1 = 1 : num_cand
    for iter_2 = (iter_1 + 1) : num_cand 
        tmp_X = candidate_stat{iter_1};
        tmp_Y = candidate_stat{iter_2};
        tmp_is_valid = ~isnan(tmp_X) & ~isnan(tmp_Y) & is_valid_Q;
        tmp_X = tmp_X(tmp_is_valid);
        tmp_Y = tmp_Y(tmp_is_valid);
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
colormap jet
daspect([1,1,1]);
ax.XAxis.TickValues = 1 : num_cand;
ax.XAxis.TickLabelRotation = 0;
ax.XAxis.TickLabels = candidate_name;
ax.YAxis.TickValues = 1 : num_cand;
ax.YAxis.TickLabels = candidate_name;
ax.YAxis.TickLabelInterpreter = 'none';
ax.XAxis.TickLabelInterpreter = 'none';
ax.FontSize = 14;
ax.FontWeight = 'bold';
ax.Title.String = 'Local anisotropy measure correlation matrix';
ax.Title.FontSize = 18;
tmp_fp = DataManager.fp_visualization_image(dataset_name, stack, 'Local_anisotropy_measure_correlation_matrix.png', 'Paper');
fun_print_image_in_several_formats(fig_hdl, tmp_fp);
% export_fig(tmp_fp, '-png', '-eps', fig_hdl);
%% Plot the 2D histogram for showing the correlation between the isotropy measures
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* 2;
ax_hdl = axes(fig_hdl);
tmp_x = wb_240_cube_stat_str.wb_ai_cap_vw.fa_z(is_valid_Q);
tmp_y = wb_240_cube_stat_str.wb_ai_cap_vw.fa_p(is_valid_Q);
tmp_x_edge = -5 : 1 : 20;
tmp_y_edge = [0, 10 .^ (-5 : 0.5 : 0)];
tmp_y(tmp_y == 0) = 5e-5;
histogram2(ax_hdl, tmp_x, tmp_y, tmp_x_edge, tmp_y_edge, 'DisplayStyle', 'tile');
ax_hdl.XLabel.String = 'Capillary FA_z';
ax_hdl.YLabel.String = 'Capillary FA_p';
ax_hdl.YScale = 'log';
ax_hdl.ColorScale = 'log';
ax_hdl.FontSize = 14;
ax_hdl.FontWeight = 'bold';
tmp_fp = DataManager.fp_visualization_image(dataset_name, stack, 'Local_vw_isotropy_vs_isotropy_hist2.png', 'Paper');
% fun_print_image_in_several_formats(fig_hdl, tmp_fp);

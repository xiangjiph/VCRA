set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
skl_grid_name = '240_cube_auto';
mask_name = '240_cube_recon';
output_graph_grid_name = sprintf('%s_analyzed', skl_grid_name);
image_grid_name = '240_cube';

grid_c_version = '240_cube_combined_5_o_2';
grid_c_info = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c_info.grid_ori;

layer_list = 1 : grid_info.num_grid_layer;
num_layer_to_process = numel(layer_list);
internal_offset = 16;
%% Load whole brain reconstruction analysis
% Load all the computed features and max projection of the reconstructed images to 3 directions. Takes a few minutes
wholebrain_stat_str = fun_analysis_load_whole_brain_features(grid_info, mask_name);
%% Compute the regional anisotropy ( unweighted and weighted by volume ) 
anisotropy_stat_all_uw = fun_analysis_get_anisotropy_stat_from_link_feature_cell(wholebrain_stat_str.link_features, inf, false);
anisotropy_stat_cap_uw = fun_analysis_get_anisotropy_stat_from_link_feature_cell(wholebrain_stat_str.link_features, 3.5, false);
% The volume - weighted version takes too long if not distributed. 
tic
anisotropy_stat_all_w = fun_analysis_get_anisotropy_stat_from_link_feature_cell(wholebrain_stat_str.link_features, inf, true);
toc
tic
anisotropy_stat_cap_w = fun_analysis_get_anisotropy_stat_from_link_feature_cell(wholebrain_stat_str.link_features, 3.5, true);
toc
%% 
region_anisotropy_data = struct;
region_anisotropy_data.dataset_name = dataset_name;
region_anisotropy_data.stack = stack;
region_anisotropy_data.grid_name = image_grid_name;
region_anisotropy_data.capillary_max_raidus = 3.5;
region_anisotropy_data.all_unweighted = anisotropy_stat_all_uw;
region_anisotropy_data.all_volume_weighted = anisotropy_stat_all_w;
region_anisotropy_data.capillary_unweighted = anisotropy_stat_cap_uw;
region_anisotropy_data.capillary_volume_weighted = anisotropy_stat_cap_w;
tmp_fp = fullfile(DataManager.fp_metadata_folder(dataset_name, stack), sprintf('%s_%s_cube_anistropy_data.mat', dataset_name, stack));
save(tmp_fp, '-struct', 'region_anisotropy_data');
%% Correlation between vessel anisotropy and FA average 
mri_voxel_size = [43, 43, 43];
avg_DT_FA = nan(grid_info.grid_size);
avg_DT_WM = nan(grid_info.grid_size);
avg_DT_AD = nan(grid_info.grid_size);

[avg_DT_vec1, avg_DT_vec2, avg_DT_vec3] = deal(nan(grid_info.grid_size));
sample_fun_hld = @mean;
for iter_bbox = 1 : grid_info.num_valid_cube
    tmp_grid_ind = grid_info.bbox_grid_ind_list(iter_bbox);
    tmp_bbox_mmll = grid_info.bbox_xyz_mmll_list(iter_bbox, :);
    % Convert to the coordinate of DTMRI data
    tmp_dt_bbox_mmll = ceil(tmp_bbox_mmll ./ [mri_voxel_size, mri_voxel_size]);
    tmp_dt_FA = crop_bbox3(im_FA_moved, tmp_dt_bbox_mmll);
    tmp_dt_FA = tmp_dt_FA(tmp_dt_FA > 0);
    if any(isfinite(tmp_dt_FA))
        avg_DT_FA(tmp_grid_ind) = sample_fun_hld(tmp_dt_FA);
    end
    
    tmp_dt_WM = crop_bbox3(im_WM_moved, tmp_dt_bbox_mmll);
    tmp_dt_WM = tmp_dt_WM(tmp_dt_WM > 0);
    if any(isfinite(tmp_dt_WM))
        avg_DT_WM(tmp_grid_ind) = sample_fun_hld(tmp_dt_WM);
    end
    
    tmp_dt_AD = crop_bbox3(im_AD_moved, tmp_dt_bbox_mmll);
    tmp_dt_AD = tmp_dt_AD(tmp_dt_AD > 0);
    if any(isfinite(tmp_dt_AD))
        avg_DT_AD(tmp_grid_ind) = sample_fun_hld(tmp_dt_AD);
    end
    
    tmp_vec_1 = crop_bbox3(im_vec_1, tmp_dt_bbox_mmll);
    tmp_vec_2 = crop_bbox3(im_vec_2, tmp_dt_bbox_mmll);
    tmp_vec_3 = crop_bbox3(im_vec_3, tmp_dt_bbox_mmll);
    
    tmp_vec_1 = tmp_vec_1(:);
    tmp_vec_2 = tmp_vec_2(:);
    tmp_vec_3 = tmp_vec_3(:);
    
    tmp_anisotropy = fun_analysis_get_link_anisotropy(cat(2, tmp_vec_1, tmp_vec_2, tmp_vec_3), true);
    
    if ~isempty(tmp_anisotropy.svd_min2max) && ~isnan(tmp_anisotropy.svd_min2max)    
        avg_DT_vec1(tmp_grid_ind) = tmp_anisotropy.svd_max_vec(1);
        avg_DT_vec2(tmp_grid_ind) = tmp_anisotropy.svd_max_vec(2);
        avg_DT_vec3(tmp_grid_ind) = tmp_anisotropy.svd_max_vec(3);
    end
end
fprintf('Finish extracing region statistics from registration data\n');
%%
test_anisotropy_str = anisotropy_stat_cap_w;
candidate_stat = {test_anisotropy_str.svd_1, test_anisotropy_str.svd_1_z, ...
    test_anisotropy_str.min2max, test_anisotropy_str.min2max_z, ...
    test_anisotropy_str.fractional_anisotropy, test_anisotropy_str.fractional_anisotropy_z, ...
    avg_DT_FA, avg_DT_WM, avg_DT_AD};

candidate_name = {'SVD1', 'SVD1_z', ...
    'SVDm2x', 'SVDm2x_z', ...
    'Vessel FA', 'Vessel FA_z', ...
    'FA', 'WM', 'AD'};

% test_anisotropy_str = anisotropy_stat_cap_w;
% candidate_stat = {test_anisotropy_str.svd_1, test_anisotropy_str.svd_1_z, ...
%     test_anisotropy_str.min2max, test_anisotropy_str.min2max_z, ...
%     avg_DT_FA, avg_DT_WM, avg_DT_AD};
% 
% candidate_name = {'SVD1', 'SVD1_z', ...
%     'SVDm2x', 'SVDm2x_z', ...
%     'FA', 'WM', 'AD'};

num_cand = numel(candidate_stat);
corr_mat = eye(num_cand);
preselect_Q = test_anisotropy_str.svd_1_z > 2;
% preselect_Q = test_anisotropy_str.num_data > 100;
for iter_1 = 1 : num_cand
    for iter_2 = (iter_1 + 1) : num_cand 
        tmp_X = candidate_stat{iter_1};
        tmp_X = tmp_X(preselect_Q);
        [tmp_X, tmp_selected_idx] = fun_analysis_select_data_by_percentile(tmp_X, 5, 95, false);
        tmp_Y = candidate_stat{iter_2};
        tmp_Y = tmp_Y(preselect_Q);
        tmp_Y = tmp_Y(tmp_selected_idx);
        [tmp_Y, tmp_selected_idx] = fun_analysis_select_data_by_percentile(tmp_Y, 5, 95, false);
        tmp_X = tmp_X(tmp_selected_idx);
        % histogram2(tmp_Y, tmp_X, 'DisplayStyle', 'tile')
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
ax.XAxis.TickLabelRotation = 90;
ax.XAxis.TickLabels = candidate_name;
ax.YAxis.TickValues = 1 : num_cand;
ax.YAxis.TickLabels = candidate_name;
ax.FontSize = 14;
ax.FontWeight = 'Bold';
%%
figure;
histogram2(test_anisotropy_str.min2max(preselect_Q), avg_DT_FA(preselect_Q), 'DisplayStyle', 'tile')
%% Local orientation
test_ori_vec = test_anisotropy_str.ori_vec;

test_ori_vec_1 = squeeze(test_ori_vec(1, :, :, :));
test_ori_vec_2 = squeeze(test_ori_vec(2, :, :, :));
test_ori_vec_3 = squeeze(test_ori_vec(3, :, :, :));

tmp_flip_Q = (test_ori_vec_3 < 0);
test_ori_vec_1(tmp_flip_Q) = - test_ori_vec_1(tmp_flip_Q);
test_ori_vec_2(tmp_flip_Q) = - test_ori_vec_2(tmp_flip_Q);
test_ori_vec_3(tmp_flip_Q) = - test_ori_vec_3(tmp_flip_Q);

tmp_flip_Q = (avg_DT_vec3 < 0);
avg_DT_vec1(tmp_flip_Q) = - avg_DT_vec1(tmp_flip_Q);
avg_DT_vec2(tmp_flip_Q) = - avg_DT_vec2(tmp_flip_Q);
avg_DT_vec3(tmp_flip_Q) = - avg_DT_vec3(tmp_flip_Q);

% test_ori_vec_1 = abs(test_ori_vec_1);
% test_ori_vec_2 = abs(test_ori_vec_2);
% test_ori_vec_3 = abs(test_ori_vec_3);
% 
% avg_DT_vec1 = abs(avg_DT_vec1);
% avg_DT_vec2 = abs(avg_DT_vec2);
% avg_DT_vec3 = abs(avg_DT_vec3);

% Compute the inner product with DTMRI data
cos_vessel_DTMRI = test_ori_vec_1 .* avg_DT_vec1 + ...
    test_ori_vec_2 .* avg_DT_vec2 + test_ori_vec_3 .* avg_DT_vec3;

histogram(cos_vessel_DTMRI(preselect_Q))
%%
% fig_hdl = figure;
% subplot(1,3,1);
% histogram2(test_ori_vec_1(preselect_Q), avg_DT_vec1(preselect_Q), 'DisplayStyle', 'tile');
% subplot(1,3,2);
% histogram2(test_ori_vec_2(preselect_Q), avg_DT_vec2(preselect_Q), 'DisplayStyle', 'tile');
% subplot(1,3,3);
% histogram2(test_ori_vec_3(preselect_Q), avg_DT_vec3(preselect_Q), 'DisplayStyle', 'tile');
%%
fig_hld = figure;
subplot(2,3,1);
histogram(test_ori_vec_1);
subplot(2,3,2);
histogram(test_ori_vec_2);
subplot(2,3,3);
histogram(test_ori_vec_3);
subplot(2,3,4);
histogram(avg_DT_vec1);
subplot(2,3,5);
histogram(avg_DT_vec2);
subplot(2,3,6);
histogram(avg_DT_vec3);
%% Remove the region with few vessels
% fig_hdl = figure;
% 
% tmp_selected_Q = anisotropy_stat_cap_uw.num_data > 150;
% tmp_X = anisotropy_stat_cap_uw.svd_1_z;
% 
% tmp_X = tmp_X(tmp_selected_Q);
% 
% tmp_Y = avg_DT_FA;
% tmp_Y = tmp_Y(tmp_selected_Q);
% tmp_corr = corrcoef(tmp_X, tmp_Y);
% ax_1 = subplot(1,3,1);
% histogram2(ax_1, tmp_X, tmp_Y, 'DisplayStyle', 'tile')
% ax_1.YLabel.String = 'FA';
% 
% tmp_Y = avg_DT_WM;
% tmp_Y = tmp_Y(tmp_selected_Q);
% ax_2 = subplot(1,3,2);
% histogram2(ax_2, tmp_X, tmp_Y, 'DisplayStyle', 'tile')
% ax_2.YLabel.String = 'WM';
% 
% tmp_Y = avg_DT_AD;
% tmp_Y = tmp_Y(tmp_selected_Q);
% ax_3 = subplot(1,3,3);
% histogram2(ax_3, tmp_X, tmp_Y, 'DisplayStyle', 'tile')
% ax_3.YLabel.String = 'AD';


%% Correlation between vessel eigenvector direction and DRMRI direction 
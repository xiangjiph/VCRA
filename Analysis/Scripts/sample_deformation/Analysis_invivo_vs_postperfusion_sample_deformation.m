DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML20190124';
wb_im_fp = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
wb_im = DataManager.load_single_tiff(wb_im_fp);
% fun_vis_grid(dataset_name, stack, '240_cube');

wb_im_horizontal = permute(wb_im, [3, 2, 1]);
%% Moving max proj
% wb_im_h_movmax = movmax(wb_im_horizontal, 5, 3);
% implay(wb_im_h_movmax);
region_sub = [72, 512, 455];
region_sub_1um = region_sub .* 16;
%% Determine the covering 1072 cubes
grid_version = '240_cube';
grid_c = DataManager.load_grid(dataset_name, stack, '240_cube_combined_5_o_1');
grid_c_label = find(all(grid_c.bbox_xyz_mmll_pixel_list(:, 1:3) < region_sub_1um, 2) & ...
    all(grid_c.bbox_xyz_mmxx_pixel_list(:, 4:6) > region_sub_1um, 2));
grid_c_label = grid_c_label(1);
%%
grid_c_bbox_sub = grid_c.bbox_grid_sub_list(grid_c_label, :);
grid_c_bbox_mmxx_grid = grid_c.bbox_xyz_mmxx_grid_list(grid_c_label, :);

vessel_image = DataManager.load_blocks_files('image', dataset_name, stack, grid_version, ...
    grid_c_bbox_mmxx_grid(1):grid_c_bbox_mmxx_grid(4), ...
    grid_c_bbox_mmxx_grid(2):grid_c_bbox_mmxx_grid(5), ...
    grid_c_bbox_mmxx_grid(3):grid_c_bbox_mmxx_grid(6), 'uint16');

% vessel_mask_0 = DataManager.load_blocks_files('mask', dataset_name, stack, grid_version, ...
%     grid_c_bbox_mmxx_grid(1):grid_c_bbox_mmxx_grid(4), ...
%     grid_c_bbox_mmxx_grid(2):grid_c_bbox_mmxx_grid(5), ...
%     grid_c_bbox_mmxx_grid(3):grid_c_bbox_mmxx_grid(6), 'logical');

vessel_mask = DataManager.load_blocks_files('mask', dataset_name, stack, grid_version, ...
    grid_c_bbox_mmxx_grid(1):grid_c_bbox_mmxx_grid(4), ...
    grid_c_bbox_mmxx_grid(2):grid_c_bbox_mmxx_grid(5), ...
    grid_c_bbox_mmxx_grid(3):grid_c_bbox_mmxx_grid(6), 'logical', '240_cube_recon_v2');

vessel_skel = DataManager.load_blocks_files('skel', dataset_name, stack, grid_version, ...
    grid_c_bbox_mmxx_grid(1):grid_c_bbox_mmxx_grid(4), ...
    grid_c_bbox_mmxx_grid(2):grid_c_bbox_mmxx_grid(5), ...
    grid_c_bbox_mmxx_grid(3):grid_c_bbox_mmxx_grid(6), 'single', '240_cube_re');

% DataManager.visualize_itksnap(vessel_image, vessel_mask)
%% Local the invivo data 
dataset_name_iv = 'DKLab';
stack_iv = 'Rui20180622_perfusion';
grid_info_iv = DataManager.load_grid(dataset_name_iv, stack_iv, grid_version);
invivo_im = DataManager.load_blocks_files('data', dataset_name_iv, stack_iv, ...
    grid_version, 1 : grid_info_iv.grid_size(1), 1 : grid_info_iv.grid_size(2), ...
    1 : grid_info_iv.grid_size(3));
vg_iv = DataManager.load_graph_in_block(dataset_name_iv, stack_iv, 'merged', 0, 0, 0);
%% Find the bounding box of the in vivo imaged volume in whole brain coordinate 
% Crop the skeleton and try rigid registration 
crop_center = [490 610 593];
crop_radius = 350;
crop_bbox_min = crop_center - [300, 350, 300];
crop_bbox_max = crop_center + [350, 350, 350];
crop_bbox_mmll = [crop_bbox_min, crop_bbox_max - crop_bbox_min + 1];

vessel_mask_cropped = crop_bbox3(vessel_mask, crop_bbox_mmll);
% volumeViewer(vessel_mask_cropped)
vessel_skel_cropped = crop_bbox3(vessel_skel, crop_bbox_mmll);
vg_pf = fun_skeleton_to_graph(vessel_skel_cropped);
vg_pf = fun_graph_add_radius(vg_pf, vessel_skel_cropped);
%% Compute the radius of the vessel segment
vg_pf.link.features = fun_analysis_get_link_features(vg_pf);
vg_iv.link.features = fun_analysis_get_link_features(vg_iv);
%% 
% Use the penetrating vessel to align first
pf_selected_Q = vg_pf.link.features.dt_median >= 3;
pf_selected_ind = cat(1, vg_pf.link.cc_ind{pf_selected_Q});
pf_selected_label = full(vg_pf.link.map_ind_2_label(pf_selected_ind));

iv_selected_Q = vg_iv.link.features.dt_median >= 3.5;
iv_selected_ind = cat(1, vg_iv.link.cc_ind{iv_selected_Q});
iv_selected_label = full(vg_iv.link.map_ind_2_label(iv_selected_ind));

downsample_step = 2;
skel_1_block_size = vg_iv.num.mask_size;
skel_2_block_size = vg_pf.num.mask_size;

skel_1_sub = fun_ind2sub(skel_1_block_size, iv_selected_ind(1:downsample_step:end));
% Correct the voxel size on x-y direction 
skel_1_sub(:, 1:2) = skel_1_sub(:, 1:2) .* 656 / 700;
skel_2_sub = fun_ind2sub(skel_2_block_size, pf_selected_ind(1:downsample_step:end));

skel_1_label = iv_selected_label(1:downsample_step:end);
skel_2_label = pf_selected_label(1:downsample_step:end);
%% Visualization
skel_2_com = mean(skel_2_sub, 1);
skel_2_r = skel_2_sub - skel_2_com;
rot_x_agl = -(20 / 180 * pi);
rot_y_agl = -pi/2;
rot_z_agl = 0;
Rx = [1, 0, 0; 0, cos(rot_x_agl), -sin(rot_x_agl); 0, sin(rot_x_agl), cos(rot_x_agl)];
Ry = [cos(rot_y_agl), 0, sin(rot_y_agl); 0, 1, 0; -sin(rot_y_agl), 0, cos(rot_y_agl)];
Rz = [cos(rot_z_agl), -sin(rot_z_agl), 0; sin(rot_z_agl), cos(rot_z_agl), 0; 0, 0, 1];

skel_2_new = [1, 0, 0; 0, -1, 0; 0, 0, 1] * Rx * Ry * skel_2_r';
skel_2_new = skel_2_new' + skel_2_com;
skel_pre_regid_T = [1, 0, 0; 0, -1, 0; 0, 0, 1] * Rx * Ry;
% skel_2_sub = skel_2_sub(:, [3,2,1]);
% Flip x and flip y
% skel_2_sub_max = max(skel_2_sub, [], 1);
% skel_2_sub(:, 1) = skel_2_sub_max(1) - skel_2_sub(:, 1);
% skel_2_sub(:, 2) = skel_2_sub_max(2) - skel_2_sub(:, 2);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
tiledlayout(fig_hdl, 1, 3);
ax_hdl_1 = nexttile;
scatter3(ax_hdl_1, skel_1_sub(:, 1), skel_1_sub(:, 2), skel_1_sub(:, 3));
ax_hdl_1.XLabel.String = 'X';
ax_hdl_1.YLabel.String = 'Y';
ax_hdl_1.ZLabel.String = 'Z';
ax_hdl_1.DataAspectRatio = [1,1,1];
ax_hdl_2 = nexttile;
scatter3(ax_hdl_2, skel_2_sub(:, 1), skel_2_sub(:, 2), skel_2_sub(:, 3));
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.XLabel.String = 'X';
ax_hdl_2.YLabel.String = 'Y';
ax_hdl_2.ZLabel.String = 'Z';
ax_hdl_2 = nexttile;
scatter3(ax_hdl_2, skel_2_new(:, 1), skel_2_new(:, 2), skel_2_new(:, 3));
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.XLabel.String = 'X';
ax_hdl_2.YLabel.String = 'Y';
ax_hdl_2.ZLabel.String = 'Z';
%%
skel_sample_opt = struct;
skel_sample_opt.total_num_descriptor = 1e4;
skel_sample_opt.rm_descriptor_by_pdist2_Q = false;
skel_sample_opt.max_disp_pixel_yxz = [20, 20, 20];

[des_1_sub, des_2_sub_shift, des_1_label, des_2_label] = fun_feature_match_sample_descriptor(skel_1_sub, ...
    skel_2_new, skel_1_label, skel_2_label, skel_sample_opt);
%% Coherent point drift for penetrating vessels 
cpd_opt.method = 'affine';
cpd_opt.rot = 0;
cpd_opt.scale = 1;
% cpd_opt.beta = 6;            % the width of Gaussian kernel (smoothness)
% cpd_opt.lambda = 16;          % regularization weight
cpd_opt.viz = 1;              % show every iteration
cpd_opt.outliers = 0.9;       % use 0.7 noise weight
cpd_opt.max_it = 100;        % Maxinum number of iteration
cpd_opt.fgt = 0;              % do not use FGT (default)
cpd_opt.tol = 1e-4;         % Error tolorance ( how is it defined? )
cpd_opt.normalize = 1;        % normalize to unit variance and zero mean before registering (default)
cpd_opt.corresp = 1;          % compute correspondence vector at the end of registration (not being estimated by default)
% matchparams_af.opt = opt;
% matchparams_af.projectionThr = 5;
% matchparams_af.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
% matchparams_af.debug = false;
% matchparams_af.viz = true;
[Transform, C] = cpd_register(des_1_sub, des_2_sub_shift, cpd_opt);
dist_1_to_2t = pdist2(des_1_sub, Transform.Y);
[mutual_nn_idx_1, mutual_nn_idx_2, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_1_to_2t, false);

X_ = des_1_sub(mutual_nn_idx_1,:);
Y_ = des_2_sub_shift(mutual_nn_idx_2,:);
tY_= Transform.Y(mutual_nn_idx_2,:);
%% Decomposition of Affine matrix
[matQ, matR] = qr(Transform.R);
scaleMatrix = zeros(size(matR));
scaleMatrix([1, 5, 9]) = diag(matR);
shearMat = scaleMatrix \ matR;
%% Transform the entire perfusion dataset according to the affine registration of the penetrating vessels
iv_ind = vg_iv.link.pos_ind;
iv_label = full(vg_iv.link.map_ind_2_label(iv_ind));

pf_ind = vg_pf.link.pos_ind;
pf_label = full(vg_pf.link.map_ind_2_label(pf_ind));

downsample_step = 2;
skel_1_block_size = vg_iv.num.mask_size;
skel_2_block_size = vg_pf.num.mask_size;

skel_1_sub = fun_ind2sub(skel_1_block_size, iv_ind(1:downsample_step:end));
skel_1_sub(:, 1:2) = skel_1_sub(:, 1:2) .* 656 / 700;
skel_2_sub = fun_ind2sub(skel_2_block_size, pf_ind(1:downsample_step:end));

skel_1_label = iv_label(1:downsample_step:end);
skel_2_label = pf_label(1:downsample_step:end);
% Transform skel_2_sub according to the transformation computed from
% penetrating vessels
skel_2_sub_T = (skel_pre_regid_T * (skel_2_sub - skel_2_com)' + skel_2_com')';
skel_2_sub_T = (Transform.R * skel_2_sub_T' + Transform.t)';

skel_sample_opt = struct;
skel_sample_opt.total_num_descriptor = 2e4;
skel_sample_opt.rm_descriptor_by_pdist2_Q = true;
skel_sample_opt.max_disp_pixel_yxz = [15, 15, 15];

[des_1_sub, des_2_sub_shift, des_1_label, des_2_label] = fun_feature_match_sample_descriptor(skel_1_sub, ...
    skel_2_sub_T, skel_1_label, skel_2_label, skel_sample_opt);
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
% tiledlayout(fig_hdl, 1, 2);
% ax_hdl_1 = nexttile;
ax_hdl_1 = axes(fig_hdl);
scatter3(ax_hdl_1, des_1_sub(:, 1), des_1_sub(:, 2), des_1_sub(:, 3));
ax_hdl_1.XLabel.String = 'X';
ax_hdl_1.YLabel.String = 'Y';
ax_hdl_1.ZLabel.String = 'Z';
ax_hdl_1.DataAspectRatio = [1,1,1];
% ax_hdl_2 = nexttile;
hold(ax_hdl_1, 'on');
scatter3(ax_hdl_1, des_2_sub_shift(:, 1), des_2_sub_shift(:, 2), des_2_sub_shift(:, 3));
% ax_hdl_2.DataAspectRatio = [1,1,1];
% ax_hdl_2.XLabel.String = 'X';
% ax_hdl_2.YLabel.String = 'Y';
% ax_hdl_2.ZLabel.String = 'Z';
%%
cpd_opt.method = 'affine';
cpd_opt.rot = 1;
cpd_opt.scale = 1;
% cpd_opt.beta = 6;            % the width of Gaussian kernel (smoothness)
% cpd_opt.lambda = 16;          % regularization weight
cpd_opt.viz = 1;              % show every iteration
cpd_opt.outliers = 0.7;       % use 0.7 noise weight
cpd_opt.max_it = 100;        % Maxinum number of iteration
cpd_opt.fgt = 0;              % do not use FGT (default)
cpd_opt.tol = 1e-5;         % Error tolorance ( how is it defined? )
cpd_opt.normalize = 1;        % normalize to unit variance and zero mean before registering (default)
cpd_opt.corresp = 1;          % compute correspondence vector at the end of registration (not being estimated by default)
% matchparams_af.opt = opt;
% matchparams_af.projectionThr = 5;
% matchparams_af.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
% matchparams_af.debug = false;
% matchparams_af.viz = true;
[Transform_all_skel, C] = cpd_register(des_1_sub, des_2_sub_shift, cpd_opt);
dist_1_to_2t = pdist2(des_1_sub, Transform_all_skel.Y);

[mutual_nn_idx_1, mutual_nn_idx_2, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_1_to_2t, false);

% within_error_dist_Q = mutual_nn_dist < 5;
% mutual_nn_idx_1 = mutual_nn_idx_1(within_error_dist_Q);
% mutual_nn_idx_2 = mutual_nn_idx_2(within_error_dist_Q);

X_ = des_1_sub(mutual_nn_idx_1,:);
Y_ = des_2_sub_shift(mutual_nn_idx_2,:);
tY_= Transform_all_skel.Y(mutual_nn_idx_2,:);
X_label = des_1_label(mutual_nn_idx_1);
Y_label = des_2_label(mutual_nn_idx_2);

%%
[matQ, matR] = qr(Transform_all_skel.R * Transform.R);
matQ = [-1, 0, 0; 0, -1, 0; 0, 0, 1] * matQ;
matR = [-1, 0, 0; 0, -1, 0; 0, 0, 1] * matR;

scaleMatrix = zeros(size(matR));
scaleMatrix([1, 5, 9]) = diag(matR);
shearMat = scaleMatrix \ matR;
%% Select matched pairs
th_inconsistancy = 0.1;

X_cc_idx = fun_bin_data_to_idx_list(X_label);
X_cc_num = numel(X_cc_idx);
X_cc_size = cellfun(@numel, X_cc_idx);
% Define inconsistancy as the number of labels / number of voxels
X_cc_inconsistancy = zeros(X_cc_num, 1);
for iter1 = 1 : X_cc_num
    X_cc_inconsistancy(iter1) = numel(unique(Y_label(X_cc_idx{iter1})))/X_cc_size(iter1);
end
matched_cc_x = X_cc_inconsistancy <= th_inconsistancy;


Y_cc_idx = fun_bin_data_to_idx_list(Y_label);
Y_cc_num = numel(Y_cc_idx);
Y_cc_size = cellfun(@numel, Y_cc_idx);
% Define inconsistancy as the number of labels / number of voxels
Y_cc_inconsistancy = zeros(Y_cc_num, 1);
for iter1 = 1 : Y_cc_num
    Y_cc_inconsistancy(iter1) = numel(unique(Y_label(Y_cc_idx{iter1})))/Y_cc_size(iter1);
end
matched_cc_y = Y_cc_inconsistancy <= th_inconsistancy;
match_x_Q = false(numel(X_label),1);
match_x_Q(cat(2, X_cc_idx{matched_cc_x})) = true;

match_y_Q = false(numel(Y_label),1);
match_y_Q(cat(2, Y_cc_idx{matched_cc_y})) = true;
matched_Q = match_x_Q & match_y_Q;

if any(matched_cc_x)
    X_ = X_(matched_Q,:);
    Y_ = Y_(matched_Q,:);
    tY_ = tY_(matched_Q, :);
    dist_X_tY = sqrt(sum((X_ - tY_) .^ 2, 2));
else
    X_ = [];
    Y_ = [];
end

if size(X_,1) < 3
    rate = 0;
    consistent_rate = 0; % Too less points matched
else
    consistent_rate = size(X_,1) / numel(X_label);
end
%%
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) * 1.5 .* [2,1];
tiledlayout(fig_hdl, 1, 2);
ax_hdl = nexttile;
scatter3(ax_hdl, X_(:, 1), X_(:, 2), X_(:, 3), '.');
hold(ax_hdl, 'on');
scatter3(ax_hdl, tY_(:, 1), tY_(:, 2), tY_(:, 3), 'o');
ax_hdl.DataAspectRatio = [1,1,1];
ax_hdl.XLim(1) = 0;
ax_hdl.YLim(1) = 0;
ax_hdl.ZLim(1) = 0;
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
leg_hdl = legend(ax_hdl, 'Post perfusion', 'Final');
leg_hdl.Location = 'northeast';
ax_hdl_2 = nexttile;
histogram(ax_hdl_2, dist_X_tY, 0 : 0.5 : 10, 'Normalization', 'count');
leg_string = sprintf('Number of pairs: %d\nLinear scale factor %.4f\nVolume scale factor %.4f\nDiagonal element of scale matrix:\n    [%.4f, %.4f, %.4f]\nOff-diagonal element of shear matrix:\n    [%.4f, %.4f, %.4f]\nRotation matrix:\n    [%.2f, %.2f, %.2f;\n     %.2f, %.2f, %.2f;\n     %.2f, %.2f, %.2f]\nMean displacement: %.2f (\\mum)\nMedian displacement %.2f (\\mum)', ...
    size(X_, 1), (det(scaleMatrix))^(1/3), det(scaleMatrix), abs(diag(scaleMatrix)), ...
    shearMat([4, 7, 8]), matQ, mean(dist_X_tY), median(dist_X_tY));
leg_hdl_2 = legend(ax_hdl_2, leg_string, 'Location', 'northeast');
ax_hdl_2.XLabel.String = 'Residual displacement (\mum)';
ax_hdl_2.YLabel.String = 'Counts';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');
ax_hdl.ZDir = 'reverse';

fig_fp = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'ML20190124'), 'Post_perfusion_vs_final', ...
    sprintf('WholeBrain_ML20190124_post_perfusion_vs_final_affine_registration_euiheon_calibration.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);
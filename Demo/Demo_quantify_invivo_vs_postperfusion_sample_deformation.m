clc;clear;close all;
DataManager = FileManager;
vg_fn_fp = './Demo/Data/Demo_vessel_graph.mat';
vg_pf_fp = './Demo/Data/Demo_post_perfusion_graph.mat';

vg_fn = DataManager.load_data(vg_fn_fp);
vg_pf = DataManager.load_data(vg_pf_fp);
%% Compute the radius of the vessel segment
vg_fn.link.features = fun_analysis_get_link_features(vg_fn);
vg_pf.link.features = fun_analysis_get_link_features(vg_pf);
%% Use the penetrating vessel to align first
pf_selected_Q = vg_fn.link.features.dt_median >= 3.5;
pf_selected_ind = cat(1, vg_fn.link.cc_ind{pf_selected_Q});
pf_selected_label = full(vg_fn.link.map_ind_2_label(pf_selected_ind));

iv_selected_Q = vg_pf.link.features.dt_median >= 3.5;
iv_selected_ind = cat(1, vg_pf.link.cc_ind{iv_selected_Q});
iv_selected_label = full(vg_pf.link.map_ind_2_label(iv_selected_ind));

downsample_step = 2;
skel_1_block_size = vg_pf.num.mask_size;
skel_2_block_size = vg_fn.num.mask_size;

skel_1_sub = fun_ind2sub(skel_1_block_size, iv_selected_ind(1:downsample_step:end));
% Correct the voxel size on x-y direction 
skel_1_sub(:, 1:2) = skel_1_sub(:, 1:2) .* 656 / 700;
skel_2_sub = fun_ind2sub(skel_2_block_size, pf_selected_ind(1:downsample_step:end));

skel_1_label = iv_selected_label(1:downsample_step:end);
skel_2_label = pf_selected_label(1:downsample_step:end);
%% Estimate initial transformation
% Estimate the rough rigid transformation parapeters by visualizing vessel
% skeletons in both vessel graphs

% Change the following transformation parameters: 
skel_2_com = mean(skel_2_sub, 1);
skel_2_r = skel_2_sub - skel_2_com;
rot_x_agl = -(20 / 180 * pi);
rot_y_agl = -pi/2;
rot_z_agl = 0;
Rx = [1, 0, 0; 0, cos(rot_x_agl), -sin(rot_x_agl); 0, sin(rot_x_agl), cos(rot_x_agl)];
Ry = [cos(rot_y_agl), 0, sin(rot_y_agl); 0, 1, 0; -sin(rot_y_agl), 0, cos(rot_y_agl)];
Rz = [cos(rot_z_agl), -sin(rot_z_agl), 0; sin(rot_z_agl), cos(rot_z_agl), 0; 0, 0, 1];
% Transform skeleton voxels:
skel_2_new = [1, 0, 0; 0, -1, 0; 0, 0, 1] * Rx * Ry * skel_2_r';
skel_2_new = skel_2_new' + skel_2_com;
skel_pre_regid_T = [1, 0, 0; 0, -1, 0; 0, 0, 1] * Rx * Ry;
% Visualization
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 1];
tiledlayout(fig_hdl, 1, 3);
ax_hdl_1 = nexttile;
scatter3(ax_hdl_1, skel_1_sub(:, 1), skel_1_sub(:, 2), skel_1_sub(:, 3));
ax_hdl_1.XLabel.String = 'X';
ax_hdl_1.YLabel.String = 'Y';
ax_hdl_1.ZLabel.String = 'Z';
ax_hdl_1.DataAspectRatio = [1,1,1];
ax_hdl_1.Title.String = 'Skeleton in graph 1';
ax_hdl_2 = nexttile;
scatter3(ax_hdl_2, skel_2_sub(:, 1), skel_2_sub(:, 2), skel_2_sub(:, 3));
ax_hdl_2.DataAspectRatio = [1,1,1];
ax_hdl_2.XLabel.String = 'X';
ax_hdl_2.YLabel.String = 'Y';
ax_hdl_2.ZLabel.String = 'Z';
ax_hdl_2.Title.String = 'Skeleton in graph 2';
ax_hdl_3 = nexttile;
scatter3(ax_hdl_3, skel_2_new(:, 1), skel_2_new(:, 2), skel_2_new(:, 3));
ax_hdl_3.DataAspectRatio = [1,1,1];
ax_hdl_3.XLabel.String = 'X';
ax_hdl_3.YLabel.String = 'Y';
ax_hdl_3.ZLabel.String = 'Z';
ax_hdl_3.Title.String = 'Skeleton in graph 2 after transformation';
%% Sample skeleton voxels
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

[Transform, C] = cpd_register(des_1_sub, des_2_sub_shift, cpd_opt);
dist_1_to_2t = pdist2(des_1_sub, Transform.Y);
[mutual_nn_idx_1, mutual_nn_idx_2, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_1_to_2t, false);

X_ = des_1_sub(mutual_nn_idx_1,:);
Y_ = des_2_sub_shift(mutual_nn_idx_2,:);
tY_= Transform.Y(mutual_nn_idx_2,:);
%% Decompose Affine matrix to rotation matrix and shear matrix
[matQ, matR] = qr(Transform.R);
scaleMatrix = zeros(size(matR));
scaleMatrix([1, 5, 9]) = diag(matR);
shearMat = scaleMatrix \ matR;
%% Transform the entire perfusion dataset according to the affine registration of the penetrating vessels
iv_ind = vg_pf.link.pos_ind;
iv_label = full(vg_pf.link.map_ind_2_label(iv_ind));

pf_ind = vg_fn.link.pos_ind;
pf_label = full(vg_fn.link.map_ind_2_label(pf_ind));

downsample_step = 2;
skel_1_block_size = vg_pf.num.mask_size;
skel_2_block_size = vg_fn.num.mask_size;

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

% Keep the top 400 um skeleton
max_depth_um = 350;
skel_1_to_keep_Q = skel_1_sub(:, 3) <= max_depth_um;
skel_1_sub = skel_1_sub(skel_1_to_keep_Q, :);
skel_1_label = skel_1_label(skel_1_to_keep_Q);

skel_2_to_keep_Q = skel_2_sub_T(:, 3) <= max_depth_um;
skel_2_sub_T = skel_2_sub_T(skel_2_to_keep_Q, :);
skel_2_label = skel_2_label(skel_2_to_keep_Q);

warning('The following function requires at least 32GB of RAM');
[des_1_sub, des_2_sub_shift, des_1_label, des_2_label] = fun_feature_match_sample_descriptor(skel_1_sub, ...
    skel_2_sub_T, skel_1_label, skel_2_label, skel_sample_opt);
%% Visualize the transformed vessel skeleton
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
ax_hdl_1 = axes(fig_hdl);
scatter3(ax_hdl_1, des_1_sub(:, 1), des_1_sub(:, 2), des_1_sub(:, 3));
ax_hdl_1.XLabel.String = 'X';
ax_hdl_1.YLabel.String = 'Y';
ax_hdl_1.ZLabel.String = 'Z';
ax_hdl_1.DataAspectRatio = [1,1,1];

hold(ax_hdl_1, 'on');
scatter3(ax_hdl_1, des_2_sub_shift(:, 1), des_2_sub_shift(:, 2), des_2_sub_shift(:, 3));
%% Coherent Point Drift on the entire vessel skeleton (including microvessels)
cpd_opt.method = 'affine';
cpd_opt.rot = 1;
cpd_opt.scale = 1;
cpd_opt.viz = 1;              % show every iteration
cpd_opt.outliers = 0.7;       % use 0.7 noise weight
cpd_opt.max_it = 100;        % Maxinum number of iteration
cpd_opt.fgt = 0;              % do not use FGT (default)
cpd_opt.tol = 1e-5;         % Error tolorance ( how is it defined? )
cpd_opt.normalize = 1;        % normalize to unit variance and zero mean before registering (default)
cpd_opt.corresp = 1;          % compute correspondence vector at the end of registration (not being estimated by default)
[Transform_all_skel, C] = cpd_register(des_1_sub, des_2_sub_shift, cpd_opt);
dist_1_to_2t = pdist2(des_1_sub, Transform_all_skel.Y);

[mutual_nn_idx_1, mutual_nn_idx_2, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_1_to_2t, false);

X_ = des_1_sub(mutual_nn_idx_1,:);
Y_ = des_2_sub_shift(mutual_nn_idx_2,:);
tY_= Transform_all_skel.Y(mutual_nn_idx_2,:);
X_label = des_1_label(mutual_nn_idx_1);
Y_label = des_2_label(mutual_nn_idx_2);
%% Decompose the Affine matrix
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
%% Visualization
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
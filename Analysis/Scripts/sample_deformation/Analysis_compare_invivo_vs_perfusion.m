clc;clear;close all;
dataset_name = 'DKLab';
stack_invivo = 'Rui20180622_invivo';
stack_perfusion = 'Rui20180622_perfusion';
grid_c_version = 'merged';
DataManager = FileManager;
vg_iv = DataManager.load_graph_in_block(dataset_name, stack_invivo, grid_c_version, 0, 0, 0);
vg_pf = DataManager.load_graph_in_block(dataset_name, stack_perfusion, grid_c_version, 0, 0, 0);
%% Initialize the relative displacement by masked FFT ??
%% Direct point cloud registration
downsample_step = 2;
invivo_skel_sub = fun_ind2sub(vg_iv.num.mask_size, vg_iv.link.pos_ind(1:downsample_step:end));
perfusion_skel_sub = fun_ind2sub(vg_pf.num.mask_size, vg_pf.link.pos_ind(1:downsample_step:end));

skel_sample_opt = struct;
skel_sample_opt.total_num_descriptor = 1e4;
skel_sample_opt.rm_descriptor_by_pdist2_Q = true;
skel_sample_opt.max_disp_pixel_yxz = [20, 20, 20];

[des_1_sub, des_2_sub_shift, des_1_label, des_2_label] = fun_feature_match_sample_descriptor(invivo_skel_sub, ...
    perfusion_skel_sub, vg_iv.link.label(1:downsample_step:end), vg_pf.link.label(1:downsample_step:end), skel_sample_opt);

des_1_ind = sub2ind(vg_iv.num.mask_size, des_1_sub(:,1), des_1_sub(:,2), des_1_sub(:,3));
X_ind_2_label = sparse(des_1_ind, ones(length(des_1_ind),1), des_1_label, prod(vg_iv.num.mask_size),1);

des_2_ind = sub2ind(vg_pf.num.mask_size, des_2_sub_shift(:,1), des_2_sub_shift(:,2), des_2_sub_shift(:,3));
Y_ind_2_label = sparse(des_2_ind, ones(length(des_2_ind),1), des_2_label, prod(vg_pf.num.mask_size),1);
%% Coherent point drift
cpd_opt.method = 'rigid';
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
[Transform, C] = cpd_register(des_1_sub, des_2_sub_shift, cpd_opt);

dist_1_to_2t = pdist2(des_1_sub, Transform.Y);
[mutual_nn_idx_1, mutual_nn_idx_2, mutual_nn_dist] = fun_find_col_row_co_minimum(dist_1_to_2t, false);

X_ = des_1_sub(mutual_nn_idx_1,:);
Y_ = des_2_sub_shift(mutual_nn_idx_2,:);
tY_= Transform.Y(mutual_nn_idx_2,:);
%% Select matched pairs
th_inconsistancy = 0.1;

X_ind = des_1_ind(mutual_nn_idx_1);
Y_ind = des_2_ind(mutual_nn_idx_2);
X_label = full(X_ind_2_label(X_ind));
Y_label = full(Y_ind_2_label(Y_ind));

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
match_x_Q = false(numel(X_ind),1);
match_x_Q(cat(2, X_cc_idx{matched_cc_x})) = true;

match_y_Q = false(numel(X_ind),1);
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
    consistent_rate = size(X_,1) / numel(X_ind);
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
ax_hdl.ZDir = 'reverse';
ax_hdl.XLabel.String = 'X (\mum)';
ax_hdl.YLabel.String = 'Y (\mum)';
ax_hdl.ZLabel.String = 'Z (\mum)';
ax_hdl.FontSize = 14;
leg_hdl = legend(ax_hdl, 'In vivo', 'Post perfusion');
leg_hdl.Location = 'northeast';
ax_hdl_2 = nexttile;
histogram(ax_hdl_2, dist_X_tY, 0 : 0.5 : 10, 'Normalization', 'count');
leg_string = sprintf('Number of pairs: %d\nRigid scale factor %.4f\nMean displacement: %.2f (\\mum)\nMedian displacement %.2f (\\mum)', ...
    size(X_, 1), Transform.s, mean(dist_X_tY), median(dist_X_tY));
leg_hdl_2 = legend(ax_hdl_2, leg_string, 'Location', 'northeast');
ax_hdl_2.XLabel.String = 'Residual displacement (\mum)';
ax_hdl_2.YLabel.String = 'Counts';
ax_hdl_2.FontSize = 14;
box(ax_hdl_2, 'off');

fig_fp = fullfile(DataManager.fp_visualization_folder('WholeBrain', 'ML20190124'), 'Invivo_vs_postPerfusion', ...
    sprintf('WholeBrain_ML20190124_invivo_postperfusion_skel_displacement_after_regid_T.png'));
fun_print_image_in_several_formats(fig_hdl, fig_fp);



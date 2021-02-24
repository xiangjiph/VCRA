function [X_stable,Y_stable,rate, pixshift] = fun_searchpair_vessel_edges(descriptor_1,descriptor_2,pixshiftinit)
% fun_searchpair_vessel_edges apply Coherent Point Drift for two point
% cloud (surface) registration. Since surface point cloud is dense, the
% point clouds are downsampled by mergering nearby points and reolaced by
% their mean position( fun_stitching_merge)surface_voxels ). The resulting
% point clouds are fed into CPD. The output voxel pairs are further
% selected by compute the local displacement standard deviation. Only local
% blocks with standard deviation smaller than a threshold are kept. 
% Input: 
%   descriptor_1_sub: N-by-3 numerical array, position of the points
%   descriptor_2_sub: M-by-3 numerical array, position of the points in the
%   other point cloud. 
%   pixshiftinit: 3-by-1 numerical array, initial estimation of the
%   translation displacment between two point clouds
% Output: 
%   X_stable: N-by-3 numerical array, position of the matched points
%   Y_stable: N-by-3 numerical array, position of the matched points
%   rate: numerical scalar, defined by Erhan, hard to explain... see
%   vessel_descriptorMatchforz. A value that reveals the reliable of
%   matching ( the higher, the better ) 
%   pixshift: 3-by-1 numerical array, median of the displacement of the
%   output paired point sets. 
% 
% Author: Xiang Ji, UC San Diego ( xiangji.ucsd@gmail.com )
% Date: Dec 5, 2018

% If the number of available feature points is larger than this number,
% downsample the point cloud.
% th_num_ds_desc = 15000; 
% max_disp_um = 15;
X_stable = [];
Y_stable = [];
rate = 0;
max_edge_voxel_point = 2e4;
max_edge_voxel_for_neighbor_search = 3 * max_edge_voxel_point;
% pixshift = [20, 0, 143];
pixshift =  pixshiftinit;
[~, iadj] = max(pixshiftinit);
switch iadj
    case 1
        max_disp_pixel_yxz = [15, 10, 5];
    case 2
        max_disp_pixel_yxz = [10, 15, 5];
    case 3
        max_disp_pixel_yxz = [30, 30, 20];
end
%% Transform the descriptor according to the estimated shift first
if isempty(descriptor_1) || isempty(descriptor_2)
    return;
end
descriptor_1_sub = descriptor_1(:,1:3);
descriptor_2_sub = descriptor_2(:,1:3);
desc_1_gradient = descriptor_1(:, 4);
desc_2_gradient = descriptor_2(:, 4);
desc_2_sub_shifted = bsxfun(@plus, descriptor_2_sub, pixshift);
desc_2_sub_shifted_max = max(desc_2_sub_shifted, [], 1);
desc_2_sub_shifted_min = min(desc_2_sub_shifted, [], 1);

desc_1_sub_max = max(descriptor_1_sub, [], 1);
desc_1_sub_min = min(descriptor_1_sub, [], 1);
% Add offset to remove the boundary produced by cutting and region without
% signal
overlap_bbox_max = min(desc_2_sub_shifted_max, desc_1_sub_max) - 10;
overlap_bbox_min = max(desc_2_sub_shifted_min, desc_1_sub_min) + 10;
desc_1_selected_Q = all(bsxfun(@ge, descriptor_1_sub, overlap_bbox_min) & bsxfun(@le, descriptor_1_sub, overlap_bbox_max), 2);
desc_2_selected_Q = all(bsxfun(@ge, desc_2_sub_shifted, overlap_bbox_min) & bsxfun(@le, desc_2_sub_shifted, overlap_bbox_max), 2);

desc_1_ds = cat(2, descriptor_1_sub(desc_1_selected_Q,:), desc_1_gradient(desc_1_selected_Q));
desc_2_ds = cat(2, desc_2_sub_shifted(desc_2_selected_Q,:), desc_2_gradient(desc_2_selected_Q));

if isempty(desc_1_ds) || isempty(desc_2_ds)
    return;
end
%% Merge edge voxels
merge_box_size = [6,6,2];
% Threshold for local displacement standard deviation
% Correspond to 0.5 um
std_th_1 = 1.5;
std_th_2 = 1.5;
std_th_3 = 0.5;
desc_1_ds = fun_stitching_merge_surface_voxels(desc_1_ds, merge_box_size);
desc_2_ds = fun_stitching_merge_surface_voxels(desc_2_ds, merge_box_size);
%% Descriptor pair selection: 
% The distance between two point set (after shifted by initial estimation)
% should not be too large
% Limit the number of voxels before pdist2 to 3 * max_edge_voxel_point

if isempty(desc_1_ds) || isempty(desc_2_ds)
    return;
else
    if size(desc_1_ds, 1) > max_edge_voxel_for_neighbor_search
        % Further select voxels with large gradients
        [~, desc_1_idx] = sort(desc_1_ds(:, 4), 'descend');
        desc_1_ds = desc_1_ds(desc_1_idx(1:max_edge_voxel_for_neighbor_search), :);
    end
    if size(desc_2_ds, 1) > max_edge_voxel_for_neighbor_search
        [~, desc_2_idx] = sort(desc_2_ds(:, 4), 'descend');
        desc_2_ds = desc_2_ds(desc_2_idx(1:max_edge_voxel_for_neighbor_search), :);
    end
end  
tmp_pdist = pdist2(single(desc_1_ds(:,1)), single(desc_2_ds(:,1)));
%     tmp_pdist3 = (tmp_pdist.^2) ./3;
tmp_pdist_reasonable = tmp_pdist < max_disp_pixel_yxz(1);
tmp_pdist = pdist2(single(desc_1_ds(:,2)), single(desc_2_ds(:,2)));
%     tmp_pdist3 = tmp_pdist3 + (tmp_pdist.^2) ./3;
tmp_pdist_reasonable = tmp_pdist_reasonable & tmp_pdist < max_disp_pixel_yxz(2);
tmp_pdist = pdist2(single(desc_1_ds(:,3)), single(desc_2_ds(:,3)));
%     tmp_pdist3 = sqrt(tmp_pdist3 + (tmp_pdist.^2));
tmp_pdist_reasonable = tmp_pdist_reasonable & tmp_pdist < max_disp_pixel_yxz(3);
clearvars tmp_pdist
desc_1_close_neighbor_Q = any(tmp_pdist_reasonable, 2);
desc_2_close_neighbor_Q = any(tmp_pdist_reasonable, 1)';
clearvars tmp_pdist_reasonable
desc_1_ds = desc_1_ds(desc_1_close_neighbor_Q, :);
desc_2_ds = desc_2_ds(desc_2_close_neighbor_Q, :);

if isempty(desc_1_ds) || isempty(desc_2_ds)
    return;
else
    if size(desc_1_ds, 1) > max_edge_voxel_point 
        % Further select voxels with large gradients
        [~, desc_1_idx] = sort(desc_1_ds(:, 4), 'descend');
        desc_1_ds = desc_1_ds(desc_1_idx(1:max_edge_voxel_point), 1:3);
    else
        desc_1_ds = desc_1_ds(:, 1:3);
    end
    if size(desc_2_ds, 1) > max_edge_voxel_point
        [~, desc_2_idx] = sort(desc_2_ds(:, 4), 'descend');
        desc_2_ds = desc_2_ds(desc_2_idx(1:max_edge_voxel_point), 1:3);
    else
        desc_2_ds = desc_2_ds(:, 1:3);
    end
end  
if size(desc_1_ds, 1) < 3 || size(desc_2_ds, 1) < 3
    return;
end
%% Point cloud registration 
% The purpose of point cloud registration is to find the correspondance.
% Will it be better to use affine transformation for the x-y matching? 
projectionThr = 5;
optimopts = statset('nlinfit');
optimopts.RobustWgtFun = 'bisquare';
opt.method='nonrigid_lowrank';
opt.beta=6;            % the width of Gaussian kernel (smoothness)
opt.lambda=16;          % regularization weight
opt.viz=false;              % show every iteration
opt.outliers=0.9;       % use 0.7 noise weight
min_num_voxel = min(size(desc_1_ds, 1), size(desc_2_ds, 1));
opt.numeig = min(100, min_num_voxel);       % Number of eigenvectors for low rank appriximation 
opt.eigfgt = 0;         % Use fast gaussian transformation for computing eigenvectors
opt.max_it = 50;        % Maxinum number of iteration 
opt.fgt=0;              % do not use FGT (default)
opt.tol = 1e-3;         % Error tolorance ( how is it defined? )
opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
matchparams.optimopts = optimopts;
matchparams.opt = opt;
matchparams.projectionThr = projectionThr;
matchparams.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
matchparams.debug = false;
matchparams.viz = false;
[rate, X_, Y_, tY_] = vessel_descriptorMatchforz(desc_1_ds, desc_2_ds, pixshift, matchparams);
%% Affine transformation 
% projectionThr = 5;
% opt.method='affine';
% opt.beta=6;            % the width of Gaussian kernel (smoothness)
% opt.lambda=16;          % regularization weight
% opt.viz=1;              % show every iteration
% opt.outliers=0.9;       % use 0.7 noise weight
% opt.max_it = 50;        % Maxinum number of iteration 
% opt.fgt=0;              % do not use FGT (default)
% opt.tol = 1e-3;         % Error tolorance ( how is it defined? )
% opt.normalize=1;        % normalize to unit variance and zero mean before registering (default)
% opt.corresp=1;          % compute correspondence vector at the end of registration (not being estimated by default)
% matchparams.optimopts = optimopts;
% matchparams.opt = opt;
% matchparams.projectionThr = projectionThr;
% matchparams.model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
% matchparams.debug = false;
% matchparams.viz = true;
% tic
% [rate_af, X_af, Y_af, tY_af] = vessel_descriptorMatchforz(desc_1_ds, desc_2_ds, pixshift, matchparams);
% toc
%% Matched point selection 
%  If the matching is bad, return empay matching directly
if rate < 0.8 || isempty(X_) || isempty(Y_)
    X_stable = [];
    Y_stable = [];
    rate = 0;
    return;
end
% Check the consistancy of the local transformation 
% Compute the standard deviation of the displacement of the matched pairs
% in small blocks ( say, of size 10x10x10 um ) 
% displacement between X_ and Y_ for each local block 
disp_X_Y = X_ - Y_;
sub_min = min(X_, [], 1);
sub_max = max(X_, [], 1);
sub_bbox = sub_max - sub_min + 1;
local_bbox_size = [30,30,10];
grid_size = ceil(sub_bbox ./ local_bbox_size);
grid_sub = ceil((1 + bsxfun(@minus, X_, sub_min))./ local_bbox_size);
grid_ind = sub2ind(grid_size, grid_sub(:,1), grid_sub(:,2), grid_sub(:,3));
[grid_idx_cell, grid_idx_list]= fun_bin_data_to_idx_list(grid_ind);
num_valid_bbox = numel(grid_idx_list);
grid_idx_cell_array = cell(grid_size);
num_pair = zeros(grid_size);
% mean_dis_1 = zeros(grid_size);
% mean_dis_2 = zeros(grid_size);
% mean_dis_3 = zeros(grid_size);
std_dis_1 = zeros(grid_size) - 3;
std_dis_2 = zeros(grid_size) - 3;
std_dis_3 = zeros(grid_size) - 3;
for iter_bbox = 1 : num_valid_bbox
    dis_list = disp_X_Y(grid_idx_cell{iter_bbox},:);
    tmp_grid_ind = grid_idx_list(iter_bbox);
    tmp_num_pair = size(dis_list, 1);
    num_pair(tmp_grid_ind) = tmp_num_pair;
    grid_idx_cell_array{tmp_grid_ind} = grid_idx_cell{iter_bbox};
    if tmp_num_pair > 1
%         mean_dis_1(tmp_grid_ind) = mean(dis_list(:,1));
%         mean_dis_2(tmp_grid_ind) = mean(dis_list(:,2));
%         mean_dis_3(tmp_grid_ind) = mean(dis_list(:,3));
        std_dis_1(tmp_grid_ind) = std(dis_list(:,1));
        std_dis_2(tmp_grid_ind) = std(dis_list(:,2));
        std_dis_3(tmp_grid_ind) = std(dis_list(:,3));
    elseif tmp_num_pair == 1
%         mean_dis_1(tmp_grid_ind) = dis_list(1);
%         mean_dis_2(tmp_grid_ind) = dis_list(:,2);
%         mean_dis_3(tmp_grid_ind) = dis_list(:,3);
    end
end
% Control local displacement variance
grid_low_std_Q = (std_dis_1 <= std_th_1) & (std_dis_2 <= std_th_2) & (std_dis_3 <= std_th_3) & ...
    ( num_pair > 1 );
grid_low_std_voxel_idx = cat(2, grid_idx_cell_array{grid_low_std_Q});
X_stable = X_(grid_low_std_voxel_idx,:);
Y_stable = Y_(grid_low_std_voxel_idx,:);
% The initial estimation of the stage position is quite close. If the
% displacement between Y_ and tY_ is too large, it must be a wrong
% matching
% Control global displacement variance
disp_X_Y = X_stable - Y_stable;
disp_X_Y_med = median(disp_X_Y);
disp_X_Y_dev = disp_X_Y - disp_X_Y_med;
disp_X_Y_tol = min(15, max(5,std(single(disp_X_Y_dev),1) * 3));
disp_inlier = all(abs(disp_X_Y_dev) < disp_X_Y_tol, 2);
X_stable = round(X_stable(disp_inlier, :));
Y_stable = round(Y_stable(disp_inlier, :));
pixshift = round(median(X_stable - Y_stable));
%% Visualize matched points
% figure;
% subplot(1,3,1)
% scatter3(desc_1_selected(:, 1), desc_1_selected(:, 2), desc_1_selected(:, 3));
% hold on
% scatter3(desc_2_selected(:, 1), desc_2_selected(:, 2), desc_2_selected(:, 3));
% xlabel 'x';
% ylabel 'y';
% zlabel 'z';
% legend('Tile 1', 'Tile 2');
% title('Descriptors in the overlapping region');
% subplot(1,3,2)
% scatter3(X_(:, 1), X_(:, 2), X_(:, 3));
% hold on
% scatter3(Y_(:, 1) + pixshiftinit(1), Y_(:, 2) + pixshiftinit(2), Y_(:, 3) + pixshiftinit(3));
% xlabel 'x';
% ylabel 'y';
% zlabel 'z';
% legend('Tile 1', 'Tile 2');
% title('Matched Descriptors in the overlapping region');
% subplot(1,3,3)
% scatter3(X_stable(:,1), X_stable(:,2), X_stable(:,3))
% hold on 
% scatter3(Y_stable(:,1) + pixshiftinit(1), Y_stable(:,2) + pixshiftinit(2), Y_stable(:,3) + pixshiftinit(3))
% legend('Tile 1', 'Tile 2');
% xlabel('x');
% ylabel('y');
% zlabel('z');
% title('Matched edge voxel with low local variance in displacement')
%% Visualization of the displacement field
% vis_sec = 1;
% figure;
% subplot(2,3,1)
% imagesc(mean_dis_1(:,:,vis_sec));
% title('Mean displacement in X')
% colorbar
% subplot(2,3,2)
% imagesc(mean_dis_2(:,:,vis_sec));
% title('Mean displacement in Y')
% colorbar
% subplot(2,3,3)
% imagesc(mean_dis_3(:,:,vis_sec));
% title('Mean displacement in Z')
% colorbar
% subplot(2,3,4)
% imagesc(std_dis_1(:,:,vis_sec));
% title('STD of displacement in X');
% colorbar
% subplot(2,3,5)
% imagesc(std_dis_2(:,:,vis_sec));
% title('STD of displacement in Y')
% colorbar
% subplot(2,3,6)
% imagesc(std_dis_3(:,:,vis_sec));
% title('STD of displacement in Z')
% colorbar
% figure;
% scatter3(X_stable(:,1), X_stable(:,2), X_stable(:,3))
% hold on 
% scatter3(Y_stable(:,1), Y_stable(:,2), Y_stable(:,3))
% legend('Tile 1', 'Tile 2');
% xlabel('X');
% ylabel('Y');
% zlabel('Z');
% title('Matched edge voxel with low local variance in displacement')
% view(2)
%% 
end

%% Sub function
% function voxel_sub = fun_stitching_merge_surface_voxels(voxel_sub, merge_block_size)
% % fun_stitching_merge_surface_voxels merges the input voxel subscript list
% % by computing the average position of the voxels within the block, whose
% % size is specified by merge_block_size
% % Input: 
% %   voxel_sub: N-by-3, coordinate of the voxel position in 3D space
% %   merge_block_size: 1-by-3 numerical vector, size of the block for merging. 
% % Output: 
% %   voxel_sub: N'-by-3 numerical array after merging
% sub_min = min(voxel_sub, [], 1);
% sub_max = max(voxel_sub, [], 1);
% image_size = sub_max - sub_min + 1;
% downsampled_image_size = ceil(image_size ./ merge_block_size);
% 
% num_edge_voxel = zeros(downsampled_image_size);
% mean_edge_pos_1 = zeros(downsampled_image_size);
% mean_edge_pos_2 = zeros(downsampled_image_size);
% mean_edge_pos_3 = zeros(downsampled_image_size);
% bbox_sub = ceil((1 + bsxfun(@minus, voxel_sub, sub_min))./ merge_block_size);
% for iter1 = 1 : size(bbox_sub, 1)
%     num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
%         num_edge_voxel(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + 1;
%     mean_edge_pos_1(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
%         mean_edge_pos_1(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 1);
%     mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
%         mean_edge_pos_2(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 2);
%     mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) = ...
%         mean_edge_pos_3(bbox_sub(iter1,1), bbox_sub(iter1,2), bbox_sub(iter1,3)) + voxel_sub(iter1, 3);
% end
% voxel_sub = cat(2, mean_edge_pos_1(:) ./ max(1, num_edge_voxel(:)), ...
%     mean_edge_pos_2(:) ./ max(1, num_edge_voxel(:)), mean_edge_pos_3(:) ./ max(1, num_edge_voxel(:)));
% 
% voxel_sub = voxel_sub(all(voxel_sub > 0, 2),:);
% end
%% Sub function
function [rate,X_,Y_,tY_] = vessel_descriptorMatchforz(X,Y,pixshift,params)
%DESCRIPTORMATCH Summary of this function goes here
%
% [OUTPUTARGS] = DESCRIPTORMATCH(INPUTARGS) Explain usage here
%
% Inputs:
%   X, Y: two 2D real, double marices, specifying the position of the point
%   in the point cloud.
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/23 14:09:29 $	$Revision: 0.1 $
% Copyright: HHMI 2016
opt = params.opt;
projectionThr = params.projectionThr;
%% Initial match based on point drift
[Transform, ~] = cpd_register(X,Y,opt);
%% check if match is found
% Compute the pairwise euclidean distance between two input array
pD = pdist2(X,Transform.Y);
[aa1,bb1] = min(pD,[],1);
[~,bb2] = min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)' < projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);
tY_= Transform.Y(keeptheseY(disttrim),:);
% Rate is the ratio of the number of matched pair of distance less than
% projectionThr over the total number of matched pairs
rate = sum(disttrim)/length(disttrim);
% [pixshift rate]
if rate < .5 % dont need to continue
    [X_,Y_] = deal([]);
    return
end
Y_ = bsxfun(@minus, Y_, pixshift);
end
%% Sub function 
function [bin_cell_array, varargout] = fun_bin_data_to_idx_list(data)
% fun_bin_data_to_idx_list bin the data according to their values and
% output the corresponding index list
% Input: 
%   data: numerical vector
% Output: 
%   bin_cell_array: cell array, each cell constains a vector, whose
%   components are the indices of the component of data that have the same
%   value. 
%   varargout: unique data value 
num_data = numel(data);
if ~issorted(data)
    [data, idx_list ]= sort(data, 'ascend');
else
    idx_list = 1 : num_data;
end

bin_size = 0;
bin_data = data(1);
bin_idx = zeros(1, round(num_data/2));
est_num_bin = 500;
bin_value_list = zeros(est_num_bin,1);
bin_value_list(1) = data(1);
bin_cell_array = cell(est_num_bin,1);
num_bin = 0;
for idx = 1 : num_data
    tmp_data = data(idx);
    
    if tmp_data == bin_data
        bin_size = bin_size + 1;
        bin_idx(bin_size) = idx_list(idx);
    else
        num_bin = num_bin + 1;
        bin_cell_array{num_bin} = bin_idx(1 : bin_size);
        bin_data = tmp_data;
        bin_value_list(num_bin + 1) = bin_data;
        bin_idx(1) = idx_list(idx);
        bin_size = 1;
    end
end
num_bin = num_bin + 1;
bin_cell_array{num_bin} = bin_idx(1 : bin_size);
bin_cell_array(num_bin + 1 : end) = [];
bin_value_list = bin_value_list(1 : num_bin);
if nargout > 1
    varargout{1} = bin_value_list;
end
end
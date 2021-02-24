function grid_cube = fun_get_grid_from_mask(mask_3d, block_size, block_overlap_l, dataset_size, valid_grid_min_volume_ratio, parforQ)
% fun_get_grid_from_mask generate 3D grid with overlaps in each dimension.
% The grid will be used for dividing the entire data set into 3D image
% stacks for segmentation.
% Input:
%   mask_3d: 3D logical array, mask of the region of interest
%   block_size: numerical scalar, length of the image stack (cube).
%   block_overlap_l: numerical scalar, number of overlapping voxels in each
%   dimension ( with the 6 neighbors of the central grid)
%   dataset_size: size of the entire grid. Assume the x-y aspectial ratio
%   in mask_3d is the same as in the dataset.
%   valid_grid_min_volume_ratio: numerical scalar between 0 and 1. Minimum
%   volume ratio of the region of interest interesct with the grid that
%   makes the grid valid. Default value is 0
%
% Output:
%   grid_cube: structure with fields define below.
%
% Implemented by Xiang Ji on 04/17/2019

% Modified by Xiang Ji on 04/18/2019
% 1. The name of field 'bbox_xy_linear_idx_mat' doesn't make sense. Change
% it to 'bbox_xy_label_mat', which is the label of valid cube in layer,
% from 1 to the number of valid cube in layer.
% 2. Add more fields

if nargin < 5
    valid_grid_min_volume_ratio = 0;
    parforQ = true;
elseif nargin < 6
    parforQ = true;
end
%%
grid_cube = struct;
[target_voxel_size, dataset_name, stack, data_type, grid_name, mask_version] = deal([]);
% Information about the data set
grid_cube.dataset_name = dataset_name;
grid_cube.stack = stack;
grid_cube.mask_version = mask_version;
grid_cube.version = grid_name;
grid_cube.data_type = data_type;
grid_cube.voxel_size_um = target_voxel_size;
%% Determine the blocks that cover the mask
disp('Determining the blocks that cover the mask');
mask_size = size(mask_3d);
grid2D = fun_generate_grid(block_size, block_overlap_l, dataset_size(1:2));
tic
clear grid_xy_info
grid_xy_info(mask_size(3)) = struct;
need_rescale_mask_Q = ~all(mask_size(:) == dataset_size(:));
if parforQ
    parfor test_section = 1 : mask_size(3)
        if need_rescale_mask_Q
            brain_mask_xy = imresize(mask_3d(:,:, test_section), dataset_size(1:2));
        else
            brain_mask_xy = mask_3d(:, :, test_section);
        end
        % For each bouding box, check if the brain is located inside the bounding box
        n_bbox = size(grid2D.mmxx, 1);
        list_bbox_contains_sample = ones(1,n_bbox, 'logical');
        mat_bbox_contains_sample = zeros(grid2D.size, 'logical');
        mat_bbox_RIO_ratio = zeros(grid2D.size);
        for test_bbox_id = 1 : n_bbox
            bbox_data = crop_bbox2(brain_mask_xy, grid2D.mmll(test_bbox_id,:), 'default');
            mat_bbox_RIO_ratio(test_bbox_id) = nnz(bbox_data(:))/ prod(grid2D.ll(test_bbox_id, :));
            list_bbox_contains_sample(test_bbox_id) = mat_bbox_RIO_ratio(test_bbox_id) > 0;
            mat_bbox_contains_sample(test_bbox_id) = list_bbox_contains_sample(test_bbox_id);
        end
        list_valid_bbox_idx = find(list_bbox_contains_sample);
        grid_xy_info(test_section).section = test_section;
        grid_xy_info(test_section).bbox_idx= list_valid_bbox_idx;
        grid_xy_info(test_section).bbox_mmll = grid2D.mmll(list_valid_bbox_idx, :);
        grid_xy_info(test_section).bbox_mmxx = grid2D.mmxx(list_valid_bbox_idx, :);
        grid_xy_info(test_section).bbox_mm_idx_pos = grid2D.mm_idx_pos(list_valid_bbox_idx,:);
        grid_xy_info(test_section).bbox_idx_mat = mat_bbox_contains_sample;
        grid_xy_info(test_section).bbox_area_ratio_mat = mat_bbox_RIO_ratio;
        grid_xy_info(test_section).bbox_area_ratio = mat_bbox_RIO_ratio(list_valid_bbox_idx);
    end
    pool_obj = gcp('nocreate');
    delete(pool_obj);
else
    for test_section = 1 : mask_size(3)
        if need_rescale_mask_Q
            brain_mask_xy = imresize(mask_3d(:,:, test_section), dataset_size(1:2));
        else
            brain_mask_xy = mask_3d(:, :, test_section);
        end
        % For each bouding box, check if the brain is located inside the bounding box
        n_bbox = size(grid2D.mmxx, 1);
        list_bbox_contains_sample = ones(1,n_bbox, 'logical');
        mat_bbox_contains_sample = zeros(grid2D.size, 'logical');
        mat_bbox_RIO_ratio = zeros(grid2D.size);
        for test_bbox_id = 1 : n_bbox
            bbox_data = crop_bbox2(brain_mask_xy, grid2D.mmll(test_bbox_id,:), 'default');
            mat_bbox_RIO_ratio(test_bbox_id) = nnz(bbox_data(:))/ prod(grid2D.ll(test_bbox_id, :));
            list_bbox_contains_sample(test_bbox_id) = mat_bbox_RIO_ratio(test_bbox_id) > 0;
            mat_bbox_contains_sample(test_bbox_id) = list_bbox_contains_sample(test_bbox_id);
        end
        list_valid_bbox_idx = find(list_bbox_contains_sample);
        grid_xy_info(test_section).section = test_section;
        grid_xy_info(test_section).bbox_idx= list_valid_bbox_idx;
        grid_xy_info(test_section).bbox_mmll = grid2D.mmll(list_valid_bbox_idx, :);
        grid_xy_info(test_section).bbox_mmxx = grid2D.mmxx(list_valid_bbox_idx, :);
        grid_xy_info(test_section).bbox_mm_idx_pos = grid2D.mm_idx_pos(list_valid_bbox_idx,:);
        grid_xy_info(test_section).bbox_idx_mat = mat_bbox_contains_sample;
        grid_xy_info(test_section).bbox_area_ratio_mat = mat_bbox_RIO_ratio;
        grid_xy_info(test_section).bbox_area_ratio = mat_bbox_RIO_ratio(list_valid_bbox_idx);
    end
end
toc
%% Construct grid in 3D
disp('Constructing the 3D grid');
tic
gridZ = fun_generate_grid(block_size, block_overlap_l, dataset_size(3));
grid_size_3d = [grid2D.size, gridZ.size];
num_layer = gridZ.size;
% Grid Size
grid_cube.grid_size = grid_size_3d;
grid_cube.num_grid_layer = num_layer;
grid_cube.data_size = dataset_size;
grid_cube.data_xy_size = dataset_size(1:2);
grid_cube.data_sec_num = dataset_size(3);
% Cube size
grid_cube.block_xy_size = block_size;
grid_cube.block_z_size = block_size;
grid_cube.block_size = ones(1,3) * block_size;
grid_cube.block_overlap = block_overlap_l;
% Valid cube bounding box
grid_cube.num_bbox_xy = zeros(num_layer, 1);
grid_cube.layer = (1:num_layer).';
grid_cube.bbox_xy_pos_ind_in_layer = cell(num_layer, 1);
grid_cube.bbox_xy_pos_sub_in_layer = cell(num_layer, 1);

grid_cube.bbox_xy_mmll = cell(num_layer, 1);
grid_cube.bbox_xy_mmxx = cell(num_layer, 1);

grid_cube.bbox_grid_sub = cell(num_layer, 1);
grid_cube.bbox_grid_ind = cell(num_layer, 1);

grid_cube.bbox_xy_valid_mat = cell(num_layer, 1);
grid_cube.bbox_xy_label_mat = cell(num_layer, 1);

grid_cube.bbox_z_mmll = cell(num_layer, 1);
grid_cube.bbox_z_mmxx = cell(num_layer, 1);
grid_cube.bbox_xyz_mmll = cell(num_layer, 1);
grid_cube.bbox_xyz_mmxx = cell(num_layer, 1);
grid_cube.bbox_xyz_mmll_in_layer = cell(num_layer, 1);
grid_cube.bbox_xyz_mmxx_in_layer = cell(num_layer, 1);

grid_cube.bbox_volume_ratio = cell(num_layer, 1);
grid_cube.bbox_volume_ratio_array = zeros(grid_size_3d);

mask_z_downsample_rate = dataset_size(3) / mask_size(3);

for layer_idx = 1 : num_layer
    sec_1 = max(floor(( gridZ.mmxx(layer_idx,1) - 1 ) / mask_z_downsample_rate ) + 1, 1);
    sec_2 = min(floor(( gridZ.mmxx(layer_idx,2) - 1 ) / mask_z_downsample_rate ) + 1, mask_size(3));
    
    cube_volume_ratio = mean(cat(3, grid_xy_info(sec_1: sec_2).bbox_area_ratio_mat),3);
    cube_contains_sampleQ = cube_volume_ratio > valid_grid_min_volume_ratio;
    
    list_valid_cube_idx = find(cube_contains_sampleQ);
    num_valid_cube_in_layer = nnz(cube_contains_sampleQ);
    % Linear index of the location of the valid cubes in 2D layer
    grid_cube.bbox_xy_pos_ind_in_layer{layer_idx} = list_valid_cube_idx;
    % Subscript of the location of the valid cubes in 2D layer
    grid_cube.bbox_xy_pos_sub_in_layer{layer_idx} = grid2D.mm_idx_pos(list_valid_cube_idx,:);
    % bounding box parameter in layer
    grid_cube.bbox_xy_mmll{layer_idx} = grid2D.mmll(list_valid_cube_idx, :);
    grid_cube.bbox_xy_mmxx{layer_idx} = grid2D.mmxx(list_valid_cube_idx, :);
    % Grid position linear indices and subscripts for valid cubes
    tmp_sub = cat(2, grid2D.mm_idx_pos(list_valid_cube_idx,:), layer_idx .* ones(num_valid_cube_in_layer,1));
    grid_cube.bbox_grid_sub{layer_idx} = tmp_sub;
    grid_cube.bbox_grid_ind{layer_idx} = sub2ind(grid_size_3d, tmp_sub(:, 1), tmp_sub(:, 2), tmp_sub(:, 3));
    % 2D information in grid matrix form
    grid_cube.bbox_xy_valid_mat{layer_idx} = cube_contains_sampleQ;
    grid_cube.bbox_volume_ratio_array(:, :, layer_idx) = cube_volume_ratio;
    grid_cube.bbox_volume_ratio{layer_idx} = cube_volume_ratio(list_valid_cube_idx);
    % Label of valid cube in layer
    grid_cube.bbox_xy_label_mat{layer_idx} = zeros(size(cube_contains_sampleQ));
    grid_cube.bbox_xy_label_mat{layer_idx}(cube_contains_sampleQ) = 1 : num_valid_cube_in_layer;
    % Number of block in plane
    grid_cube.num_bbox_xy(layer_idx) =  length(list_valid_cube_idx);
    grid_cube.bbox_z_mmll{layer_idx} = gridZ.mmll(layer_idx, :);
    grid_cube.bbox_z_mmxx{layer_idx} = gridZ.mmxx(layer_idx, :);
    
    grid_cube.bbox_xyz_mmll{layer_idx} = cat(2,grid_cube.bbox_xy_mmll{layer_idx}(:,1:2),ones(grid_cube.num_bbox_xy(layer_idx),1)*gridZ.ul(layer_idx),...
        grid_cube.bbox_xy_mmll{layer_idx}(:,3:4),ones(grid_cube.num_bbox_xy(layer_idx) ,1)*gridZ.mmll(layer_idx,2) );
    
    grid_cube.bbox_xyz_mmxx{layer_idx} = cat(2,grid_cube.bbox_xy_mmll{layer_idx}(:,1:2),ones(grid_cube.num_bbox_xy(layer_idx),1)*gridZ.ul(layer_idx),...
        grid_cube.bbox_xy_mmxx{layer_idx}(:,3:4),ones(grid_cube.num_bbox_xy(layer_idx) ,1)*gridZ.mmxx(layer_idx,2) );
    
    grid_cube.bbox_xyz_mmll_in_layer{layer_idx} = cat(2,grid_cube.bbox_xy_mmll{layer_idx}(:,1:2),ones(grid_cube.num_bbox_xy(layer_idx),1),...
        grid_cube.bbox_xy_mmll{layer_idx}(:,3:4),ones(grid_cube.num_bbox_xy(layer_idx) ,1)*gridZ.ll(layer_idx) );
    
    grid_cube.bbox_xyz_mmxx_in_layer{layer_idx} = cat(2,grid_cube.bbox_xy_mmxx{layer_idx}(:,1:2),ones(grid_cube.num_bbox_xy(layer_idx),1),...
        grid_cube.bbox_xy_mmxx{layer_idx}(:,3:4), ones(grid_cube.num_bbox_xy(layer_idx) ,1)*gridZ.ll(layer_idx) );
end
% Merge all the grid position and bounding box into list
grid_cube.num_valid_cube = sum(grid_cube.num_bbox_xy);
grid_cube.bbox_grid_ind_list = cat(1, grid_cube.bbox_grid_ind{:});
assert(grid_cube.num_valid_cube == numel(grid_cube.bbox_grid_ind_list), 'Number of boudning box does not match');
grid_cube.bbox_grid_sub_list = cat(1, grid_cube.bbox_grid_sub{:});
grid_cube.bbox_grid_label_array = nan(grid_size_3d);
grid_cube.bbox_grid_label_array(grid_cube.bbox_grid_ind_list) = 1 : grid_cube.num_valid_cube;
grid_cube.bbox_xyz_mmll_list = cat(1, grid_cube.bbox_xyz_mmll{:});
grid_cube.bbox_xyz_mmxx_list = cat(1, grid_cube.bbox_xyz_mmxx{:});
% Extra fields
grid_cube.gridZ = gridZ;
grid_cube.grid2D = grid2D;
grid_cube.downsample_rate = 1; % Preserved for downsampling grid
toc
end
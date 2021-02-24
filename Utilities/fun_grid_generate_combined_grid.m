function combined_grid = fun_grid_generate_combined_grid(grid_info, combined_size, overlap_size, saveQ)
% fun_generate_combined_grid combines the input grid into grid of larger
% size, which can be used for refining segmentation, annotation, reduce the
% boundary effect of skeletonization. 
% Input: 
%   grid_info: structure, original grid used for segmentation and visualization 
%   combined_size: numerical scalar, number of grids combined in each
%   dimension
%   overlap_size: number of grid the 6-neighbor combined grids pair share
%   in each dimension
% Output: 
%   combined_grid: structure, see below. 
% Debug: debug_generate_combined_bbox_grid.m
% Implemented by Xiang Ji on 12/03/2018
% 
% Modified by Xiang Ji on 03/25/2019
% 1. Update from Generate_combined_bbox_grid.m
% 2. Add overlapping information into the version of the grid 
% Modified by Xiang Ji on 04/18/2019
% 1. Change the name of sub_grid_ind to sub_grid_label and add sub_grid_ind,
% which is the global linear indices of the subgrid cube in the entire
% subgrid ( 240_cube )
if nargin < 4
    saveQ = false;
end
grid_info.grid_size = [grid_info.grid2D.grid_size, grid_info.gridZ.grid_size];
all_grid_sub_list = cat(1, grid_info.bbox_grid_sub{:});
all_grid_ind_list = sub2ind(grid_info.grid_size, all_grid_sub_list(:,1), all_grid_sub_list(:,2), all_grid_sub_list(:,3)); 
all_grid_label_array = zeros(grid_info.grid_size);
all_grid_label_array(all_grid_ind_list) = 1 : numel(all_grid_ind_list);
all_grid_xyz_mmll = cat(1, grid_info.bbox_xyz_mmll{:});
all_grid_xyz_mmxx = cat(1, grid_info.bbox_xyz_mmxx{:});

if isscalar(combined_size)
    tmp_grid_version = sprintf('%s_combined_%d_o_%d', grid_info.version, combined_size, overlap_size);
    combined_size = ones(1,3) * combined_size;
else
    tmp_grid_version = sprintf('%s_c%d_%d_%d_o_%d', grid_info.version, combined_size, overlap_size);
end
%% Generate combined grid for 
combined_grid_paraXY = fun_generate_grid(combined_size(1:2), overlap_size, grid_info.grid2D.grid_size);
combined_grid_paraZ = fun_generate_grid(combined_size(3), overlap_size, grid_info.gridZ.grid_size);

combined_grid = struct;
combined_grid.dataset_name = grid_info.dataset_name;
combined_grid.stack = grid_info.stack;
combined_grid.version = tmp_grid_version;
combined_grid.num_grid_layer = combined_grid_paraZ.grid_size;
combined_grid.block_size = combined_size;
combined_grid.block_overlap = overlap_size;
combined_grid.data_type = grid_info.data_type;
combined_grid.data_size = grid_info.data_size;
combined_grid.voxel_size_um = grid_info.voxel_size_um;

combined_grid.grid_size = [combined_grid_paraXY.grid_size, combined_grid_paraZ.grid_size];
combined_grid.num_valid_cube = [];
combined_grid.num_bbox_xy = zeros(1, combined_grid.num_grid_layer);

combined_grid.layer = 1 : combined_grid_paraZ.grid_size;
combined_grid.sub_grid_label_array = cell(combined_grid.grid_size);
combined_grid.sub_grid_valid_idx = cell(combined_grid.grid_size);
combined_grid.sub_grid_label = cell(combined_grid.grid_size);
combined_grid.sub_grid_ind = cell(combined_grid.grid_size);
combined_grid.sub_grid_sub = cell(combined_grid.grid_size);
combined_grid.sub_grid_bbox_mmll = cell(combined_grid.grid_size);
combined_grid.sub_grid_bbox_mmxx = cell(combined_grid.grid_size);
combined_grid.sub_grid_sub_list = all_grid_sub_list;
combined_grid.sub_grid_ind_array = all_grid_label_array;
combined_grid.bbox_xyz_mmll_pixel = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_xyz_mmxx_pixel = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_xyz_mmll_grid = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_xyz_mmxx_grid = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_xy_valid_mat = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_xy_label_mat = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_xyz_label_array = [];
combined_grid.bbox_grid_sub = cell(1, combined_grid_paraZ.num_block);
combined_grid.bbox_grid_sub_list = [];
combined_grid.bbox_grid_ind_list = [];
combined_grid.bbox_xyz_mmxx_grid_list = [];
combined_grid.bbox_xyz_mmll_grid_list = [];
combined_grid.bbox_xyz_mmxx_pixel_list = [];
combined_grid.bbox_xyz_mmll_pixel_list = [];

combined_grid.gridZ = combined_grid_paraZ;
combined_grid.grid2D = combined_grid_paraXY;
combined_grid.grid_ori = grid_info; 
%% Record information: 
% 1. Parameters for generating the grid
% 2. Global coordinate at 1um isotropic resolution of the minimum
% coordinate in the combined grid; Bounding bbox of the global coordinate
% 3. List of grid subscripts in the combined grid 

for layer_idx = 1 : combined_grid_paraZ.grid_size
    sub_grid_layer_min = combined_grid_paraZ.mmxx(layer_idx, 1);
    sub_grid_layer_max = combined_grid_paraZ.mmxx(layer_idx, 2);
    global_sub_3_min = grid_info.gridZ.mmxx(sub_grid_layer_min, 1);
    global_sub_3_max = grid_info.gridZ.mmxx(sub_grid_layer_max, 2);
    
    
    tmp_bbox_xy_label_mat = zeros(combined_grid_paraXY.grid_size);
    tmp_num_valid_combined_grid = 0;
    tmp_valid_combined_grid_bbox_mmll = zeros(combined_grid_paraXY.num_block, 6);
    tmp_valid_combined_grid_bbox_mmxx = zeros(combined_grid_paraXY.num_block, 6);
    tmp_bbox_mmxx_grid = zeros(combined_grid_paraXY.num_block, 6);
    for bbox_idx = 1 : prod(combined_grid_paraXY.grid_size)
        % grid subscript in the combined grid
        idx_1 = combined_grid_paraXY.mm_idx_pos(bbox_idx,1);
        idx_2 = combined_grid_paraXY.mm_idx_pos(bbox_idx,2);
        % grid subscript in the original grid
        sub_grid_idx_1_min = combined_grid_paraXY.mmxx(bbox_idx,1);
        sub_grid_idx_2_min = combined_grid_paraXY.mmxx(bbox_idx,2);
        sub_grid_idx_1_max = combined_grid_paraXY.mmxx(bbox_idx,3);
        sub_grid_idx_2_max = combined_grid_paraXY.mmxx(bbox_idx,4);
        % grid index in the original grid
        sub_grid_ind_min = sub2ind(grid_info.grid2D.grid_size, sub_grid_idx_1_min, sub_grid_idx_2_min);
        sub_grid_ind_max = sub2ind(grid_info.grid2D.grid_size, sub_grid_idx_1_max, sub_grid_idx_2_max);
        % minimum/maximum coordinate in the data set
        global_sub_1_min = grid_info.grid2D.mmxx(sub_grid_ind_min, 1);
        global_sub_2_min = grid_info.grid2D.mmxx(sub_grid_ind_min, 2);
        global_sub_1_max = grid_info.grid2D.mmxx(sub_grid_ind_max, 3);
        global_sub_2_max = grid_info.grid2D.mmxx(sub_grid_ind_max, 4);
        
        % bounding box in the original grid 
        grid_bbox = [combined_grid_paraXY.mmll(bbox_idx,1:2), combined_grid_paraZ.mmll(layer_idx,1), ...
            combined_grid_paraXY.mmll(bbox_idx,3:4), combined_grid_paraZ.mmll(layer_idx,2)];
        sub_grid_label_array = crop_bbox3(all_grid_label_array, grid_bbox, 'default');
        if any(sub_grid_label_array(:))
            tmp_num_valid_combined_grid = tmp_num_valid_combined_grid + 1;
            tmp_bbox_xy_label_mat(bbox_idx) = tmp_num_valid_combined_grid;
            combined_grid.sub_grid_label_array{idx_1, idx_2, layer_idx} = sub_grid_label_array;
            % Linear indices of the position of the subgrid in 5x5x5 array
            combined_grid.sub_grid_valid_idx{idx_1, idx_2, layer_idx} = find(sub_grid_label_array);
            % 
            combined_grid.sub_grid_label{idx_1, idx_2, layer_idx} = sub_grid_label_array(combined_grid.sub_grid_valid_idx{idx_1, idx_2, layer_idx});
            % Global grid ind of the valid 
            combined_grid.sub_grid_ind{idx_1, idx_2, layer_idx} = all_grid_ind_list(combined_grid.sub_grid_label{idx_1, idx_2, layer_idx}, :);            
            combined_grid.sub_grid_sub{idx_1, idx_2, layer_idx} = all_grid_sub_list(combined_grid.sub_grid_label{idx_1, idx_2, layer_idx}, :);
            combined_grid.sub_grid_bbox_mmll{idx_1, idx_2, layer_idx} = all_grid_xyz_mmll(combined_grid.sub_grid_label{idx_1, idx_2, layer_idx}, :);
            
            % bounding boxes of the original grid that are in this combined
            % grid element
            tmp_bbox_mmxx = all_grid_xyz_mmxx(combined_grid.sub_grid_label{idx_1, idx_2, layer_idx}, :);
            combined_grid.sub_grid_bbox_mmxx{idx_1, idx_2, layer_idx} = tmp_bbox_mmxx;
            
            tmp_combined_grid_bbox_mmxx = [global_sub_1_min, global_sub_2_min, global_sub_3_min, ...
                global_sub_1_max, global_sub_2_max, global_sub_3_max];
            tmp_combined_grid_bbox_mmll = tmp_combined_grid_bbox_mmxx;
            tmp_combined_grid_bbox_mmll(4:6) = tmp_combined_grid_bbox_mmll(4:6) - tmp_combined_grid_bbox_mmll(1:3) + 1;
            tmp_valid_combined_grid_bbox_mmll(tmp_num_valid_combined_grid, :) = tmp_combined_grid_bbox_mmll;
            tmp_valid_combined_grid_bbox_mmxx(tmp_num_valid_combined_grid, :) = tmp_combined_grid_bbox_mmxx;
            
            tmp_bbox_mmxx_grid(tmp_num_valid_combined_grid, :) = [sub_grid_idx_1_min, sub_grid_idx_2_min, sub_grid_layer_min,...
            sub_grid_idx_1_max, sub_grid_idx_2_max, sub_grid_layer_max];
            
        end
    end
    tmp_bbox_mmxx_grid = tmp_bbox_mmxx_grid(1:tmp_num_valid_combined_grid, :);
    tmp_bbox_mmll_grid = tmp_bbox_mmxx_grid;
    tmp_bbox_mmll_grid(:, 4:6) = tmp_bbox_mmxx_grid(:, 4:6) - tmp_bbox_mmxx_grid(:, 1:3) + 1;    
    
    combined_grid.num_bbox_xy(layer_idx) = tmp_num_valid_combined_grid;
    combined_grid.bbox_xyz_mmxx_pixel{layer_idx} = tmp_valid_combined_grid_bbox_mmxx(1 : tmp_num_valid_combined_grid, :);
    combined_grid.bbox_xyz_mmll_pixel{layer_idx} = tmp_valid_combined_grid_bbox_mmll(1 : tmp_num_valid_combined_grid, :);
    combined_grid.bbox_xyz_mmxx_grid{layer_idx} = tmp_bbox_mmxx_grid;
    combined_grid.bbox_xyz_mmll_grid{layer_idx} = tmp_bbox_mmll_grid;   
    
    combined_grid.bbox_xy_label_mat{layer_idx} = tmp_bbox_xy_label_mat;
    combined_grid.bbox_xy_valid_mat{layer_idx} = tmp_bbox_xy_label_mat > 0;
    [valid_combined_grid_idx_1, valid_combined_grid_idx_2] = find(combined_grid.bbox_xy_valid_mat{layer_idx});
    combined_grid.bbox_grid_sub{layer_idx} = [valid_combined_grid_idx_1, valid_combined_grid_idx_2, ones(numel(valid_combined_grid_idx_1), 1) * layer_idx];
end
combined_grid.bbox_grid_sub_list = cat(1, combined_grid.bbox_grid_sub{:});
combined_grid.bbox_grid_ind_list = sub2ind(combined_grid.grid_size, ...
    combined_grid.bbox_grid_sub_list(:, 1), combined_grid.bbox_grid_sub_list(:, 2), ...
    combined_grid.bbox_grid_sub_list(:, 3));
combined_grid.bbox_xyz_mmxx_grid_list = cat(1, combined_grid.bbox_xyz_mmxx_grid{:});
combined_grid.bbox_xyz_mmll_grid_list = cat(1, combined_grid.bbox_xyz_mmll_grid{:});
combined_grid.bbox_xyz_mmxx_pixel_list = cat(1, combined_grid.bbox_xyz_mmxx_pixel{:});
combined_grid.bbox_xyz_mmll_pixel_list = cat(1, combined_grid.bbox_xyz_mmll_pixel{:});

grid_valid_ind = find(cat(3, combined_grid.bbox_xy_valid_mat{:}));
num_valid_grid = numel(grid_valid_ind);
combined_grid.num_valid_cube = num_valid_grid;
bbox_xyz_label_array = zeros(combined_grid.grid_size);
bbox_xyz_label_array(grid_valid_ind) = 1 : num_valid_grid;
combined_grid.bbox_xyz_label_array = bbox_xyz_label_array;
combined_grid.bbox_xyz_label_mat = cell(1, combined_grid.num_grid_layer);
for idx = 1 : combined_grid.num_grid_layer
    combined_grid.bbox_xyz_label_mat{idx} = bbox_xyz_label_array(:,:,idx);
end
%%
if saveQ
    DataManager = FileManager;
    DataManager.write_grid_info(combined_grid, combined_grid.dataset_name, combined_grid.stack, combined_grid.version);
end
function grid_cube = fun_get_downsampled_grid3D_from_grid(grid_cube)

DataManager = FileManager;
required_input_field = {'dataset_name', 'stack', 'mask_version', 'enhance_image_version', 'block_size_1x', 'block_overlap_1x', 'downsample_rate', 'grid_name'};
grid_cube = rmfield(grid_cube, setxor(fieldnames(grid_cube), required_input_field));

if grid_cube.downsample_rate > 1 + eps
    if ~isfield(grid_cube, 'grid_name')
        grid_cube.grid_name = sprintf('%d_cube_downsampled_%dx', grid_cube.block_size_1x(1), grid_cube.downsample_rate);
    else
        grid_cube.grid_name = sprintf('%s_downsampled_%dx', grid_cube.grid_name, grid_cube.downsample_rate);
    end
else
    if ~isfield(grid_cube, 'grid_name')
        grid_cube.grid_name = sprintf('%d_cube', grid_cube.block_size_1x(1));
    end
end

input_check = isfield(grid_cube, required_input_field);
if ~all(input_check)
    need_field_list = required_input_field(~input_check);
    for input_idx = 1 : length(need_field_list)
        warning('The input structure should have field %s', need_field_list{input_idx});
    end
    error('The input structure should contains the following fields: dataset_name, stack, mask_version, enhance_image_version, block_size_1x, block_overlap_x1, downsample_rate');
end



% Clear other field of the structure;


grid_cube.block_size = grid_cube.block_size_1x - ceil(grid_cube.block_overlap_1x * (1 - 1/grid_cube.downsample_rate));
if length(grid_cube.block_size) == 1
    grid_cube.block_xy_size = grid_cube.block_size;
    grid_cube.block_z_size = grid_cube.block_size;
elseif length(grid_cube.block_size) == 3
    % Assume the grid block is a square in xy plane
    grid_cube.block_xy_size = grid_cube.block_size(1);
    grid_cube.block_z_size = grid_cube.block_size(3);
else
    error('Unsupported grid block size');
end

grid_cube.block_overlap = ceil(grid_cube.block_overlap_1x/grid_cube.downsample_rate);
% Load data
mask_info = DataManager.load_mask_info(grid_cube.dataset_name, grid_cube.stack, grid_cube.mask_version);

if ~isfield(mask_info, 'voxel_size_1x')
    warning('Voxel size not specified in mask structure. Assume isotropic');
    grid_cube.voxel_size_1x = [1,1,1]; % Assume isotropic in xy plane. Normalize with respect to voxel size in z direction
end
grid_cube.block_size_resized = [grid_cube.block_size(1:2), round(grid_cube.block_size(3)/grid_cube.voxel_size_1x(1))];
grid_cube.block_z_overlap_resized = round(grid_cube.block_overlap / grid_cube.voxel_size_1x(1));

warning('Assume the brain mask file is written by the Tiff library');
brain_mask = DataManager.load_brain_mask(grid_cube.dataset_name, grid_cube.stack, grid_cube.mask_version);
% figure;
% imshow(brain_mask(:,:,1));
% Load enhanced image and downsample it
data_xy_size = ceil(double(mask_info.brain_mask.ori_bbox2_xy_size)/grid_cube.downsample_rate);
data_sec_num = ceil(double( mask_info.brain_mask.ori_L3/grid_cube.downsample_rate));

%% Generate 2D grid to collect data

grid2D = fun_generate_grid(grid_cube.block_xy_size, grid_cube.block_overlap, data_xy_size);

disp('Record xy bouding boxes that contains the brain.');
grid_xy = struct;
grid_xy(mask_info.brain_mask.L3) = struct;
n_bbox = size(grid2D.mmxx, 1);

tic
parfor test_section = 1 : mask_info.brain_mask.L3
    % test_section = 1;
    im_test_dNx = brain_mask(:,:, test_section);
    % Resize the mask to match the grid
    im_resize = imresize(im_test_dNx, data_xy_size);% Cannot rescale the whole BW image stack because it requires a lot of memory
    % For each bouding box, check if the brain is located inside the
    % bounding box
    list_bbox_contains_sample = ones(1,n_bbox, 'logical');
    mat_bbox_contains_sample = zeros(grid2D.size, 'logical');
    mat_bbox_RIO_ratio = zeros(grid2D.size);
    %     bbox_image_record = cell(grid2D.size);
    for test_bbox_id = 1 : n_bbox
        bbox_data = crop_bbox2(im_resize, grid2D.mmll(test_bbox_id,:), 'default');
        % Compute the ratio of the RIO in the block:
        mat_bbox_RIO_ratio(test_bbox_id) = nnz(bbox_data(:))/ prod(grid2D.ll(test_bbox_id, :));
        % Check if the block contains the region of interest
        list_bbox_contains_sample(test_bbox_id) = mat_bbox_RIO_ratio(test_bbox_id) > 0;
        %         if mat_bbox_RIO_ratio(test_bbox_id)
        %             bbox_image_record{test_bbox_id} = bbox_data;
        %         end
        mat_bbox_contains_sample(test_bbox_id) = list_bbox_contains_sample(test_bbox_id);
    end
    list_valid_bbox_idx = find(list_bbox_contains_sample);
    grid_xy(test_section).grid_size = grid2D.size;
    grid_xy(test_section).section = test_section;
    grid_xy(test_section).bbox_idx= list_valid_bbox_idx;
    grid_xy(test_section).bbox_mmll = grid2D.mmll(list_valid_bbox_idx, :);
    grid_xy(test_section).bbox_mmxx = grid2D.mmxx(list_valid_bbox_idx, :);
    grid_xy(test_section).bbox_mm_idx_pos = grid2D.mm_idx_pos(list_valid_bbox_idx,:);
    grid_xy(test_section).bbox_idx_mat = mat_bbox_contains_sample;
    grid_xy(test_section).bbox_area_ratio_mat = mat_bbox_RIO_ratio;
    grid_xy(test_section).bbox_area_ratio = mat_bbox_RIO_ratio(list_valid_bbox_idx);
    %     grid_xy(test_section).bbox_image = bbox_image_record;
end
toc
%%
disp('Combine 2D bouding box into 3d bounding box');
% gridZ generate grid for the downsampled stack
gridZ = fun_generate_grid(grid_cube.block_z_size,grid_cube.block_overlap, data_sec_num);
num_layer = gridZ.size;


grid_cube.data_xy_size = data_xy_size;
grid_cube.data_sec_num = data_sec_num;

grid_cube.num_grid_layer = num_layer;
grid_cube.num_bbox_xy = zeros(num_layer,1);


grid_cube.mask_info = mask_info;
grid_cube.grid2D = grid2D;
% Initialization
grid_cube.layer = 1: num_layer;
grid_cube.bbox_xy_idx = cell(1, num_layer);
grid_cube.bbox_xy_mmll = cell(1, num_layer);
grid_cube.bbox_xy_mmxx = cell(1, num_layer);
grid_cube.bbox_xy_mm_id_pos = cell(1, num_layer);
grid_cube.bbox_grid_sub = cell(1, num_layer);
grid_cube.bbox_xy_valid_mat = cell(1, num_layer);
grid_cube.bbox_xy_linear_idx_mat = cell(1, num_layer);
grid_cube.bbox_z_mmll = cell(1, num_layer);
grid_cube.bbox_z_mmxx = cell(1, num_layer);
grid_cube.bbox_xyz_mmll = cell(1, num_layer);
grid_cube.bbox_xyz_mmxx = cell(1, num_layer);
grid_cube.bbox_xyz_mmll_in_layer = cell(1, num_layer);
grid_cube.bbox_xyz_mmxx_in_layer = cell(1, num_layer);
grid_cube.bbox_volume_ratio = cell(1, num_layer);
grid_cube.bbox_volume_ratio_mat = cell(1, num_layer);
grid_cube.gridZ = gridZ;
% grid_cube.bbox_mask = cell(1, num_layer);
grid_cube.bbox_on_boundaryQ = cell(1, num_layer);
grid_cube.bbox_on_boundary_idx = cell(1, num_layer);
disp('Generating 3D grid information');
for layer_idx = 1 : num_layer
    % Determine which cube in this layer contains the sample.
    sec_1 = max(floor(( gridZ.mmxx(layer_idx,1) - 1 ) * grid_cube.downsample_rate / mask_info.scale_ratio_z ) + 1,1);
    sec_2 = min(floor(( gridZ.mmxx(layer_idx,2) - 1 ) * grid_cube.downsample_rate / mask_info.scale_ratio_z ) + 1, mask_info.brain_mask.L3);
    %     layer_mask1 = logical(imresize3(uint8(cat(3, grid_xy(sec_1:sec_2).bbox_image)), [grid_cube.data_xy_size, gridZ.ll(layer_idx)]));
    
    cube_volume_ratio = mean(cat(3, grid_xy(sec_1: sec_2).bbox_area_ratio_mat),3);
    cube_contains_sampleQ = cube_volume_ratio>0;
    list_valid_cube_idx = find(cube_contains_sampleQ);
    num_valid_cube_in_layer = nnz(cube_contains_sampleQ);
    % Linear index in layer
    grid_cube.bbox_xy_idx{layer_idx} = list_valid_cube_idx;
    % bounding box parameter in layer
    grid_cube.bbox_xy_mmll{layer_idx} = grid2D.mmll(list_valid_cube_idx, :);
    grid_cube.bbox_xy_mmxx{layer_idx} = grid2D.mmxx(list_valid_cube_idx, :);
    % 2D index in layer
    grid_cube.bbox_xy_mm_id_pos{layer_idx} = grid2D.mm_idx_pos(list_valid_cube_idx,:);
    grid_cube.bbox_grid_sub{layer_idx} = cat(2, grid2D.mm_idx_pos(list_valid_cube_idx,:), layer_idx .* ones(num_valid_cube_in_layer,1));
    % 2D information in grid matrix form
    grid_cube.bbox_xy_valid_mat{layer_idx} = cube_contains_sampleQ;
    grid_cube.bbox_volume_ratio_mat{layer_idx} = cube_volume_ratio;
    grid_cube.bbox_volume_ratio{layer_idx} = cube_volume_ratio(list_valid_cube_idx);
    grid_cube.bbox_xy_linear_idx_mat{layer_idx} = zeros(size(cube_contains_sampleQ));
    grid_cube.bbox_xy_linear_idx_mat{layer_idx}(cube_contains_sampleQ) = 1 : num_valid_cube_in_layer;
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
    
    % valid blocks on the boundary of the sample
    grid_cube.bbox_on_boundaryQ{layer_idx} = grid_cube.bbox_volume_ratio{layer_idx} < 1;
    grid_cube.bbox_on_boundary_idx{layer_idx} = find(grid_cube.bbox_on_boundaryQ{layer_idx});
    % Collect mask for each blocks on the boundary. Leave the internal
    % block cell empty
    %     grid_cube.bbox_mask{layer_idx} = cell(1, length(list_valid_cube_idx));
    %     for bbox_idx = 1 : num_valid_cube_in_layer
    %         if grid_cube.bbox_on_boundaryQ{layer_idx}(bbox_idx)
    %             num_section = sec_2 - sec_1 + 1;
    %             sec_list = sec_1 : sec_2;
    %             tmpMask = false([grid_cube.bbox_xy_mmll{layer_idx}(bbox_idx, 3:4), num_section]);
    %             tmpPos1 = grid_cube.bbox_xy_mm_id_pos{layer_idx}(bbox_idx,1);
    %             tmpPos2 = grid_cube.bbox_xy_mm_id_pos{layer_idx}(bbox_idx,2);
    %             for sec_idx = 1 : num_section
    %                 if ~isempty(grid_xy(sec_list(sec_idx)).bbox_image{tmpPos1, tmpPos2})
    %                     tmpMask(:,:,sec_idx) = grid_xy(sec_list(sec_idx)).bbox_image{tmpPos1, tmpPos2};
    %                 end
    %             end
    % %             grid_cube.bbox_mask{layer_idx}{bbox_idx} = logical(imresize3(uint8(tmpMask),grid_cube.bbox_xyz_mmll{layer_idx}(bbox_idx, 4:6)));
    %               grid_cube.bbox_mask{layer_idx}{bbox_idx} = tmpMask;
    %         end
    %     end
end
DataManager.write_grid_info(grid_cube, grid_cube.dataset_name, grid_cube.stack, grid_cube.grid_name);
disp('Grid information saved');
end
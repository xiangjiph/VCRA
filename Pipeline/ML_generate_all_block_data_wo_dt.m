%% To do
% 1. Take care of the 0-valued voxels in the rendered image that is due to
% stitching error. Find the region and pad it with the median value of its
% neighbors.
% 2. Add image preprocessing. Saturation lead to line shift artefact in the
% rendered data
set_env;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
% stack = 'ML_2018_08_15';
grid_version = '240_cube';
grid_cube = DataManager.load_grid(dataset_name, stack, grid_version);
data_octree = grid_cube.octree;
tile_size = data_octree.block_size;
grid_voxel_size = grid_cube.voxel_size_um;

target_tile_size = data_octree.target_block_size;
target_tile_num_voxel = prod(target_tile_size);
target_voxel_scale = data_octree.target_voxel_scale;
render_dataset_size = data_octree.data_size;
render_tile_size = data_octree.block_size;

process_layer_list = 1 : grid_cube.num_grid_layer;
large_vessel_min_r_um = 10;
max_background_int = 1.4e4;
fill_0_w_neighbor_backgroundQ = true;
log_file_fp = DataManager.fp_task_file(dataset_name, stack, 'Generate_240_cube', 'Running_log.txt');
%% Check the integrety of the octree files
% for iter_tile = 1 : data_octree.num_block
%     tmp_ind = data_octree.grid_pos_ind;
%     tmp_fp = data_octree.filepath{tmp_ind};
%     if ~isfile(tmp_fp)
%         warning('File %s does not exist.', tmp_fp);
%     end
% end
% fprintf('Finish checking.\n');
% tile_exist_Q = ~cellfun(@isempty, data_octree.filepath);
% implay(tile_exist_Q);
%%
fprintf('%s\nStart generating 240 image cubes\n', datestr(now));
record_t_start = tic;
diary(log_file_fp);
for layer_idx = 23 : process_layer_list(end)
    %     layer_idx = 50;
    diary('on');
    tmp_tic = tic;
    fprintf('%s\tProcessing layer %d\n', datestr(now), layer_idx);
    
    if grid_cube.num_bbox_xy(layer_idx) == 0
        fprintf('Grid layer %d has no valid 240 cubes\n', layer_idx);
        continue;
    end    
    layer_bbox_z_l = grid_cube.bbox_z_mmll{layer_idx}(2);
    % 240-grid bounding boxes in the layer
    layer_bbox_xyz_mmmxxx_list = grid_cube.bbox_xyz_mmxx{layer_idx};
    % bounding box of all the bounding boxes in the layer - in pixel 
    layer_bbox_xyz_mmmxxx = [min(layer_bbox_xyz_mmmxxx_list(:, 1:3),[],1), max(layer_bbox_xyz_mmmxxx_list(:, 4:6),[],1)];
    % Convert the bounding box in pixel to in micron 
    % Pad voxels in each dimension to avoid rounding error
    load_bbox_xyz_mmmxxx = layer_bbox_xyz_mmmxxx;
    load_bbox_xyz_mmmxxx(4:5) = min(grid_cube.data_size(1:2), load_bbox_xyz_mmmxxx(4:5) + 1);
    
    render_min_sub_str = fun_mouselight_grid_sub_to_octree_sub(layer_bbox_xyz_mmmxxx(1:3), grid_voxel_size, data_octree);
    render_max_sub_str = fun_mouselight_grid_sub_to_octree_sub(load_bbox_xyz_mmmxxx(4:6), ...
        grid_voxel_size, data_octree);
    
    tile_cell_array_size = render_max_sub_str.render_grid_sub - render_min_sub_str.render_grid_sub + 1;
    num_tile_xy = tile_cell_array_size(1) * tile_cell_array_size(2);
    num_tile_layer = tile_cell_array_size(3);
    image_data_type = grid_cube.data_type;
    
    [tile_grid_sub_1, tile_grid_sub_2] = ndgrid(render_min_sub_str.render_grid_sub(1) : render_max_sub_str.render_grid_sub(1), ...
        render_min_sub_str.render_grid_sub(2) : render_max_sub_str.render_grid_sub(2));
    tile_images = cell(tile_cell_array_size(1:2));
    % raw_data_info = cell(tile_cell_array_size);
    % block_data_info = cell(grid_cube.grid2D.size);
    parpool(16);
    parfor (tile_idx = 1 : num_tile_xy, 16)
        num_thread_parfor = maxNumCompThreads(3);
        
        tmp_tile_sub1 = tile_grid_sub_1(tile_idx);
        tmp_tile_sub2 = tile_grid_sub_2(tile_idx);
        tmp_tile_image_cell = cell(num_tile_layer, 1);
        for tile_layer = 1 : num_tile_layer
            tmp_tile_sub3 = render_min_sub_str.render_grid_sub(3) + tile_layer - 1;
           
            if tile_layer == 1
                tmp_section_list = render_min_sub_str.render_in_tile_sub(3) : tile_size(3);
            elseif tile_layer == num_tile_layer
                tmp_section_list = 1 : render_max_sub_str.render_in_tile_sub(3);
            else
                % load all
                tmp_section_list = 1 : tile_size(3);
            end
            tmp_fp = data_octree.filepath{tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3};
           
            if ~isempty(tmp_fp)
                if ~isfile(tmp_fp)
                    warning('Tile (%d %d %d) is missing. Use blank array\n', tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3);
                    tmp_image = zeros([render_tile_size(1:2), numel(tmp_section_list)], image_data_type);
                else
                    fprintf('Processing tile (%d %d %d)\n', tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3);
                    tmp_image = DataManager.load_single_tiff(tmp_fp, tmp_section_list);
                end
            else
                tmp_image = zeros([render_tile_size(1:2), numel(tmp_section_list)], image_data_type);
            end
            tmp_tile_image_cell{tile_layer} = tmp_image;
        end
        % Concatenate image stacks
        tmp_image = cat(3, tmp_tile_image_cell{:});
        
        if fill_0_w_neighbor_backgroundQ
            is_zero_Q = (tmp_image == 0);
            if any(is_zero_Q(:))
                tmp_neighbor = imdilate(is_zero_Q, strel('cube', 3));
                tmp_neighbor = tmp_neighbor & ~ is_zero_Q;
                tmp_neighbor_int = tmp_image(tmp_neighbor);
                tmp_neighbor_int = tmp_neighbor_int(tmp_neighbor_int < max_background_int);
                if isempty(tmp_neighbor_int)
                    tmp_neighbor_int_medial = 11500;
                else
                    tmp_neighbor_int_medial = median(tmp_neighbor_int);
                end
                tmp_image(is_zero_Q) = tmp_neighbor_int_medial;
            end
        end        
        tmp_im_size = size(tmp_image);
        if tmp_im_size(3) == target_tile_size(3)
            tmp_image_rz_size = target_tile_size;
            tmp_image_rz_size(3) = layer_bbox_z_l;
        else
            fprintf('The size of the concatenated image stack is (%d, %d, %d). Compute the size...\n', ...
                tmp_im_size);
            tmp_image_rz_size = round(tmp_im_size ./ target_voxel_scale);
        end
        tmp_image_rz = imresize3(tmp_image, tmp_image_rz_size);
        tile_images{tile_idx} = tmp_image_rz;
    end     
    poolobj = gcp('nocreate');
    delete(poolobj);
    tile_images = cell2mat(tile_images);
    %% Crop the tile images and output first, before doing distance transform
    disp('Crop images and write blocks');
    for block_idx = 1 : grid_cube.num_bbox_xy(layer_idx)
        block_bbox_local = grid_cube.bbox_xyz_mmll{layer_idx}(block_idx,:);
        block_bbox_local(1) = block_bbox_local(1) - render_min_sub_str.grid_space_bbox_mm(1) + 1;
        block_bbox_local(2) = block_bbox_local(2) - render_min_sub_str.grid_space_bbox_mm(2) + 1;
        block_bbox_local(3) = 1;
        tmp_data = crop_bbox3(tile_images, block_bbox_local, 'default');
        tmp_idx1 = grid_cube.bbox_xy_pos_sub_in_layer{layer_idx}(block_idx,1);
        tmp_idx2 = grid_cube.bbox_xy_pos_sub_in_layer{layer_idx}(block_idx,2);
        fprintf('Writing block (%d %d) in layer %d\n', tmp_idx1, tmp_idx2, layer_idx);
        DataManager.write_block_data_file(tmp_data, grid_cube.dataset_name, grid_cube.stack, grid_cube.version, ...
            tmp_idx1, tmp_idx2, layer_idx);
    end
    clearvars tile_images tmp_data
    fprintf('Finish processing layer %d. Elapsed time is %f seconds.\n', layer_idx, toc(tmp_tic));
    diary('off');
end
diary('on');
fprintf('Finish task. Elapsed time is %f seconds\n%s', toc(record_t_start), datestr(now));
diary('off');
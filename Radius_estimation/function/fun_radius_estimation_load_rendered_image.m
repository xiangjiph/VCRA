function [tile_images, tile_image_bbox_r_mm] = fun_radius_estimation_load_rendered_image(data_octree, ...
    bbox_mmxx, voxel_size, parallel_load_Q)

if nargin < 4
    parallel_load_Q = false;
end

persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
current_machine_data_root_directory = DataManager.ROOT_PATH;
render_min_sub_str = fun_mouselight_grid_sub_to_octree_sub(bbox_mmxx(1:3), ...
    voxel_size, data_octree);
render_max_sub_str = fun_mouselight_grid_sub_to_octree_sub(bbox_mmxx(4:6), ...
    voxel_size, data_octree);
% Load the rendered data
tile_cell_array_size = render_max_sub_str.render_grid_sub - render_min_sub_str.render_grid_sub + 1;
num_tile_xy = tile_cell_array_size(1) * tile_cell_array_size(2);
num_tile_layer = tile_cell_array_size(3);

[tile_grid_sub_1, tile_grid_sub_2] = ndgrid(render_min_sub_str.render_grid_sub(1) : render_max_sub_str.render_grid_sub(1), ...
    render_min_sub_str.render_grid_sub(2) : render_max_sub_str.render_grid_sub(2));

tile_images = cell(tile_cell_array_size(1:2));
if parallel_load_Q
    parfor (tile_idx = 1 : num_tile_xy, 8)
        %     num_thread_parfor = maxNumCompThreads(24);
        tmp_tile_sub1 = tile_grid_sub_1(tile_idx);
        tmp_tile_sub2 = tile_grid_sub_2(tile_idx);
        tmp_tile_image_cell = cell(num_tile_layer, 1);
        for tile_layer = 1 : num_tile_layer
            tmp_tile_sub3 = render_min_sub_str.render_grid_sub(3) + tile_layer - 1;
            
            if tile_layer == 1
                tmp_section_list = render_min_sub_str.render_in_tile_sub(3) : data_octree.block_size(3);
            elseif tile_layer == num_tile_layer
                tmp_section_list = 1 : render_max_sub_str.render_in_tile_sub(3);
            else
                % load all
                tmp_section_list = 1 : data_octree.block_size(3);
            end
            if isfield(data_octree, 'relative_filepath')
                tmp_fp = fullfile(current_machine_data_root_directory, ...
                    data_octree.relative_filepath{tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3});
            else
                tmp_fp = data_octree.filepath{tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3};
            end
            
            if ~isempty(tmp_fp)
                if ~isfile(tmp_fp)
                    warning('Tile (%d %d %d) is missing. Use blank array\n', tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3);
                    tmp_image = zeros([data_octree.block_size(1:2), numel(tmp_section_list)], data_octree.data_type);
                else
                    tmp_image = DataManager.load_single_tiff(tmp_fp, tmp_section_list);
                end
            else
                tmp_image = zeros([data_octree.block_size(1:2), numel(tmp_section_list)], data_octree.data_type);
            end
            tmp_tile_image_cell{tile_layer} = tmp_image;
        end
        % Concatenate image stacks
        tile_images{tile_idx} = cat(3, tmp_tile_image_cell{:});
    end
else
    for tile_idx = 1 : num_tile_xy
        tmp_tile_sub1 = tile_grid_sub_1(tile_idx);
        tmp_tile_sub2 = tile_grid_sub_2(tile_idx);
        tmp_tile_image_cell = cell(num_tile_layer, 1);
        for tile_layer = 1 : num_tile_layer
            tmp_tile_sub3 = render_min_sub_str.render_grid_sub(3) + tile_layer - 1;
            
            if tile_layer == 1
                tmp_section_list = render_min_sub_str.render_in_tile_sub(3) : data_octree.block_size(3);
            elseif tile_layer == num_tile_layer
                tmp_section_list = 1 : render_max_sub_str.render_in_tile_sub(3);
            else
                % load all
                tmp_section_list = 1 : data_octree.block_size(3);
            end
            if isfield(data_octree, 'relative_filepath')
                tmp_fp = fullfile(current_machine_data_root_directory, ...
                    data_octree.relative_filepath{tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3});
            else
                tmp_fp = data_octree.filepath{tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3};
            end
          
            if ~isempty(tmp_fp)
                if ~isfile(tmp_fp)
                    warning('Tile (%d %d %d) is missing. Use blank array\n', tmp_tile_sub1, tmp_tile_sub2, tmp_tile_sub3);
                    tmp_image = zeros([data_octree.block_size(1:2), numel(tmp_section_list)], data_octree.data_type);
                else
                    tmp_image = DataManager.load_single_tiff(tmp_fp, tmp_section_list);
                end
            else
                tmp_image = zeros([data_octree.block_size(1:2), numel(tmp_section_list)], data_octree.data_type);
            end
            tmp_tile_image_cell{tile_layer} = tmp_image;
        end
        % Concatenate image stacks
        tile_images{tile_idx} = cat(3, tmp_tile_image_cell{:});
    end    
end
tile_images = cell2mat(tile_images);
tile_image_bbox_r_mm = render_min_sub_str.render_space_bbox_mm;
tile_image_bbox_r_mm(3) = tile_image_bbox_r_mm(3) + render_min_sub_str.render_in_tile_sub(3) - 1;
end
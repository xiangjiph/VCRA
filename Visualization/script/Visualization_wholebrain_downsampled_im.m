dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
DataManager = FileManager;
render_tile_info = DataManager.load_dataset_info(dataset_name, stack);

vis_octree_level = 5;
vis_octree_info = render_tile_info.octree(vis_octree_level);

max_z_proj_um = 100;
max_z_proj_num_sec = round(max_z_proj_um / vis_octree_info.voxel_size(3));
max_proj_step_sec = round(max_z_proj_num_sec / 2);
max_z_proj_setp_um = round(max_z_proj_um / 2);

proj_grid = fun_generate_grid(max_z_proj_num_sec, max_proj_step_sec, vis_octree_info.data_size(3));
im_save_folder_name = sprintf('Movmax_proj_z_%d_sec_step_%d_sec', max_z_proj_um, max_z_proj_setp_um);
im_save_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), im_save_folder_name);
num_core = 8;
%%
fprintf('Start generating downsampled max projection visualization...\n');
task_tic = tic;
for iter_z = 1 : proj_grid.num_block
    tmp_tic = tic;
    tmp_z_0 = proj_grid.mmxx_array(1, iter_z);
    tmp_z_1 = proj_grid.mmxx_array(2, iter_z);
    tmp_z_in_tile_0 = mod(tmp_z_0, vis_octree_info.block_size(3));
    if tmp_z_in_tile_0 == 0
        tmp_z_in_tile_0 = vis_octree_info.block_size(3);
    end
    tmp_z_in_tile_1 = mod(tmp_z_1, vis_octree_info.block_size(3));
    if tmp_z_in_tile_1 == 0
        tmp_z_in_tile_1 = vis_octree_info.block_size(3);
    end
    tmp_z_grid_sub_0 = ceil(tmp_z_0 / vis_octree_info.block_size(3));
    tmp_z_grid_sub_1 = ceil(tmp_z_1 / vis_octree_info.block_size(3));
    tmp_load_fp_cell = vis_octree_info.filepath(:, :, tmp_z_grid_sub_0 : tmp_z_grid_sub_1);
    tmp_load_fp_valid_Q = ~cellfun(@isempty, tmp_load_fp_cell);
    tmp_load_fp_valid_Q = any(tmp_load_fp_valid_Q, 3);
    tmp_cell_size = size(tmp_load_fp_cell);
    tmp_num_layer = tmp_z_grid_sub_1 - tmp_z_grid_sub_0 + 1;
    tmp_load_data_cell = cell(tmp_cell_size(1:2));
    tmp_num_cell = numel(tmp_load_data_cell);
%     for iter_cell = 1 : tmp_num_cell
   parfor (iter_cell = 1 : tmp_num_cell, num_core)
        if ~tmp_load_fp_valid_Q(iter_cell)
            tmp_load_data_cell{iter_cell} = zeros(vis_octree_info.block_size(1:2), ...
                vis_octree_info.data_type);
        else
            tmp_cell_data = cell(tmp_num_layer, 1);
            tmp_grid_sub = fun_ind2sub(tmp_cell_size(1:2), iter_cell);
            switch tmp_num_layer
                case 1
                    tmp_cell_data{1} = max(DataManager.load_single_tiff(...
                        tmp_load_fp_cell{iter_cell}, tmp_z_in_tile_0 : tmp_z_in_tile_1), [], 3);
                case 2
                    if ~isempty(tmp_load_fp_cell{tmp_grid_sub(1), tmp_grid_sub(2), 1})
                        tmp_cell_data{1} = max(DataManager.load_single_tiff(...
                            tmp_load_fp_cell{tmp_grid_sub(1), tmp_grid_sub(2), 1},...
                            tmp_z_in_tile_0 : vis_octree_info.block_size(3)), [], 3);
                    end
                    if ~isempty(tmp_load_fp_cell{tmp_grid_sub(1), tmp_grid_sub(2), 2})
                        tmp_cell_data{2} = max(DataManager.load_single_tiff(...
                            tmp_load_fp_cell{tmp_grid_sub(1), tmp_grid_sub(2), 2},...
                            1 : tmp_z_in_tile_1), [], 3);
                    end
                otherwise
                    error('To be implemented');
            end
            tmp_cell_data = cat(3, tmp_cell_data{:});
            tmp_cell_data = max(tmp_cell_data, [], 3);
            tmp_load_data_cell{iter_cell} = tmp_cell_data;
        end
    end
    tmp_load_data_cell = cell2mat(tmp_load_data_cell);
    tmp_save_file_path = fullfile(im_save_folder, ...
        sprintf('%s_%s_%s_section_%d.tiff', dataset_name, stack, ...
        im_save_folder_name, iter_z));
    tmp_load_data_cell = im2uint8(tmp_load_data_cell);
    DataManager.write_tiff_stack(tmp_load_data_cell, tmp_save_file_path);
    fprintf('Finish writing %s. Elapsed time is %f seconds.\n', ...
        tmp_save_file_path, toc(tmp_tic));
end
fprintf('Finish generating downsampled visualization. Elapsed time is %f seconds.\n', ...
    toc(task_tic));
function raw_data_grid = fun_mouselight_get_raw_data_grid(raw_data_path, channel, write_fp)
% fun_mouselight_get_raw_data_grid collect raw data file path and
% their corresponding miscroscope stage information and organize them into
% a structure; 
% Input: 
%   raw_data_path: string, root path of the raw data folder, which contains folders
%   YYYY-MM-DD/**/*.tif
%   channel: numerical scalar. image channel. Default value is 0.
%   write_fp: string. If not empty, write the raw data grid into that
%   directory
if nargin < 2
    channel = 0;
    write_fp = [];
elseif nargin < 3
    write_fp = [];
end
% Gather raw images file paths
match_filename_pattern = sprintf('**/*.%d.tif', channel);
disp('Collect data file list');
raw_data_info = dir(fullfile(raw_data_path, match_filename_pattern));
num_file = numel(raw_data_info);
% Test if the fov_x_size_um and x_size_um are redundant 
disp('Load microscope information');
% Might take a few minute without parallelism
for tile_idx = 1 : numel(raw_data_info)
    tmp_stage_info = fun_io_load_microscope_info(raw_data_info(tile_idx).folder);
    tmp_field_name  = fieldnames(tmp_stage_info);
    for field_idx = 1 : numel(tmp_field_name)
        raw_data_info(tile_idx).(strip(tmp_field_name{field_idx})) = tmp_stage_info.(tmp_field_name{field_idx});
    end
end
%% Determine the grid position of the valid tile
disp('Generate raw data grid');
raw_data_grid = struct;
grid_sub_1 = [raw_data_info.y];
grid_sub_2 = [raw_data_info.x];
grid_sub_3 = [raw_data_info.z];
stage_x_um = [raw_data_info.x_mm];
stage_y_um = [raw_data_info.y_mm];
stage_z_um = [raw_data_info.z_mm];
grid_origin = [min(grid_sub_1), min(grid_sub_2), min(grid_sub_3)];
grid_sub_max = [max(grid_sub_1), max(grid_sub_2), max(grid_sub_3)];
grid_size = grid_sub_max - grid_origin + 1;
grid_sub = bsxfun(@minus, [grid_sub_1;grid_sub_2;grid_sub_3], grid_origin') + 1;

grid_imaged_linear_idx = zeros(grid_size);
raw_data_grid.stage_xyz_um = zeros([3, grid_size]);
raw_data_grid.grid_size = grid_size;
raw_data_grid.num_tile = num_file;
raw_data_grid.image_filepath_array = cell(raw_data_grid.grid_size);

for tile_idx = 1 : num_file
    if ~grid_imaged_linear_idx(grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx))
        grid_imaged_linear_idx(grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx)) = tile_idx;
        raw_data_grid.image_filepath_array{grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx)} = fullfile(raw_data_info(tile_idx).folder, ...
            raw_data_info(tile_idx).name);
        raw_data_grid.stage_xyz_um(:, grid_sub(1,tile_idx), grid_sub(2, tile_idx), grid_sub(3, tile_idx)) = [stage_x_um(tile_idx),...
            stage_y_um(tile_idx), stage_z_um(tile_idx)];
    else
        error('Duplicated tile grid coordinate');
    end
end
raw_data_grid.linear_idx_array = grid_imaged_linear_idx;
raw_data_grid.grid_pos_ind = find(raw_data_grid.linear_idx_array>0);
raw_data_grid.grid_pos_sub = fun_ind2sub(grid_size, raw_data_grid.grid_pos_ind);
raw_data_grid.image_filepath = raw_data_grid.image_filepath_array(raw_data_grid.grid_pos_ind);

raw_data_grid.image_filepath_array(raw_data_grid.grid_pos_ind) = raw_data_grid.image_filepath;
raw_data_grid.FOV_size_um = [raw_data_info(1).fov_y_size_um, raw_data_info(1).fov_x_size_um, 251];
raw_data_grid.block_size = [1536, 1024, 251];
raw_data_grid.voxel_size = raw_data_grid.FOV_size_um ./ raw_data_grid.block_size;
raw_data_grid.overlap_size_um = [raw_data_info(1).fov_y_overlap_um, ...
    raw_data_info(1).fov_x_overlap_um, raw_data_info(1).fov_z_overlap_um];
raw_data_grid.overlap_size = round(raw_data_grid.overlap_size_um ./ raw_data_grid.voxel_size);
raw_data_grid.stage_grid_xyz = [grid_sub_2; grid_sub_1; grid_sub_3]';

if ~isempty(write_fp) && ischar(write_fp)
    tmp_folder = fileparts(write_fp);
    if ~isfolder(tmp_folder)
        mkdir(tmp_folder);
    end
    save(write_fp, '-struct', 'raw_data_grid');
    fprintf('Write structure to %s\n', write_fp);
end

end
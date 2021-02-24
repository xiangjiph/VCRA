function data_info = fun_mouselight_get_render_data_info(octree_root_path, channel, local_save_fp)

if nargin < 3
    local_save_fp = [];
end
save_fp = fullfile(octree_root_path, sprintf('octree_data_info_ch%d.mat', channel));
% Load information from jl file
calculated_parameter_fp = fullfile(octree_root_path, 'calculated_parameters.jl');
if ~isfile(calculated_parameter_fp)
    error('%s does not exist', calculated_parameter_fp);
end
fileID = fopen(calculated_parameter_fp, 'r');
parameter_file_txt = textscan(fileID, '%s');
fclose(fileID);
% Hard code here
num_level = str2num(parameter_file_txt{1}{8}) + 1;
voxel_size = str2num(parameter_file_txt{1}{20});
tmp_folder_level = strsplit(octree_root_path, '/');
dataset_name = tmp_folder_level{5};
stack = tmp_folder_level{6};
% render_tile_size = str2num(parameter_file_txt{1}{16});
% num_channel = str2num(parameter_file_txt{1}{12});
%% 
match_template = sprintf('**/*.%d.tif', channel);
raw_data_info = dir(fullfile(octree_root_path, match_template));
% save('./Utilities/mouselight_file_str.mat','raw_data_info');
for idx = 1 : numel(raw_data_info)
    test_folder_path = raw_data_info(idx).folder;
    octree_coordinate = strsplit(erase(test_folder_path, octree_root_path), '/');
    octree_coordinate = cellfun(@str2num,  octree_coordinate(2:end));
    raw_data_info(idx).octree_coordinate = octree_coordinate;
    raw_data_info(idx).level = length(octree_coordinate);
end
clear raw_data_at_level
raw_data_at_level(num_level) = struct;
for idx = 1 : num_level
    raw_data_at_level(idx).ind = find([raw_data_info.level] == (idx - 1));
end
clear data_info
data_info.dataset_name = dataset_name;
data_info.stack = stack;
data_info.num_level = num_level;
data_info.octree(num_level) = struct;
for idx = 1 : num_level
    downsample_power = idx - 1;
    data_info.octree(idx).filepath = cell(2^downsample_power, 2^downsample_power, 2^downsample_power);    
    data_info.octree(idx).grid_size = ones(1,3) * (2^downsample_power);
    data_info.octree(idx).voxel_size = voxel_size .* 2^(num_level - downsample_power -1);
    data_info.octree(idx).downsample_times = 2^(num_level - idx);
end
disp('Getting filepaths for each level of octree');
for tmp_level = 1 : num_level
    tmp_num_tile = numel(raw_data_at_level(tmp_level).ind);
    tmp_cell_array_size = size(data_info.octree(tmp_level).filepath);
    
    data_info.octree(tmp_level).grid_pos_ind_array = zeros(tmp_cell_array_size);
    for idx = 1 : tmp_num_tile
        raw_data_idx = raw_data_at_level(tmp_level).ind(idx);
        tmp_octree_coordinate = raw_data_info(raw_data_idx).octree_coordinate;
        [~, tmp_ind] = fun_io_octree_coordinate_to_array_sub(tmp_octree_coordinate);
        data_info.octree(tmp_level).filepath{tmp_ind} = sprintf('%s/%s', raw_data_info(raw_data_idx).folder, raw_data_info(raw_data_idx).name);
        data_info.octree(tmp_level).grid_pos_ind_array(tmp_ind) = idx;
    end
    data_info.octree(tmp_level).grid_pos_ind = find(data_info.octree(tmp_level).grid_pos_ind_array>0);
    tmp_grid_sub = zeros(tmp_num_tile,3);
    [tmp_grid_sub(:,1), tmp_grid_sub(:,2), tmp_grid_sub(:,3)] = ind2sub(tmp_cell_array_size, data_info.octree(tmp_level).grid_pos_ind);
    data_info.octree(tmp_level).grid_pos_sub = tmp_grid_sub;
    data_info.octree(tmp_level).num_block = tmp_num_tile;
    data_info.octree(tmp_level).exist_mat = ~cellfun(@isempty, data_info.octree(tmp_level).filepath);
    data_info.octree(tmp_level).valid_ind = find(data_info.octree(tmp_level).exist_mat);
    % Get the size of the block at each level 
    tmp_im_info = imfinfo(data_info.octree(tmp_level).filepath{data_info.octree(tmp_level).valid_ind(1)});
    tmp_size_3 = numel(tmp_im_info);
    if tmp_size_3 > 1
        tmp_im_info = tmp_im_info(1);
    end
    if strcmp(tmp_im_info.SampleFormat, 'Unsigned integer')
        switch tmp_im_info.BitDepth
            case 16
                tmp_data_type = 'uint16';
            case 8
                tmp_data_type = 'uint8';
            case 32
                tmp_data_type = 'uint32';
        end
    else
        error('Unrecognized data type.');
    end
    data_info.octree(tmp_level).block_size = [tmp_im_info.Height, tmp_im_info.Width, tmp_size_3];
    data_info.octree(tmp_level).data_type = tmp_data_type;
    data_info.octree(tmp_level).data_size = data_info.octree(tmp_level).block_size .* data_info.octree(tmp_level).grid_size;    
end
save(save_fp, '-struct', 'data_info');
if ~isempty(local_save_fp)
    [tmp_folder, ~, ~] = fileparts(local_save_fp);
    if ~isfolder(tmp_folder)
        mkdir(tmp_folder);
    end
    save(local_save_fp, '-struct', 'data_info');
end
end
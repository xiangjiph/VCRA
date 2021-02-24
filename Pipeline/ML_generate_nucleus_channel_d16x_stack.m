%% Collect octree file organization information
set_env;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
data_info_fp = DataManager.fp_dataset_info(dataset_name, stack);
if isfile(data_info_fp)
    data_info = DataManager.load_dataset_info(dataset_name, stack);
else
    disp('Collect octree tile organization information');
    rendered_data_path = DataManager.fp_raw_data_folder(dataset_name, stack);
    % The file exist on the server
    if ~isfolder(rendered_data_path)
        rendered_data_path = strrep(rendered_data_path, DataManager.ROOTPATH, DataManager.ServerRootPath);
    end
    channel = 0;
    tic
    data_info = fun_mouselight_get_render_data_info(rendered_data_path, 0, data_info_fp);
    toc
end
%% Load downsampled nucleus image stacks
mask_image_level = 3;
disp('Load downsampled images to compute the mask');
image_overview = cell(data_info.octree(mask_image_level).grid_size);
num_octree_tile = prod(data_info.octree(mask_image_level).grid_size);
for idx = 1 : num_octree_tile
    fprintf('Loading octree tile %d\n', idx);
    tmp_filepath = data_info.octree(mask_image_level).filepath{idx};
    % Get the directory of the nucleus channel
    if ~isempty(tmp_filepath)
        assert(isfile(tmp_filepath), sprintf('%s does not exist', tmp_filepath));
        tmp_filepath = strrep(tmp_filepath, 'default.0.tif', 'default.1.tif');
        assert(isfile(tmp_filepath), sprintf('%s does not exist', tmp_filepath));
        image_overview{idx} = DataManager.load_single_tiff(tmp_filepath);
    else
        image_overview{idx} = zeros(data_info.octree(mask_image_level).block_size, data_info.octree(mask_image_level).data_type);
    end
end
image_overview = cell2mat(image_overview);
%% Write downsampled 16x image and grid visualization to disk
% Generate the 16um isotropic voxel size image
overview_rz_size = round(data_info.octree(mask_image_level).data_size .* data_info.octree(mask_image_level).voxel_size ./ 16);
image_rz = imresize3(image_overview, overview_rz_size);
d16x_image_fp = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x_ch1.tiff');
DataManager.write_tiff_stack(image_rz, d16x_image_fp);
fprintf('Finish writing the d16x image to %s\n', d16x_image_fp);
%% Collect octree file organization information
set_env;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
% stack = 'ML_2019_01_24';
% stack = 'ML20190124';
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
%%
linear_scaling_factor = 1.0926;
% Rename the folder 
processed_folder_name_old = DataManager.fp_processed_data(dataset_name, stack);
processed_folder_name_new = sprintf('%s_0', processed_folder_name_old);
assert(isfolder(processed_folder_name_old) && ~isfolder(processed_folder_name_new))
system(sprintf('mv %s %s', processed_folder_name_old, processed_folder_name_new));
% Modify the voxel size in the data_info and resave it processed_data
for  iter_level = 1 : data_info.num_level
    data_info.octree(iter_level).voxel_size = data_info.octree(iter_level).voxel_size .* linear_scaling_factor;
end
mkdir(processed_folder_name_old);
save(data_info_fp, '-struct', 'data_info');
%% Parameters
set_voxel_size = [1,1,1];
block_size = 240;
block_overlap_l = 32;
grid_name = '240_cube';
% level 7 for 0.25 x 0.25 x 1 um;
% level 5 for 1 x 1 x 4 um. 
image_octree_level = 7; 
data_type = data_info.octree(image_octree_level).data_type;
raw_dataset_size = data_info.octree(image_octree_level).data_size;
raw_data_voxel_size = data_info.octree(image_octree_level).voxel_size;

raw_tile_size_pxl = data_info.octree(image_octree_level).block_size;
raw_tile_size_um = raw_tile_size_pxl .* raw_data_voxel_size;
% Downsampling strategy: concatenate raw image tiles in z and rescale the
% tile to as close to [1,1,1] um voxel size as possible. 
% 1. Compute the number of raw tile sections 
num_raw_z_sec_for_block = block_size * set_voxel_size(3) / raw_data_voxel_size(3);
raw_tile_target_size_3 = round(num_raw_z_sec_for_block);
computed_voxel_size_z = raw_tile_target_size_3 * raw_data_voxel_size(3) / block_size;
% 2. Compute the target size of the raw tile on xy direction
raw_tile_target_size_12 = round(raw_tile_size_um(1:2) ./ set_voxel_size(1:2));
computed_voxel_size_xy = raw_tile_size_um(1:2) ./ raw_tile_target_size_12(1:2);
% Compute the size of the dataset
computed_voxel_size = [computed_voxel_size_xy, computed_voxel_size_z];
target_tile_size = [raw_tile_target_size_12, raw_tile_target_size_3];

voxel_size_scale = computed_voxel_size ./ raw_data_voxel_size;
target_data_size = ceil(raw_dataset_size .* raw_data_voxel_size ./ computed_voxel_size);
target_voxel_size = computed_voxel_size;

%% Load downsampled image to compute the mask 
mask_image_level = 3;
disp('Load downsampled images to compute the mask');
image_overview = cell(data_info.octree(mask_image_level).grid_size);
num_octree_tile = prod(data_info.octree(mask_image_level).grid_size);
for idx = 1 : num_octree_tile
    fprintf('Loading octree tile %d\n', idx);
    tmp_filepath = data_info.octree(mask_image_level).filepath{idx};
    if ~isempty(tmp_filepath)
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
d16x_image_fp = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
DataManager.write_tiff_stack(image_rz, d16x_image_fp);
fprintf('Finish writing the d16x image to %s\n', d16x_image_fp);
%% Load the previous annotated mask and resize it using nearest neighbor interpolation
annotated_mask_fp = fullfile(DataManager.fp_processed_data(dataset_name, stack), ...
    'whole_brain_d16x_annotated_mask.nii.gz');
annotated_mask = niftiread(annotated_mask_fp);
annotated_mask = imresize3(annotated_mask, overview_rz_size, 'Method', 'nearest');

itk_fp = DataManager.fp_constructor(dataset_name, stack, 'processed_data', 'whole_brain_d16x_annotated');
DataManager.visualize_itksnap(image_rz, annotated_mask, itk_fp, true);
%% Mask for registration - remove sinus
% Remove the sinus, clean up the mask and save for registration
% registration_mask = niftiread(annotated_mask_fp);
registration_mask = annotated_mask;
registration_mask = (registration_mask == 1);
reg_cc = bwconncomp(registration_mask);
reg_cc_size = cellfun(@numel, reg_cc.PixelIdxList);
registration_mask(cat(1, reg_cc.PixelIdxList{reg_cc_size < 1e3})) = false;
registration_mask = imclose(registration_mask, strel('sphere', 3));
DataManager.write_brain_mask(registration_mask, dataset_name, stack, 'whole_brain_d16x_annotated_registration');
%% Mask for generating grid - need large surface vessels. 
grid_mask = annotated_mask > 0;
am_cc = bwconncomp(grid_mask);
am_cc_size = cellfun(@numel, am_cc.PixelIdxList);
grid_mask(cat(1, am_cc.PixelIdxList{am_cc_size < 1e3})) = false;
grid_mask = imclose(grid_mask, strel('sphere', 3));
DataManager.visualize_itksnap(image_rz, grid_mask, itk_fp, true);
DataManager.write_brain_mask(grid_mask, dataset_name, stack, 'whole_brain_d16x_annotated');
%% Generate grid and add information
% Use annotated mask 
brain_mask = grid_mask;
% brain_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated');
grid_info = fun_get_grid_from_mask(brain_mask, block_size, block_overlap_l, target_data_size);
grid_info.dataset_name = dataset_name;
grid_info.stack = stack;
grid_info.version = grid_name;
grid_info.data_type = data_type;
grid_info.voxel_size_um = target_voxel_size;

grid_info.octree = data_info.octree(image_octree_level);
grid_info.octree.target_block_size = target_tile_size;
grid_info.octree.target_voxel_scale = voxel_size_scale;
grid_info.octree.clearing_linear_correction_factor = linear_scaling_factor;
DataManager.write_grid_info(grid_info, grid_info.dataset_name, grid_info.stack, grid_info.version);
disp('Finish writing grid information');
% fun_vis_grid(dataset_name, stack, grid_name, image_rz, 16, true);
%% Generate c5_o1 and c5_o2 grid
grid_c5o1 = fun_grid_generate_combined_grid(grid_info, 5, 1, true);
grid_c5o2 = fun_grid_generate_combined_grid(grid_info, 5, 2, true);
fun_vis_grid(dataset_name, stack, grid_c5o1.version, image_rz, 16, true);
fun_vis_grid(dataset_name, stack, grid_c5o2.version, image_rz, 16, true);
%% Generate image blocks for d16x image
% grid_name = '240_cube_d16x';
% target_data_size = data_info.octree(mask_image_level).data_size ./ voxel_size_scale;
% target_voxel_size = data_info.octree(mask_image_level).voxel_size .* voxel_size_scale;
% 
% grid_info = fun_get_grid_from_mask(brain_mask_refine, block_size, block_overlap_l, target_data_size);
% grid_info.dataset_name = dataset_name;
% grid_info.stack = stack;
% grid_info.version = grid_name;
% grid_info.data_type = data_type;
% grid_info.voxel_size_um = target_voxel_size;
% 
% grid_info.octree = data_info.octree(mask_image_level);
% grid_info.octree.target_block_size = target_tile_size;
% grid_info.octree.target_voxel_scale = voxel_size_scale;
% DataManager.write_grid_info(grid_info, grid_info.dataset_name, grid_info.stack, grid_info.version);
% % Downsample xy 
% image_stack = imresize3(image_overview, grid_info.data_size);
% fun_generate_block_data_from_image_stack(image_stack, grid_info);
% 
% disp('Wrote grid information');
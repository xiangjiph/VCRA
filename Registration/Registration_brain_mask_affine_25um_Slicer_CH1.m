%% Load the Allen Atlas
set_env;
DataManager = FileManager;
%% Process allen atlas
atlas_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_25.nrrd');
annotation_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd');
atlas_im = nrrdread(atlas_fp);
atlas_annotation = nrrdread(annotation_fp);
allen_atlas_annotation_name = readtable(fullfile(DataManager.ROOTPATH, 'Allen_atlas', '/structure_tree_safe_2017.csv'));
%% Update the table with new structure label
% Remove unused variables
allen_atlas_annotation_name = removevars(allen_atlas_annotation_name, {'atlas_id',...
    'st_level', 'ontology_id', 'hemisphere_id', 'weight', 'graph_id', ...
    'graph_order', 'color_hex_triplet', 'neuro_name_structure_id', ...
    'neuro_name_structure_id_path', 'failed', 'sphinx_id', ...
    'structure_name_facet', 'failed_facet', 'parent_structure_id'});
assert(all(allen_atlas_annotation_name.id > 0));
% Map id: 
atlas_id_old = allen_atlas_annotation_name.id;
atlas_id_old_max = max(atlas_id_old);
num_atlas_id = numel(atlas_id_old);
assert(num_atlas_id == numel(unique(atlas_id_old)), 'atlas_id is not unique');
map_atlas_id_to_new_id = sparse(atlas_id_old, 1, 1 : num_atlas_id, atlas_id_old_max, 1);
atlas_structure_path_old = allen_atlas_annotation_name.structure_id_path;
% Convert the old id in structure path to new id
atlas_structure_path = cellfun(@(x) str2double(strsplit(x, '/')), atlas_structure_path_old, 'UniformOutput', false);
atlas_structure_path = cellfun(@(x) x(2:end-1), atlas_structure_path, 'UniformOutput', false);
atlas_structure_path = cellfun(@(x) full(map_atlas_id_to_new_id(x)), atlas_structure_path, 'UniformOutput', false);
% Update ID:
allen_atlas_annotation_name.id_old = allen_atlas_annotation_name.id;
allen_atlas_annotation_name.structure_id_path = [];
allen_atlas_annotation_name.id = full(map_atlas_id_to_new_id(allen_atlas_annotation_name.id));
%% Map allen annotation index to continuous positive integer index
% Relabeled annotation
atlas_mask = atlas_annotation > 0;
if num_atlas_id <= 255
    relabel_array_class = 'uint8';
elseif num_atlas_id <= intmax('uint16')
    relabel_array_class = 'uint16';
elseif num_atlas_id <= intmax('uint32')
    relabel_array_class = 'uint32';
else
    relabel_array_class = 'uint64';
end
atlas_annotation_rl = zeros(size(atlas_annotation), relabel_array_class);
atlas_annotation_rl(atlas_mask) = cast(full(map_atlas_id_to_new_id(atlas_annotation(atlas_mask))), relabel_array_class);

%% Save processed atlas information
aa_rl_str = struct;
aa_rl_str.structure_table = allen_atlas_annotation_name;
aa_rl_str.map_old_id_to_new_id = map_atlas_id_to_new_id;
aa_rl_str.map_new_id_to_old_id = atlas_id_old;
aa_rl_str.num_atlas_id = num_atlas_id;
aa_rl_str.structure_path = atlas_structure_path;
aa_rl_str.structure_labeled_array = atlas_annotation_rl;
aa_rl_str.atlas_image = atlas_im;
aa_rl_str.voxel_size_um = [25, 25, 25];
aa_rl_str.atlas_image_fp = atlas_fp;
aa_rl_str.atlas_annotation_fp = annotation_fp;
aa_rl_str.filepath = fullfile(DataManager.fp_Dataset('Allen_atlas'), 'relabeled_structure_data.mat');
DataManager.write_data(aa_rl_str.filepath, aa_rl_str);
%% Load lowest resolution mouselight image
dataset_name = 'WholeBrain';
stack = 'ML20190124';

mask_fp = DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration');
itk_fp = fullfile(mask_fp, sprintf('%s_%s_d16x_registration', dataset_name, stack));
file_path = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
% file_path = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x_ch1.tiff');
image_stack = DataManager.load_single_tiff(file_path);
% image_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated_registration');
wb_im_voxel_size_um = [16, 16, 16];
target_im_voxel_size_um = [25, 25, 25];
%% Further annotate the mask for more precise registration
% DataManager.visualize_itksnap(image_stack, image_mask, itk_fp, true);
% Set the range of downsampled image for registration 
% switch stack
%     case 'ML_2018_08_15'
%         registration_sec = 1 : 780;
%         allen_zero_sec = 1 : 99;
%     case {'ML_2019_01_24', 'ML20190124'}
%         registration_sec = 1 : 828;
%         allen_zero_sec = 1 : 99;
%     case 'ML20200201'
%         registration_sec = 70 : 950;
%         allen_zero_sec = [];
% end
%% Clean up the mask 
image_mask = niftiread(sprintf('%s_mask.nii.gz', itk_fp)) > 0;
im_cc = bwconncomp(image_mask);
im_cc_vol_ratio = cellfun(@numel, im_cc.PixelIdxList);
im_cc_vol_ratio = im_cc_vol_ratio ./ sum(im_cc_vol_ratio);
for iter_cc = 1 : numel(im_cc_vol_ratio)
    if im_cc_vol_ratio(iter_cc) < 1e-3
        image_mask(im_cc.PixelIdxList{iter_cc}) = false;
    end
end        
%% Remove the image outside the mask 
image_stack_masked = image_stack;
empty_int = median(image_stack(image_mask & ~ imerode(image_mask, strel('sphere', 2))));
image_stack_masked(~image_mask) = empty_int;
% DataManager.visualize_itksnap(image_stack_masked, image_mask);
% Crop the image to fit the range of the Allen Atlas
% im_to_register = image_stack_masked(:, :, registration_sec);
% im_mask_to_register = image_mask(:, :, registration_sec);

im_to_register = image_stack_masked;
im_mask_to_register = image_mask;
% Resize the image and save it as nrrd file
image_size = size(im_to_register);
target_im_size = round(image_size .* wb_im_voxel_size_um ./ target_im_voxel_size_um);
im_to_register = imresize3(im_to_register, target_im_size);
im_mask_to_register = imresize3(uint8(im_mask_to_register), target_im_size) > 0;
% Permute the axes to get the segittal plane to be on the xy plane
im_to_register_sagittal = permute(im_to_register, [1, 3, 2]);
im_mask_to_register_sagittal = permute(im_mask_to_register, [1, 3, 2]);
%% Modify the label to differentiate two hemisphere
allen_im_to_register = atlas_im;
% allen_im_to_register = permute(allen_im_to_register, [1, 3, 2]);
% allen_im_to_register(:, :, allen_zero_sec) = 0;
% allen_im_to_register = permute(allen_im_to_register, [1, 3, 2]);
flip_section = 229;
% Generate the mask for one hemisphere
hemisphere_mask = atlas_annotation > 0;
hemisphere_mask(:, :, flip_section : end) = false;
annotation_size = size(atlas_annotation);
tmp_xy_ind = 1 : 2 : prod(annotation_size(1:2));
tmp_xy_sub = fun_ind2sub(annotation_size(1:2), tmp_xy_ind);
tmp_xyz_ind = sub2ind(annotation_size, tmp_xy_sub(:, 1), tmp_xy_sub(:, 2), ones(size(tmp_xy_sub, 1), 1) .* flip_section);
tmp_xyz_ind = tmp_xyz_ind(atlas_annotation(tmp_xyz_ind) > 0);
hemisphere_mask(tmp_xyz_ind) = true;
% hemisphere_mask = uint8(hemisphere_mask);
%% Save the image as nrrd for registration
atlas_im_fp = fullfile(mask_fp, sprintf('allen_image.tiff'));
DataManager.write_tiff_stack(atlas_im, atlas_im_fp); 
atlas_annotation_fp = fullfile(mask_fp, sprintf('allen_annotation_relabeled.tiff'));
DataManager.write_tiff_stack(atlas_annotation_rl, atlas_annotation_fp);

hemisphere_mask_fp = fullfile(mask_fp, sprintf('allen_atlas_hemisphere_mask.tiff'));
DataManager.write_tiff_stack(uint8(hemisphere_mask), hemisphere_mask_fp);

reg_im_fp = fullfile(mask_fp, sprintf('%s_brain_image_25um.tiff', stack));
DataManager.write_tiff_stack(im_to_register_sagittal, reg_im_fp);
reg_mask_fp = fullfile(mask_fp, sprintf('%s_brain_mask_25um.tiff', stack));
DataManager.write_tiff_stack(uint8(im_mask_to_register_sagittal), reg_mask_fp);

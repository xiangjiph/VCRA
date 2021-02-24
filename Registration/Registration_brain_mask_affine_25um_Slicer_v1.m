%% Load the Allen Atlas
set_env;
DataManager = FileManager;
atlas_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_25.nrrd');
annotation_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd');
atlas_im = nrrdread(atlas_fp);
atlas_annotation = nrrdread(annotation_fp);
%% Load lowest resolution mouselight image
dataset_name = 'WholeBrain';
stack = 'ML20200201';

file_path = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
image_stack = DataManager.load_single_tiff(file_path);
image_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated_registration');
%% Further annotate the mask for more precise registration
mask_fp = DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration');
itk_fp = fullfile(mask_fp, sprintf('%s_%s_d16x_registration', dataset_name, stack));
DataManager.visualize_itksnap(image_stack, image_mask, itk_fp, true);
% Set the range of downsampled image for registration 
switch stack
    case 'ML_2018_08_15'
        registration_sec = 1 : 780;
        allen_zero_sec = 1 : 99;
    case {'ML_2019_01_24', 'ML20190124'}
        registration_sec = 1 : 828;
        allen_zero_sec = 1 : 99;
    case 'ML20200201'
        registration_sec = 57 : 855;
        allen_zero_sec = [];
end
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
%%
% Remove the image outside the mask 
image_stack_masked = image_stack;
empty_int = median(image_stack(image_mask & ~ imerode(image_mask, strel('sphere', 2))));
image_stack_masked(~image_mask) = empty_int;
% DataManager.visualize_itksnap(image_stack_masked, image_mask);
% Resize the image and save it as nrrd file
current_voxel_size = [16, 16, 16];
target_voxel_size = [25, 25, 25];
% DataManager.visualize_itksnap(image_stack, image_mask);
% Crop the image to fit the range of the Allen Atlas
im_to_register = image_stack_masked(:, :, registration_sec);
im_mask_to_register = image_mask(:, :, registration_sec);
% Resize the image
image_size = size(im_to_register);
target_im_size = image_size .* current_voxel_size ./ target_voxel_size;
im_to_register = imresize3(im_to_register, target_im_size);
im_mask_to_register = imresize3(uint8(im_mask_to_register), target_im_size);
% Permute the axes to get the segittal plane to be on the xy plane
im_to_register_sagittal = permute(im_to_register, [1, 3, 2]);
im_mask_to_register_sagittal = permute(im_mask_to_register, [1, 3, 2]);
imshowpair(im_to_register_sagittal(:, :, 100), im_mask_to_register_sagittal(:, :, 100));
%% Generate another copy of image and mask, in which only the region with
% complete structures and good image quality are used for registration. 
switch stack
    case 'ML_2018_08_15'
        zero_sec = 1 : 130;
    case {'ML_2019_01_24', 'ML20190124'}
        zero_sec = 1 : 156;
    case 'ML20200201'
        zero_sec = [];
end
im_to_register_valid = image_stack_masked(:, :, registration_sec);
im_to_register_valid(:, :, zero_sec) = 0;
im_to_register_valid = imresize3(im_to_register_valid, target_im_size);
im_to_register_valid = permute(im_to_register_valid, [1, 3, 2]);

im_mask_to_register_valid = image_mask(:, :, registration_sec);
im_mask_to_register_valid(:, :, zero_sec) = false;
im_mask_to_register_valid = imresize3(uint8(im_mask_to_register_valid), target_im_size);
% Make the xy plane the sagittal plane
im_mask_to_register_valid = permute(im_mask_to_register_valid, [1, 3, 2]);
%% Modify the label to differentiate two hemisphere
allen_im_to_register = atlas_im;
allen_im_to_register = permute(allen_im_to_register, [1, 3, 2]);
allen_im_to_register(:, :, allen_zero_sec) = 0;
allen_im_to_register = permute(allen_im_to_register, [1, 3, 2]);
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
hemisphere_mask = uint8(hemisphere_mask);
%% Save the image as nrrd for registration
reg_im_fp = fullfile(mask_fp, sprintf('%s_brain_image_25um_cropped.tiff', stack));
DataManager.write_tiff_stack(im_to_register_sagittal, reg_im_fp);
reg_mask_fp = fullfile(mask_fp, sprintf('%s_brain_mask_25um_cropped.tiff', stack));
DataManager.write_tiff_stack(im_mask_to_register_sagittal, reg_mask_fp);
% Flip the image intensity since major vessels are dark in Allen's image
% Simply invert the intensity 
im_to_register_sagittal_invert = max(im_to_register_sagittal, [], 'all') - im_to_register_sagittal;
im_to_register_sagittal_invert(im_mask_to_register_sagittal == 0) = 0;
reg_im_fp = fullfile(mask_fp, sprintf('%s_brain_image_25um_cropped_inverted.tiff', stack));
DataManager.write_tiff_stack(im_to_register_sagittal_invert, reg_im_fp);

reg_im_valid_fp = fullfile(mask_fp, sprintf('%s_brain_image_25um_cropped_zeroed.tiff', stack));
reg_mask_valid_fp = fullfile(mask_fp, sprintf('%s_brain_mask_25um_cropped_zeroed.tiff', stack));
DataManager.write_tiff_stack(im_to_register_valid, reg_im_valid_fp);
DataManager.write_tiff_stack(im_mask_to_register_valid, reg_mask_valid_fp);

allen_im_zeroed_fp = fullfile(mask_fp, sprintf('allen_atlas_image_zeroed.tiff'));
hemisphere_mask_fp = fullfile(mask_fp, sprintf('allen_atlas_hemisphere_mask.tiff'));
DataManager.write_tiff_stack(allen_im_to_register, allen_im_zeroed_fp);
DataManager.write_tiff_stack(hemisphere_mask, hemisphere_mask_fp);

% Copy the allen atlas to the target directory
[~, atlas_fn, atlas_fn_ext] = fileparts(atlas_fp);
target_atlas_fp = fullfile(mask_fp, sprintf('%s%s', atlas_fn, atlas_fn_ext));
copyfile(atlas_fp, target_atlas_fp);

[~, atlas_annotation_fn, atlas_annotation_fn_ext] = fileparts(annotation_fp);
copyfile(annotation_fp, fullfile(mask_fp, sprintf('%s%s', atlas_annotation_fn, atlas_annotation_fn_ext)));
%% Load the registration result from 3D Slicer
mask_fp = DataManager.fp_mask_folder(dataset_name, 'ML_2018_08_15', 'whole_brain_d16x_registration');
registered_annotation = DataManager.load_data(fullfile(mask_fp, 'annotation_25_registered.nrrd'));
registered_hemisphere_mask = DataManager.load_data(fullfile(mask_fp, 'allen_atlas_hemisphere_mask_registered.nrrd'));
% Convert back to orignianl d16x coordinate
annotation_moved_coronal = permute(registered_annotation, [1, 3, 2]);
annotation_moved_coronal = imresize3(annotation_moved_coronal, image_size, 'nearest');
% Add the cropped part back to the annotation array
annotation_moved_coronal = cat(3, annotation_moved_coronal, zeros(image_size(1), image_size(2), size(image_stack, 3) - size(annotation_moved_coronal, 3), 'like', annotation_moved_coronal));

hemisphere_mask_coronal = permute(registered_hemisphere_mask, [1, 3, 2]);
hemisphere_mask_coronal = imresize3(hemisphere_mask_coronal, image_size, 'nearest') > 0;
hemisphere_mask_coronal = cat(3, hemisphere_mask_coronal, zeros(image_size(1), image_size(2), size(image_stack, 3) - size(hemisphere_mask_coronal, 3), 'like', hemisphere_mask_coronal));
% DataManager.visualize_itksnap(image_stack, annotation_moved_coronal);
%% Generate registration information
allen_atlas_annotation_name = readtable('/data/Vessel/Allen_atlas/structure_tree_safe_2017.csv');
allen_atlas_str = struct;
allen_atlas_str.registered_label_array = annotation_moved_coronal;
allen_atlas_str.registered_hemisphere_mask = hemisphere_mask_coronal;
allen_atlas_str.structure_table = allen_atlas_annotation_name;
allen_atlas_str.label_2_name = containers.Map(allen_atlas_annotation_name.id, allen_atlas_annotation_name.name);
allen_atlas_str.label_2_acronym = containers.Map(allen_atlas_annotation_name.id, allen_atlas_annotation_name.acronym);
% Convert the string of structure id path to number
allen_atlas_str.structure_id_path = cellfun(@(x) str2double(strsplit(x, '/')), allen_atlas_annotation_name.structure_id_path, 'UniformOutput', false);
allen_atlas_str.structure_id_path = cellfun(@(x) x(2:end-1), allen_atlas_str.structure_id_path, 'UniformOutput', false);
allen_atlas_str.structure_id_depth = cellfun(@numel, allen_atlas_str.structure_id_path) - 1;
allen_atlas_str.info.im_fix = 'whole brain mask @ 25um crop the first 780 sections before permutation';
allen_atlas_str.info.im_move = 'Allen atlas ccf 2017 average_template_25.nrrd';
allen_atlas_str.info.registration = '3D Slicer: rigid and affine of template to mask, followed by BSpline of template to inverted image';
allen_atlas_str.info.fp = DataManager.fp_registration_data(dataset_name, stack, 'Allen_2017_25um_nonrigid.mat');
allen_atlas_str.info.voxel_size_um = [16, 16, 16];
allen_atlas_str.fun_name_to_ind = @(name_list, name) find(~cellfun(@isempty, strfind(name_list, name)));
% Indices conversion 
allen_atlas_str.id = allen_atlas_annotation_name.id;
allen_atlas_str.num_structure = numel(allen_atlas_annotation_name.id);
allen_atlas_str.id_2_ind = sparse(allen_atlas_str.id, ones(size(allen_atlas_str.id)), ...
    1 : allen_atlas_str.num_structure, max(allen_atlas_str.id), 1);
%% Compute the center of mass for each structure
allen_atlas_str.array_size = size(annotation_moved_coronal);
allen_atlas_str.structure_sub_left = zeros(allen_atlas_str.num_structure, 3);
allen_atlas_str.structure_sub_right = zeros(allen_atlas_str.num_structure, 3);
allen_atlas_str.structure_cc_ind = cell(allen_atlas_str.num_structure, 1);
allen_atlas_str.structure_cc_ind_hemisphere = cell(allen_atlas_str.num_structure, 2);
[bin_ind, bin_label] = fun_bin_data_to_idx_list(allen_atlas_str.registered_label_array);

for iter_label = 1 : numel(bin_label)
    tmp_label = bin_label(iter_label);
    if tmp_label > 0
        tmp_ind = bin_ind{iter_label};
        tmp_numel = numel(tmp_ind);
        structure_list_ind  = full(allen_atlas_str.id_2_ind(tmp_label));
        if structure_list_ind <= 0
            warning('Label %d is not in the structure list. Size of the connected component is %d', tmp_label, tmp_numel);
        else
            is_in_left_hemisphere_Q = allen_atlas_str.registered_hemisphere_mask(tmp_ind);
            % Randomly sample one voxels in the structure
            if any(is_in_left_hemisphere_Q)
                tmp_ind_1 = tmp_ind(is_in_left_hemisphere_Q);
                tmp_ind_1 = randsample(tmp_ind_1, 1);
                tmp_sub_1 = fun_ind2sub(allen_atlas_str.array_size, tmp_ind_1);
                allen_atlas_str.structure_sub_left(structure_list_ind, :) = tmp_sub_1;
            else
                warning('Label %d is not in the left hemisphere. Number of pixel in the structure is %d', tmp_label, tmp_numel);
            end
            if ~all(is_in_left_hemisphere_Q)
                tmp_ind_2 = tmp_ind(~is_in_left_hemisphere_Q);
                tmp_ind_2 = randsample(tmp_ind_2, 1);
                tmp_sub_2 = fun_ind2sub(allen_atlas_str.array_size, tmp_ind_2);           
                allen_atlas_str.structure_sub_right(structure_list_ind, :) = tmp_sub_2;
            else
                warning('Label %d is not in the right hemisphere. Number of pixel in the structure is %d', tmp_label, tmp_numel);
            end
            allen_atlas_str.structure_cc_ind{structure_list_ind} = tmp_ind;
            allen_atlas_str.structure_cc_ind_hemisphere{structure_list_ind, 1} = tmp_ind(is_in_left_hemisphere_Q);
            allen_atlas_str.structure_cc_ind_hemisphere{structure_list_ind, 2} = tmp_ind(~is_in_left_hemisphere_Q);            
        end
    end
end
%% Merge the structure path into a matrix
allen_atlas_str.structure_id_path_mat = nan(allen_atlas_str.num_structure, ...
    max(allen_atlas_str.structure_id_depth + 1));
for iter_str = 1 : allen_atlas_str.num_structure
    allen_atlas_str.structure_id_path_mat(iter_str, 1 : (allen_atlas_str.structure_id_depth(iter_str) + 1)) = allen_atlas_str.structure_id_path{iter_str};
end
save(allen_atlas_str.info.fp, '-struct', 'allen_atlas_str');
% Write ITK
itk_name = fullfile(DataManager.fp_itksnap_data_folder(dataset_name, stack, 'registration'), sprintf('%s_brain_mask_registration_Allen_2017_25um', stack));
DataManager.visualize_itksnap(image_stack, annotation_moved_coronal, itk_name, true);
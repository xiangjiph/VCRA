%% Load the Allen Atlas
atlas_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/average_template/average_template_25.nrrd');
annotation_fp = fullfile(DataManager.ROOTPATH, 'Allen_atlas', 'download.alleninstitute.org/informatics-archive/current-release/mouse_ccf/annotation/ccf_2017/annotation_25.nrrd');
atlas_im = nrrdread(atlas_fp);
atlas_annotation = nrrdread(annotation_fp);
%% Load lowest resolution mouselight image
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';

file_path = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
image_stack = DataManager.load_single_tiff(file_path);
image_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated');

% Resize the image and save it as nrrd file
current_voxel_size = [16, 16, 16];
target_voxel_size = [25, 25, 25];
% DataManager.visualize_itksnap(image_stack, image_mask);
% Crop the image to fit the range of the Allen Atlas
im_to_register = image_stack(:, :, 1 : 780);
im_mask_to_register = image_mask(:, :, 1 : 780);
% Resize the image
image_size = size(im_to_register);
target_im_size = image_size .* current_voxel_size ./ target_voxel_size;
im_to_register = imresize3(im_to_register, target_im_size);
im_mask_to_register = imresize3(uint8(im_mask_to_register), target_im_size);
figure;imagesc(im_to_register(:, :, 100));
% Permute the axes to get the segittal plane to be on the xy plane
im_to_register_sagittal = permute(im_to_register, [1, 3, 2]);
im_mask_to_register_sagittal = permute(im_mask_to_register, [1, 3, 2]);
figure;imagesc(im_to_register_sagittal(:, :, 100));
imshowpair(im_to_register_sagittal(:, :, 100), im_mask_to_register_sagittal(:, :, 100));
%% Affine registration
% It seems that use the brian mask for registration works better
% (sadly...)
im_move = atlas_im; 
% im_fix = im_to_register_sagittal;
im_fix = im_mask_to_register_sagittal;
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 500;
affine_reg_info = imregtform(im_move, im_fix, 'affine', optimizer, metric, 'DisplayOptimization', true);
im_moved = imwarp(im_move, affine_reg_info, 'OutputView', imref3d(size(im_fix)));
annotation_moved = imwarp(atlas_annotation, affine_reg_info, 'nearest', 'OutputView', imref3d(size(im_fix)));
%% Convert back to the original coordinate
annotation_moved_coronal = permute(annotation_moved, [1, 3, 2]);
annotation_moved_coronal = imresize3(annotation_moved_coronal, image_size, 'nearest');
% Add the cropped part back to the annotation array
annotation_moved_coronal = cat(3, annotation_moved_coronal, zeros(image_size(1), image_size(2), size(image_stack, 3) - size(annotation_moved_coronal, 3), 'like', annotation_moved_coronal));
% Save registration information
itk_name = fullfile(DataManager.fp_itksnap_data_folder(dataset_name, stack, 'registration'), 'ML_2018_08_15_brain_mask_16um_affine_registration');

DataManager.visualize_itksnap(image_stack, annotation_moved_coronal, itk_name, true);
%% Save the image as nrrd for registration
reg_im_fp = DataManager.fp_metadata_file(dataset_name, stack, [stack, '_brain_image_25um_cropped_780'], '.tiff');
DataManager.write_tiff_stack(im_to_register_sagittal, reg_im_fp);
reg_mask_fp = DataManager.fp_metadata_file(dataset_name, stack, [stack, '_brain_mask_25um_cropped_780'], '.tiff');
DataManager.write_tiff_stack(im_mask_to_register_sagittal, reg_mask_fp);
% Flip the image intensity since major vessels are dark in Allen's image
% Simply invert the intensity 
im_to_register_sagittal_invert = max(im_to_register_sagittal, [], 'all') - im_to_register_sagittal;
im_to_register_sagittal_invert(im_mask_to_register_sagittal == 0) = 0;
% implay(im_to_register_sagittal_invert)
reg_im_fp = DataManager.fp_metadata_file(dataset_name, stack, [stack, '_brain_image_25um_cropped_780_inverted'], '.tiff');
DataManager.write_tiff_stack(im_to_register_sagittal_invert, reg_im_fp);
%% Load the registration result from 3D Slicer
registered_annotation = nrrdread('/home/dklab/Documents/Slicer3D/Allen_25_to_ML_2018_08_15_inverted_v2/template_to_image_BSpline_annotation_resampled_nn.nrrd');
DataManager.visualize_itksnap(im_to_register_sagittal, registered_annotation);
% Convert back to orignianl d16x coordinate
annotation_moved_coronal = permute(registered_annotation, [1, 3, 2]);
annotation_moved_coronal = imresize3(annotation_moved_coronal, image_size, 'nearest');
% Add the cropped part back to the annotation array
annotation_moved_coronal = cat(3, annotation_moved_coronal, zeros(image_size(1), image_size(2), size(image_stack, 3) - size(annotation_moved_coronal, 3), 'like', annotation_moved_coronal));
DataManager.visualize_itksnap(image_stack, annotation_moved_coronal);
%% Generate registration information
allen_atlas_annotation_name = readtable('/data/Vessel/Allen_atlas/structure_tree_safe_2017.csv');
allen_atlas_str = struct;
allen_atlas_str.registered_label_array = annotation_moved_coronal;
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
allen_atlas_str.info.voxel_size = 16;
allen_atlas_str.fun_name_to_ind = @(name_list, name) find(~cellfun(@isempty, strfind(name_list, name)));
% Indices conversion 
allen_atlas_str.id = allen_atlas_annotation_name.id;
allen_atlas_str.num_structure = numel(allen_atlas_annotation_name.id);
allen_atlas_str.id_2_ind = sparse(allen_atlas_str.id, ones(size(allen_atlas_str.id)), ...
    1 : allen_atlas_str.num_structure, max(allen_atlas_str.id), 1);
% Compute the center of mass for each structure
allen_atlas_str.array_size = size(annotation_moved_coronal);
allen_atlas_str.structure_center_sub = zeros(allen_atlas_str.num_structure, 3);
%% Compute the center of mass of the valid structures
[bin_ind, bin_label] = fun_bin_data_to_idx_list(allen_atlas_str.registered_label_array);
bin_ind = bin_ind(bin_label ~= 0);
bin_label = bin_label(bin_label ~= 0);
label_str_sub = cellfun(@(x) mean(fun_ind2sub(allen_atlas_str.array_size, ...
    x), 1, 'omitnan'), bin_ind, 'UniformOutput', false);
label_str_sub = cat(1, label_str_sub{:});

structure_sub_ind = full(allen_atlas_str.id_2_ind(bin_label));
structure_sub_ind_valid_Q = structure_sub_ind > 0;

allen_atlas_str.structure_center_sub(structure_sub_ind(structure_sub_ind_valid_Q), :) = ...
    label_str_sub(structure_sub_ind_valid_Q, :);
%% Merge the structure path into a matrix
allen_atlas_str.structure_id_path_mat = nan(allen_atlas_str.num_structure, ...
    max(allen_atlas_str.structure_id_depth + 1));
for iter_str = 1 : allen_atlas_str.num_structure
    allen_atlas_str.structure_id_path_mat(iter_str, 1 : (allen_atlas_str.structure_id_depth(iter_str) + 1)) = allen_atlas_str.structure_id_path{iter_str};
end
save(allen_atlas_str.info.fp, '-struct', 'allen_atlas_str');
%% Write ITK
itk_name = fullfile(DataManager.fp_itksnap_data_folder(dataset_name, stack, 'registration'), 'ML_2018_08_15_brain_mask_registration_Allen_2017_25um');
DataManager.visualize_itksnap(image_stack, annotation_moved_coronal, itk_name, true);
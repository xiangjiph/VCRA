set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
allen_relabeled_data = DataManager.load_data(fullfile(DataManager.fp_Dataset('Allen_atlas'), ...
    'relabeled_structure_data.mat'));
aa_structure_table = allen_relabeled_data.structure_table;
%% Load vessel image
% mask_fp = DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration');
% itk_fp = fullfile(mask_fp, sprintf('%s_%s_d16x_registration', dataset_name, stack));
% file_path_vsl = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
% image_stack_vsl = DataManager.load_single_tiff(file_path_vsl);
file_path = fullfile(DataManager.fp_processed_data(dataset_name, stack), 'whole_stack_d16x.tiff');
image_stack = DataManager.load_single_tiff(file_path);
% image_mask = DataManager.load_brain_mask(dataset_name, stack, 'whole_brain_d16x_annotated_registration');
wb_im_voxel_size_um = [16, 16, 16];
target_im_voxel_size_um = [25, 25, 25];
image_size = size(image_stack);
%% Load the registration result from 3D Slicer
mask_fp = DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registration_landmark');
registered_annotation_fp = fullfile(mask_fp, sprintf('%s_allen_annotation_relabeled_landmark_registered.nrrd', strrep(stack, '_', '')));
registered_hemisphere_mask_fp = fullfile(mask_fp, sprintf('%s_hemisphere_mask_landmark_registered.nrrd', strrep(stack, '_', '')));
registered_annotation = DataManager.load_data(registered_annotation_fp);
assert(all(registered_annotation <= allen_relabeled_data.num_atlas_id, 'all'));
if allen_relabeled_data.num_atlas_id <= intmax('uint16')
    registered_annotation = uint16(registered_annotation);
end
registered_hemisphere_mask = DataManager.load_data(registered_hemisphere_mask_fp);
% Convert back to orignianl d16x coordinate
annotation_coronal = permute(registered_annotation, [1, 3, 2]);
annotation_coronal = imresize3(annotation_coronal, image_size, 'nearest');
% Add the cropped part back to the annotation array
% annotation_moved_coronal = cat(3, annotation_moved_coronal, zeros(image_size(1), image_size(2), size(image_stack, 3) - size(annotation_moved_coronal, 3), 'like', annotation_moved_coronal));
hemisphere_mask_coronal = permute(registered_hemisphere_mask, [1, 3, 2]);
hemisphere_mask_coronal = imresize3(hemisphere_mask_coronal, image_size, 'nearest') > 0;
% hemisphere_mask_coronal = cat(3, hemisphere_mask_coronal, zeros(image_size(1), image_size(2), size(image_stack, 3) - size(hemisphere_mask_coronal, 3), 'like', hemisphere_mask_coronal));
% DataManager.visualize_itksnap(image_stack, annotation_moved_coronal);
%% Visualize annotation result with vessel image
% vis_itk_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack, 'whole_brain_d16x_registered'), ...
%     sprintf('%s_%s_d16x_registered', dataset_name, stack));
% DataManager.visualize_itksnap(image_stack_vsl, annotation_coronal, vis_itk_fp, true);
%% Generate registration information
allen_atlas_str = struct;
allen_atlas_str.registered_label_array = annotation_coronal;
allen_atlas_str.registered_hemisphere_mask = hemisphere_mask_coronal;
allen_atlas_str.structure_table = aa_structure_table;
allen_atlas_str.label_2_name = containers.Map(aa_structure_table.id, aa_structure_table.name);
allen_atlas_str.label_2_acronym = containers.Map(aa_structure_table.id, aa_structure_table.acronym);
% Convert the string of structure id path to number
allen_atlas_str.structure_id_path = allen_relabeled_data.structure_path;
allen_atlas_str.map_oldID_to_newID = allen_relabeled_data.map_old_id_to_new_id;
allen_atlas_str.map_newID_to_oldID = allen_relabeled_data.map_new_id_to_old_id;

allen_atlas_str.structure_id_depth = cellfun(@numel, allen_atlas_str.structure_id_path) - 1;
allen_atlas_str.annotation_filepath = registered_annotation_fp;
allen_atlas_str.info.fp = DataManager.fp_registration_data(dataset_name, stack, 'Allen_2017_25um_landmark.mat');
allen_atlas_str.info.voxel_size_um = [16, 16, 16];
allen_atlas_str.fun_name_to_ind = @(name_list, name) find(contains(name_list, name));
% Indices conversion 
allen_atlas_str.id = aa_structure_table.id;
allen_atlas_str.num_structure = numel(aa_structure_table.id);
allen_atlas_str.id_2_ind = sparse(allen_atlas_str.id, ones(size(allen_atlas_str.id)), ...
    1 : allen_atlas_str.num_structure, max(allen_atlas_str.id), 1);
allen_atlas_str.voxel_size = [16, 16, 16];
%% Compute the center of mass for each structure
allen_atlas_str.array_size = size(annotation_coronal);
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
itk_name = fullfile(DataManager.fp_itksnap_data_folder(dataset_name, stack, 'registration'), sprintf('%s_Allen_2017_25um_landmark_registration', stack));
DataManager.visualize_itksnap(image_stack, annotation_coronal, itk_name, true);
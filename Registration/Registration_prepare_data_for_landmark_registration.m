set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
fixed_stack = 'ML20190124';
moving_stack = 'ML20200201';

fixed_stack_im_fp = fullfile(DataManager.fp_processed_data(dataset_name, fixed_stack), ...
    'whole_stack_d16x.tiff');
moving_stack_im_fp = fullfile(DataManager.fp_processed_data(dataset_name, moving_stack), ...
    'whole_stack_d16x.tiff');
fixed_im = DataManager.load_single_tiff(fixed_stack_im_fp);
moving_stack_im = DataManager.load_single_tiff(moving_stack_im_fp);

save_folder_fp = DataManager.fp_mask_folder(dataset_name, ...
    fixed_stack, sprintf('Landmark_registration_from_%s', moving_stack));
if ~isfolder(save_folder_fp)
    mkdir(save_folder_fp);
end
%%
source_im_voxel_size = [16, 16, 16];
target_im_voxel_size = [25, 25, 25];

fixed_stack_target_size = round(size(fixed_im) .* source_im_voxel_size ./ target_im_voxel_size);
moving_stack_target_size = round(size(moving_stack_im) .* source_im_voxel_size ./ target_im_voxel_size);

fixed_im_rs = imresize3(fixed_im, fixed_stack_target_size);
moving_im_rs = imresize3(moving_stack_im, moving_stack_target_size);

% Permute the array
fixed_im_rs = permute(fixed_im_rs, [1, 3, 2]);
moving_im_rs = permute(moving_im_rs, [1, 3, 2]);
% Write to file
output_fp_fixed_im = fullfile(save_folder_fp, sprintf('%s_%s_d25x.tiff', ...
    dataset_name, fixed_stack));
DataManager.write_tiff_stack(fixed_im_rs, output_fp_fixed_im);
output_fp_moving_im = fullfile(save_folder_fp, sprintf('%s_%s_d25x.tiff', ...
    dataset_name, moving_stack));
DataManager.write_tiff_stack(moving_im_rs, output_fp_moving_im);
%% Load moving registration data
annotation_ext = 'nrrd';
hemisphere_ext = 'nrrd';
allen_annotation_registered_fp = fullfile(DataManager.fp_mask_folder(dataset_name, ...
    moving_stack, 'whole_brain_d16x_registration_landmark'), ...
    sprintf('allen_annotation_relabeled_landmark_registered.%s', annotation_ext));
allen_hemisphere_mask_registered_fp = fullfile(DataManager.fp_mask_folder(dataset_name, ...
    moving_stack, 'whole_brain_d16x_registration_landmark'), ...
    sprintf('allen_atlas_hemisphere_mask_landmark_registered.%s', hemisphere_ext));
assert(isfile(allen_annotation_registered_fp) && ...
    isfile(allen_hemisphere_mask_registered_fp), 'File does not exist');
annotation_target_fp = fullfile(save_folder_fp, sprintf('%s_%s_d25x_annotation_relabeled.%s',...
    dataset_name, moving_stack, annotation_ext));
hemisphere_mask_target_fp = fullfile(save_folder_fp, sprintf('%s_%s_d25x_hemisphere_mask.%s', ...
    dataset_name, moving_stack, hemisphere_ext));

copyfile(allen_annotation_registered_fp, annotation_target_fp);
copyfile(allen_hemisphere_mask_registered_fp, hemisphere_mask_target_fp);


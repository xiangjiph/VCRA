set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
reg_mask_version = 'whole_brain_d16x_registration';
stack_moving = 'ML20200201';
stack_fixed = 'ML_2018_08_15';

mask_moving_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack_moving, reg_mask_version), ...
    sprintf('%s_%s_d16x_registration_mask.nii.gz', dataset_name, stack_moving));
im_moving_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack_moving, reg_mask_version), ...
    sprintf('%s_%s_d16x_registration_image.nii', dataset_name, stack_moving));

mask_moving = DataManager.load_data(mask_moving_fp);
im_moving = DataManager.load_data(im_moving_fp);


mask_fixed_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack_fixed, reg_mask_version), ...
    sprintf('%s_%s_d16x_registration_mask.nii.gz', dataset_name, stack_fixed));
im_fixed_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack_fixed, reg_mask_version), ...
    sprintf('%s_%s_d16x_registration_image.nii', dataset_name, stack_fixed));

mask_fixed = DataManager.load_data(mask_fixed_fp);
im_fixed = DataManager.load_data(im_fixed_fp);
%% Remove the olfactory blub and spinal cord before registration
fixed_data_info = DataManager.load_dataset_info(dataset_name, stack_fixed);
moving_data_info = DataManager.load_dataset_info(dataset_name, stack_moving);
% implay(im_fixed)
% implay(im_moving)
im_fixed_valid_sec = 90 : 775;
im_moving_valid_sec = 202 : 856;
%% Remove image outside the mask 
im_moving_sb = im_moving .* uint16(mask_moving > 0);
im_fixed_sb = im_fixed .* uint16(mask_fixed > 0);

im_moving_sb = im_moving_sb(:, :, im_moving_valid_sec);
mask_moving_sb = mask_moving(:, :, im_moving_valid_sec);

im_fixed_sb = im_fixed_sb(:, :, im_fixed_valid_sec);
mask_fixed_sb = mask_fixed(:, :, im_fixed_valid_sec);

output_folder = DataManager.fp_mask_folder(dataset_name, stack_moving, ...
    'registration_to_ML20180815');
itk_fixed_fp = fullfile(output_folder, sprintf('%s_%s_d16x_sbc', dataset_name, ...
    stack_fixed));
DataManager.visualize_itksnap(im_fixed_sb, mask_fixed_sb, itk_fixed_fp, true);

itk_moving_fp = fullfile(output_folder, sprintf('%s_%s_d16x_sbc', dataset_name, ...
    stack_moving));
DataManager.visualize_itksnap(im_moving_sb, mask_moving_sb, itk_moving_fp, true);
%% Downsample the image
ds_ratio = 0.25;
im_mov_ds = imresize3(im_moving_sb, ds_ratio);
im_fixed_ds = imresize3(im_fixed_sb, ds_ratio);
%% Affine registration
[optimizer, metric] = imregconfig('multimodal');
optimizer.MaximumIterations = 500;
affine_reg_info = imregtform(im_mov_ds, im_fixed_ds, 'affine', optimizer, metric, 'DisplayOptimization', true);
% im_moved = imwarp(im_mov_ds, affine_reg_info, 'OutputView', imref3d(size(im_fixed_sb)));
% annotation_moved = imwarp(atlas_annotation, affine_reg_info, 'nearest', 'OutputView', imref3d(size(im_fix)));
%%
[matQ, matR] = qr(affine_reg_info.T(1:3, 1:3));
%% Analyze 3D slicer data
registration_data_fp_cell = {'Affine_WholeBrain_ML20200201_d16x_registration_image_To_WholeBrain_ML_2018_08_15_d16x_registration_image.mat', ...
    'Affine_WholeBrain_ML20200201_d16x_registration_image_To_WholeBrain_ML_2018_08_15_d16x_registration_image_ini_mask_affine.mat', ...
    'Scale_WholeBrain_ML20200201_d16x_registration_image_To_WholeBrain_ML_2018_08_15_d16x_registration_image.mat', ...
    'Affine_WholeBrain_ML20200201_d16x_sbc_image_To_WholeBrain_ML_2018_08_15_d16x_sbc_image_ini_mask_affine.mat', ...
    'Affine_WholeBrain_ML20200201_d16x_sbc_mask_To_WholeBrain_ML_2018_08_15_d16x_sbc_mask.mat', ...
    'Scale_WholeBrain_ML20200201_d16x_sbc_image_To_WholeBrain_ML_2018_08_15_d16x_sbc_image.mat'};
num_reg_data = numel(registration_data_fp_cell);
for iter_cell = 1 : num_reg_data
    tmp_fp = fullfile(DataManager.fp_mask_folder(dataset_name, stack_moving, 'registration_to_ML20180815'), ...
        registration_data_fp_cell{iter_cell});
    tmp_data = load(tmp_fp);
    tmp_data = reshape(tmp_data.AffineTransform_double_3_3(1:9), 3, 3);
    tmp_data = inv(tmp_data);
    [~, tmp_scale_mat] = qr(tmp_data);
    fprintf('%s scaling matrix\n', registration_data_fp_cell{iter_cell});
    disp(tmp_scale_mat);
    fprintf('\tDeterminant: %f\n', det(tmp_scale_mat));
    fprintf('\tAverage linear deformation: %f\n', det(tmp_scale_mat)^(1/3));
end
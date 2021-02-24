function wb_data_str = fun_analysis_load_whole_brain_data_for_regional_analysis(...
    dataset_name, stack, grid_version, reconstruction_version, skeleton_version, ...
    registration_version, wb_data_folder_name, load_skel_Q)
if nargin < 8
    load_skel_Q = true;
end
% To do list: 
% 1. Remove the hard code in this file. Maybe put the whole brian
% mask for distance transform calculation into the registration structure. 
persistent DataManager
if isempty(DataManager)
    DataManager = FileManager;
end
start_tic = tic;
wb_data_str = struct;
wb_data_str.dataset_name = dataset_name;
wb_data_str.stack = stack;
wb_data_str.grid_version = grid_version;
wb_data_str.skeleton_version = skeleton_version;
wb_data_str.reconstruction_version = reconstruction_version;
wb_data_str.registration_version = registration_version;
fprintf('Loading grid information...\n');
wb_data_str.grid_info = DataManager.load_grid(wb_data_str.dataset_name,...
    wb_data_str.stack, wb_data_str.grid_version);
fprintf('Loading registration data...\n');
% Registration data
registration_str = DataManager.load_registration_data(wb_data_str.dataset_name,...
    wb_data_str.stack, wb_data_str.registration_version);
if ~isfield(registration_str.info, 'voxel_size_um')
    warning('Missing field');
    registration_str.info.voxel_size_um = [16, 16, 16];
end
wb_data_str.registration = registration_str;
% Brain mask
fprintf('Loading whole brain mask and compute the distance transform...\n');
wb_mask = DataManager.load_data(sprintf('%s_mask.nii.gz', ...
    fullfile(DataManager.fp_mask_folder(wb_data_str.dataset_name, wb_data_str.stack, 'whole_brain_d16x_registration'),...
    sprintf('%s_%s_d16x_registration', wb_data_str.dataset_name, wb_data_str.stack)))) > 0;
wb_mask_str = struct;
wb_mask_str.ds_ratio = 16;
wb_mask_str.dt = bwdist(~wb_mask) .* wb_mask_str.ds_ratio;
wb_data_str.brain_mask = wb_mask_str;
% Node and link features in each 240 cube and skeleton. Due to the size
% of the data, it seems faster to load the data using parfor directly
% from the hard drive, instead of saving a large file
fprintf('Loading skeleton data, node and link features in all the 240 cubes...\n');
% The most time-consuming step. Not necessary for all the analysis. 
if load_skel_Q
    wb_data_str.cube_data = fun_analysis_load_whole_brain_240_cube_features_N_skel(...
        wb_data_str.grid_info, wb_data_str.reconstruction_version, wb_data_str.skeleton_version);
end

fprintf('Loading 240 cube statistics...\n');
wb_data_str.cube_stat = DataManager.load_analysis_data(wb_data_str.dataset_name, wb_data_str.stack, ...
    sprintf('%s_%s_%s_240_cube_stat_data.mat', wb_data_str.dataset_name, wb_data_str.stack, wb_data_str.reconstruction_version), wb_data_folder_name);
fprintf('Finish loading all the data. Elapsed time is %f seconds.\n', toc(start_tic));
end
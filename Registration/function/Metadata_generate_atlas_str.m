dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
registration_name = 'Allen_2017_25um_nonrigid.mat';
registration_str = DataManager.load_registration_data(dataset_name, stack, registration_name);
% Remove fields
registration_str = rmfield(registration_str, {'registered_label_array', 'registered_hemisphere_mask', 'array_size', ...
    'structure_sub_left', 'structure_sub_right', 'structure_cc_ind', 'structure_cc_ind_hemisphere'});
save('./Metadata/Allen_atlas.mat', '-struct', 'registration_str');
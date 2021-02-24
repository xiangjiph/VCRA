function registration_vis_str = fun_vis_permute_registration_str(registration_vis_str, vis_dir)

permute_field_name = {'registered_label_array', 'registered_hemisphere_mask'...
    'combined_regional_mask'};
permute_field_name = intersect(permute_field_name, fieldnames(registration_vis_str));
for iter_permute = 1 : numel(permute_field_name)
    tmp_data = registration_vis_str.(permute_field_name{iter_permute});
    switch vis_dir
        case 1
            tmp_data = permute(tmp_data, [2, 3, 1]);
        case 2
            tmp_data = permute(tmp_data, [1, 3, 2]);
        case 3
        otherwise 
            error('Visualization axis should be any one of {1, 2, 3}');
    end
    registration_vis_str.(permute_field_name{iter_permute}) = tmp_data;
end
end
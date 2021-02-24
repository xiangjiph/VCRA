function wb_cube_stat_str = fun_analysis_collect_regional_cube_stat_by_list_ind(wb_cube_stat_str, list_ind)

direct_mask_field = {'cap2vsl_length_fraction', 'cap2vsl_surf_area_fraction', ...
    'cap2vsl_vol_fraction', 'capillary_length_density_m_mm3', ...
    'capillary_surface_area_density_mm2_mm3', 'capillary_volume_density', ...
    'cube_in_brain_mask_ratio', 'link_length_density_m_mm3', ...
    'link_surface_area_density_mm2_mm3', 'link_volume_density', ...
    'mask_surface_area_density_mm2_mm3', 'mask_volume_density', ...
    'mask_volume', 'num_surface_voxel', 'num_link', 'num_cap'};
direct_mask_field = intersect(direct_mask_field, fieldnames(wb_cube_stat_str));
for iter_field = 1 : numel(direct_mask_field)
    wb_cube_stat_str.(direct_mask_field{iter_field}) = wb_cube_stat_str.(direct_mask_field{iter_field})(list_ind);
end
%% Copy 240 cube network statistics
cube_stat_field = {'link_all_stat', 'link_cap_stat', 'node_stat'};
for iter_field = 1 : numel(cube_stat_field)
    if ~isfield(wb_cube_stat_str, cube_stat_field{iter_field})
        continue;
    end
    tmp_data_str = wb_cube_stat_str.(cube_stat_field{iter_field});
    tmp_data_field = fieldnames(tmp_data_str);
    for iter_field_1 = 1 : numel(tmp_data_field)
        tmp_stat_data = tmp_data_str.(tmp_data_field{iter_field_1});
        tmp_stat_feature_name = fieldnames(tmp_stat_data);
        for iter_feature = 1 : numel(tmp_stat_feature_name)
            tmp_stat_data.(tmp_stat_feature_name{iter_feature}) = tmp_stat_data.(tmp_stat_feature_name{iter_feature})(list_ind);
        end
        tmp_data_str.(tmp_data_field{iter_field_1}) = tmp_stat_data;
    end
    wb_cube_stat_str.(cube_stat_field{iter_field}) = tmp_data_str;
end
%% Copy anisotropy data
ai_data_field = {'wb_ai_all_lw', 'wb_ai_all_vw', 'wb_ai_cap_lw', 'wb_ai_cap_vw'};
ai_data_field = intersect(ai_data_field, fieldnames(wb_cube_stat_str));
for iter_field = 1 : numel(ai_data_field)
    if ~isfield(wb_cube_stat_str, ai_data_field{iter_field})
        continue;
    end
    tmp_data_str = wb_cube_stat_str.(ai_data_field{iter_field});
    tmp_data_field = fieldnames(tmp_data_str);
    for iter_feature = 1 : numel(tmp_data_field)
        if isempty(tmp_data_str.(tmp_data_field{iter_feature}))
            tmp_data_str.(tmp_data_field{iter_feature}) = [];
            continue;
        end
        tmp_data_str.(tmp_data_field{iter_feature}) = tmp_data_str.(tmp_data_field{iter_feature})(list_ind, :);
    end
    wb_cube_stat_str.(ai_data_field{iter_field}) = tmp_data_str;
end
if isfield(wb_cube_stat_str, 'grid_info')
    wb_cube_stat_str = rmfield(wb_cube_stat_str, 'grid_info');
end
end
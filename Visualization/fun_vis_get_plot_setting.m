function plot_setting = fun_vis_get_plot_setting(compare_feature_subgroup)

switch compare_feature_subgroup
    case {'link_all'}
        plot_setting = cell(4, 0);
        plot_setting(:, end+1) = {'length', 'Vessel segment length(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'dt_median', 'Vessel radius(\mum)', 'linear', 'log'};
        plot_setting(:, end+1) = {'surface_area', 'Vessel segment surface area(\mum^2)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'volume', 'Vessel segment volume(\mum^3)', 'linear', 'linear'};
        
        plot_setting(:, end+1) = {'nearest_tissue_radius', 'Nearest tissue radius(\mum)', 'linear', 'linear'};
        
        plot_setting(:, end+1) = {'capillary_branching_order', 'Capillary branching order', 'linear', 'linear'};
        plot_setting(:, end+1) = {'shortest_loop_length', 'Length of shortest loop(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'shortest_loop_geodesic_length', 'Number of edges', 'linear', 'linear'};
        
        plot_setting(:, end+1) = {'dist_to_nearest_noncapillary_mean', 'Distance to the nearest noncapillary(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'straightness', 'Straightness', 'linear', 'linear'};
    case {'link_cap'}
        plot_setting = cell(4, 0);
        plot_setting(:, end+1) = {'length', 'Capillary segment length(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'dt_median', 'Capillary radius(\mum)', 'linear', 'log'};
        plot_setting(:, end+1) = {'surface_area', 'Capillary segment surface area(\mum^2)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'volume', 'Capillary segment volume(\mum^3)', 'linear', 'linear'};
        
        plot_setting(:, end+1) = {'nearest_tissue_radius', 'Nearest tissue radius(\mum)', 'linear', 'linear'};
        
        plot_setting(:, end+1) = {'capillary_branching_order', 'Capillary branching order', 'linear', 'linear'};
        plot_setting(:, end+1) = {'shortest_loop_length', 'Length of the shortest loop(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'shortest_loop_geodesic_length', 'Number of edges', 'linear', 'linear'};
        
        plot_setting(:, end+1) = {'dist_to_nearest_noncapillary_mean', 'Distance to the nearest noncapillary(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'straightness', 'Straightness', 'linear', 'linear'};
    case {'node'}
        plot_setting = cell(4, 0);
        plot_setting(:, end+1) = {'nearest_node_dist', 'Distance to the nearest node(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'path_to_nearest_neighbor_geodesic_length', 'Number of edges to the nearest node', 'linear', 'log'};
        plot_setting(:, end+1) = {'degree', 'degree', 'linear', 'linear'};
        plot_setting(:, end+1) = {'link_length_min', 'Length of the shortest link(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'link_length_max', 'Length of the longest link(\mum)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'link_length_median', 'Median length of the links(\mum)', 'linear', 'linear'};
    case 'network_prop'
        plot_setting = cell(4, 0);
        plot_setting(:, end+1) = {'surface_area_mm2mm3', 'Vessel surface area density(mm^2/mm^3)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'volume_ratio', 'Vessel volume ratio', 'linear', 'linear'};
        plot_setting(:, end+1) = {'length_density_m_mm3', 'Vessel length density(m/mm^3)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'length_density_cap_m_mm3', 'Capillary length density(m/mm^3)', 'linear', 'linear'};
        plot_setting(:, end+1) = {'ai_all_vw_fa', 'Vessel fractional anisotropy', 'linear', 'linear'};
        plot_setting(:, end+1) = {'ai_all_vw_min2max_z', 'Vessel isotropy score', 'linear', 'linear'};
        plot_setting(:, end+1) = {'ai_all_vw_svd1_z', 'Vessel anisotropy score', 'linear', 'linear'};
        plot_setting(:, end+1) = {'ai_cap_vw_min2max_z', 'Capillary isotropy score', 'linear', 'linear'};
        plot_setting(:, end+1) = {'ai_cap_vw_svd1_z', 'Capillary anisotropy score', 'linear', 'linear'};
        plot_setting(:, end+1) = {'ai_cap_vw_fa', 'Capillary fractional anisotropy', 'linear', 'linear'};
    otherwise
        fprintf('Unrecognized feature group. Do not add plot settings\n');
end
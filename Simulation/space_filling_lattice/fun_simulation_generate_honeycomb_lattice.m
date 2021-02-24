function honeycomb_lattice_str = fun_simulation_generate_honeycomb_lattice(bond_length_um, system_length_um, pixel_size_um)

% Initialization 
system_length_pxl = round(system_length_um / pixel_size_um);
bond_length_pxl = round(bond_length_um / pixel_size_um);

est_num_site = round(system_length_pxl / (bond_length_pxl));
%% Put a site at the origin 
point_sub_cell = cell(1, 4);

[tmp_sub_1, tmp_sub_2] = ndgrid(round(3 * bond_length_pxl * (0 : 1 : est_num_site)) , ...
    round(sqrt(3) * bond_length_pxl * (0 : 1 : est_num_site)));
tmp_sub = cat(2, tmp_sub_1(:), tmp_sub_2(:));
point_sub_cell{1} = tmp_sub;
tmp_sub(:, 1) = tmp_sub(:, 1) - bond_length_pxl;
point_sub_cell{2} = tmp_sub;

[tmp_sub_1, tmp_sub_2] = ndgrid(round((0.5 + 3 * (0 : 1 : est_num_site)) * bond_length_pxl), ...
    round((sqrt(3)/2 + sqrt(3) * (0 : 1 : est_num_site)) * bond_length_pxl));
tmp_sub = cat(2, tmp_sub_1(:), tmp_sub_2(:));
point_sub_cell{3} = tmp_sub;
tmp_sub(:, 1) = tmp_sub(:, 1) + bond_length_pxl;
point_sub_cell{4} = tmp_sub;
%%
% Add site
system_mask = false(system_length_pxl, system_length_pxl);
system_size = size(system_mask);

num_site_per_dim = est_num_site + 1;
point_valid_Q_cell = cell(1, 4);
for iter_point_set = 1 : 4
    tmp_point_sub = point_sub_cell{iter_point_set};
    tmp_is_valid_Q = all(tmp_point_sub >= 1 & tmp_point_sub <= system_length_pxl, 2);
    tmp_point_ind = sub2ind(system_size, tmp_point_sub(tmp_is_valid_Q, 1), tmp_point_sub(tmp_is_valid_Q, 2));
    system_mask(tmp_point_ind) = true;
    tmp_is_valid_Q = reshape(tmp_is_valid_Q, num_site_per_dim, num_site_per_dim);
    point_valid_Q_cell{iter_point_set} = tmp_is_valid_Q;
end
% Add bond
[valid_site_1, valid_site_2, valid_site_3, valid_site_4] = point_valid_Q_cell{:};
[site_1_sub, site_2_sub, site_3_sub, site_4_sub] = point_sub_cell{:};
for iter_1 = 1 : num_site_per_dim
    for iter_2 = 1 : num_site_per_dim
        tmp_ind = sub2ind([num_site_per_dim, num_site_per_dim], iter_1, iter_2);
        if valid_site_1(tmp_ind) && valid_site_2(tmp_ind)
            tmp_2p2_ind = fun_graph_get_p2p_line_mask(system_size, ...
                site_1_sub(tmp_ind, :), site_2_sub(tmp_ind, :));
            system_mask(tmp_2p2_ind) = true;
        end
        
        if valid_site_3(tmp_ind) && valid_site_4(tmp_ind)
            tmp_2p2_ind = fun_graph_get_p2p_line_mask(system_size, ...
                site_3_sub(tmp_ind, :), site_4_sub(tmp_ind, :));
            system_mask(tmp_2p2_ind) = true;
        end
        
        if valid_site_1(tmp_ind) && valid_site_3(tmp_ind)
            tmp_2p2_ind = fun_graph_get_p2p_line_mask(system_size, ...
                site_1_sub(tmp_ind, :), site_3_sub(tmp_ind, :));
            system_mask(tmp_2p2_ind) = true;
        end
        
        if iter_2 > 1
            tmp_ind_2 = sub2ind([num_site_per_dim, num_site_per_dim], iter_1, iter_2 - 1);
            if valid_site_1(tmp_ind) && valid_site_3(tmp_ind_2)
                tmp_2p2_ind = fun_graph_get_p2p_line_mask(system_size, ...
                    site_1_sub(tmp_ind, :), site_3_sub(tmp_ind_2, :));
                system_mask(tmp_2p2_ind) = true;
            end
        end
        
        if iter_1 < num_site_per_dim
            tmp_ind_2 = sub2ind([num_site_per_dim, num_site_per_dim], iter_1 + 1, iter_2);
            if valid_site_2(tmp_ind_2) && valid_site_4(tmp_ind)
                tmp_2p2_ind = fun_graph_get_p2p_line_mask(system_size, ...
                    site_2_sub(tmp_ind_2, :), site_4_sub(tmp_ind, :));
                system_mask(tmp_2p2_ind) = true;
            end
        end
        
        if iter_1 < num_site_per_dim && iter_2 > 1
            tmp_ind_1 = sub2ind([num_site_per_dim, num_site_per_dim], iter_1 + 1, iter_2);
            tmp_ind_2 = sub2ind([num_site_per_dim, num_site_per_dim], iter_1, iter_2 - 1);
            if valid_site_2(tmp_ind_1) && valid_site_4(tmp_ind_2)
                tmp_2p2_ind = fun_graph_get_p2p_line_mask(system_size, ...
                    site_2_sub(tmp_ind_1, :), site_4_sub(tmp_ind_2, :));
                system_mask(tmp_2p2_ind) = true;
            end
        end
    end
end
mask_skl = bwskel(system_mask);
system_graph = fun_skeleton_to_graph_2D(mask_skl);
%%
honeycomb_lattice_str = struct;
honeycomb_lattice_str.mask = mask_skl;
honeycomb_lattice_str.system_graph = system_graph;
honeycomb_lattice_str.connectivity_graph_str = fun_analysis_get_connectivity_graph(system_graph);
honeycomb_lattice_str.unweighted_graph = honeycomb_lattice_str.connectivity_graph_str.graph_uw;
end
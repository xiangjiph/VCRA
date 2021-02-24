function subregion_list_ind = fun_registration_get_subregion_list_ind(registration_str, region_id, verboseQ)

if nargin < 3
    verboseQ = true;
end

if isscalar(region_id)
    is_subregion_Q = any(registration_str.structure_id_path_mat == region_id, 2);
else
    is_subregion_Q = false(size(registration_str.structure_id_path_mat, 1), 1);
    for iter_id = 1 : numel(region_id)
        is_subregion_Q = is_subregion_Q | any(registration_str.structure_id_path_mat == region_id(iter_id), 2);
    end
end
% Display the information for the subregion
if verboseQ
    fprintf('The following are the substructures:\n');
    disp(registration_str.structure_table(is_subregion_Q, 1:5));
end
subregion_id = registration_str.id(is_subregion_Q);
% Get the array voxel indices for the subregion connected components
subregion_list_ind = full(registration_str.id_2_ind(subregion_id));
subregion_list_ind = subregion_list_ind(subregion_list_ind > 0);
end
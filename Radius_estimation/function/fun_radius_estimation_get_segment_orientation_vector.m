function ori_vec = fun_radius_estimation_get_segment_orientation_vector(voxel_sub, method)

if nargin < 2
    method = 'auto';
end
[num_vxl, vxl_dim] = size(voxel_sub);
assert(vxl_dim == 3);

switch method
    case 'auto'
        if num_vxl >= 3
            % Use endpoint to endpoint
            ori_vec = fun_analysis_get_ordered_vxl_sub_list_ori_vec(voxel_sub, 'svd');
        elseif num_vxl == 2
            ori_vec = fun_analysis_get_ordered_vxl_sub_list_ori_vec(voxel_sub, 'ep2ep');
        elseif num_vxl == 1
            ori_vec = [sqrt(3/8), sqrt(3/8), 1/2]; %Make the z-component 0 
        else
            error('Number of voxels in the link connected component is less than 1');
        end
    otherwise
        ori_vec = fun_analysis_get_ordered_vxl_sub_list_ori_vec(voxel_sub, method);
end         
end
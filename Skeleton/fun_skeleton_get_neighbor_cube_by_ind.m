function [nhood, varargout]= fun_skeleton_get_neighbor_cube_by_ind(input_array,ind_list, ind_offset_list)
% fun_vectorize_get_neighbor_cube get the 3x3x3 26-neighbor cube of the
% points in ind_listfrom the input_array. 
% Input:
%   input_array: 3d array
%   ind_list: N-by-1 array of the indices of voxel in the input_array
%   ind_offset_list: output by fun_skeleton_neighbor_add_coeff_3D
% Output:
%   nhood: N-by 27 array extracted from the neighboring 3x3x3 cube of the
%   specified voxels in the input_array.
% Note:
%   This function does not take care of the boundary. All of the voxels on the
%   face of the input array should be 0.
% Compute the indice location of the neighbors
if isscalar(ind_list)
    nhood_ind = ind_list + ind_offset_list';
else
    nhood_ind = bsxfun(@plus, ind_list, ind_offset_list');
end
% Extract neighboring voxel values 
nhood = input_array(nhood_ind);
if nargout > 1
    varargout = nhood_ind;
end

end

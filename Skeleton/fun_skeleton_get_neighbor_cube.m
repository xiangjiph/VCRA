function [nhood, varargout]= fun_skeleton_get_neighbor_cube(input_array,ind_list, connectivity, delete_center_Q)
% fun_vectorize_get_neighbor_cube get the 3x3x3 26-neighbor cube of the
% points in ind_listfrom the input_array. 
% Input:
%   input_array: 3d array
%   ind_list: N-by-1 array of the indices of voxel in the input_array
%   delete_center_Q: logical scalar. If true, delete the center point and output N-by-26 array
% Output:
%   nhood: N-by 27 array extracted from the neighboring 3x3x3 cube of the
%   specified voxels in the input_array.
% Note:
%   This function does not take care of the boundary. All of the voxels on the
%   face of the input array should be 0.
% Compute the indice location of the neighbors
persistent ind_add r_array_size r_connectivity r_delete_center_Q
array_size = size(input_array);
if nargin < 3
    connectivity = 26;
    delete_center_Q = false;
elseif nargin < 4
    delete_center_Q = false;
end

% Define persistent variables to save some time on the ind_add computation.
if isempty(r_array_size)
    r_array_size = array_size;
end

if isempty(r_connectivity)
    r_connectivity = connectivity;
end
if isempty(r_delete_center_Q)
    r_delete_center_Q = delete_center_Q;
end

if isempty(ind_add) || any(r_array_size ~= array_size) || ...
        (r_delete_center_Q ~= delete_center_Q)|| (r_connectivity ~= connectivity)
    r_array_size = array_size;
    r_connectivity = connectivity;
    r_delete_center_Q = delete_center_Q;
    ind_add = fun_skeleton_neighbor_add_coeff_3D(r_array_size, r_connectivity, r_delete_center_Q);
end
nhood_ind = bsxfun(@plus, ind_list, ind_add');
% Extract neighboring voxel values 
nhood = input_array(nhood_ind);
if nargout > 1
    varargout = nhood_ind;
end

end

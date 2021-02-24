function neighbor_add_coeff = fun_skeleton_neighbor_add_coeff_3D(block_size, connectivity, remove_center_Q)
% fun_skeleton_neighbor_add_coeff_3D compute the ind difference between the
% neighbors and the center point in a 3x3x3 cube
% Input: 
%   block_size: 3-by-1 numerical array
%   connectivity: numerical scaler(6, 18 or 26),  type of neighboring
%   voxels. Default value is 26. 
%   removal_center_Q: logical scalar, if ture, remove center line index (0)
%   from the output array. Default value is false. 
if nargin < 2
    connectivity = 26;
    remove_center_Q = false;
elseif nargin < 3
    remove_center_Q = false;
end

switch connectivity
    case 26 
        [x,y,z] = ndgrid(1:3, 1:3, 1:3);
    case 6
        cube = zeros(3,3,3);
        cube(2,2,:) = 1;
        cube(:,2,2) = 1;
        cube(2,:,2) = 1;
        [x,y,z] = ind2sub([3,3,3], find(cube(:)));
    case 18
        cube = zeros(3,3,3);
        cube(2,:,:) = 1;
        cube(:,2,:) = 1;
        cube(:,:,2) = 1;
        [x,y,z] = ind2sub([3,3,3], find(cube(:)));
    otherwise 
        error('Connectivity should be 6, 18 or 26');
end
middle_idx = connectivity/2 + 1;
neighbor_add_coeff = sub2ind(block_size, x, y, z);
neighbor_add_coeff = neighbor_add_coeff - neighbor_add_coeff (middle_idx);

if remove_center_Q
    neighbor_add_coeff(middle_idx) = [];
end
neighbor_add_coeff = neighbor_add_coeff(:);

end
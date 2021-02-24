function cc_length = fun_graph_sub_to_length(sub_list, voxel_size)
% Compute the length of the link segments. Assume the sub_list is the
% subscript of the voxel in a single connected components and the
% subscript list should be in order. In other words, the absolute neighboring
% subscript difference in each direction shoudl be no larger than 1.
% Input: 
%   sub_list: N-by-D numerical array, list of N voxel subscripts in D
%   dimension
%   voxel_size: scalar or D-by-1 numerical arraysize of the voxle
% Output: 
%   cc_length: the length of the connected component
%
[num_voxel, voxel_dim] = size(sub_list);
if nargin < 2
    voxel_size = 1;
end
if isscalar(voxel_size)
    voxel_size = voxel_size .* ones(voxel_dim,1);
end
if num_voxel == 1
    cc_length = voxel_size(1);
else
    cc_length  = sum(sqrt(sum(bsxfun(@times, diff(sub_list,1)', voxel_size).^2 ,1))) + voxel_size(1);   
end
end
function [voxel_sub, varargout] = fun_uniform_sample_points_in_space(voxel_sub, sample_grid_size, num_sample_per_block, selection_rule)
% fun_uniform_sample_points_in_space uniformly sample voxels in the space
% that the voxels live in. 
% Input: 
%   voxel_sub: N-by-D numerical array, subscript of the voxels, where N is
%   the number of voxels and D is the dimension of the voxels. 
%   sample_grid_size: D-by-1 numerical array. Grid size for the uniform
%   sampling. 
%   num_smaple_per_block: scalar or numerical vector, number of samples per sampling grid 
%   selection_rule: string
%       'random': randomly select 1 voxel in the grid 
%       'first': select the first one
% Output: 
%   voxel_sub: sampled voxel subscript list
%   sampled_ind: indices of the voxels in the input voxel list. 
% 

if nargin < 4
    selection_rule = 'random';
end
if isempty(voxel_sub)
    varargout{1} =  [];
    return;
end

bbox_1_min = min(voxel_sub, [], 1);
bbox_1_max = max(voxel_sub, [], 1);
grid_sub = ceil(bsxfun(@minus, voxel_sub, bbox_1_min - 1) ./ sample_grid_size);
grid_size = ceil((bbox_1_max - bbox_1_min + 1) ./ sample_grid_size);
grid_ind = sub2ind(grid_size, grid_sub(:,1), grid_sub(:,2), grid_sub(:,3));
block_ind = fun_bin_data_to_idx_list(grid_ind);
num_block = numel(block_ind);

sampled_ind = cell(num_block, 1);
for iter_block = 1 : num_block
    tmp_ind = block_ind{iter_block};
    if numel(tmp_ind) <= num_sample_per_block
        sampled_ind{iter_block} = tmp_ind;
    else
        switch selection_rule
            case 'random'
                sampled_ind{iter_block} = randsample(tmp_ind, num_sample_per_block);
            case 'first'            
                sampled_ind{iter_block} = tmp_ind(1 : num_sample_per_block);
            otherwise
                error('Unrecognized selection rule');
        end
    end
end
sampled_ind = cat(2, sampled_ind{:});
voxel_sub = voxel_sub(sampled_ind, :);
if nargout > 1
    varargout{1} = sampled_ind;
end
end
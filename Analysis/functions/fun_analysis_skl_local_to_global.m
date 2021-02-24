function skl = fun_analysis_skl_local_to_global(skl)
% fun_analysis_skl_local_to_global converts the skeleton indices in local
% block into the indices in the entire dataset and remove the backup. 
%
if isfield(skl, 'backup')
    skl = rmfield(skl, 'backup');
end
num_voxel = numel(skl.ind);
local_sub = fun_ind2sub(skl.block_size, skl.ind);
global_sub = local_sub + skl.global_bbox_mmxx(1:3) - 1;
if num_voxel > 1
    skl.ind = sub2ind(skl.global_block_size, global_sub(:,1), global_sub(:,2), global_sub(:,3));
elseif num_voxel == 1
    skl.ind = sub2ind(skl.global_block_size, global_sub(1), global_sub(2), global_sub(3));
end
end
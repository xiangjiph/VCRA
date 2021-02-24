function varargout = fun_vis_voxel_list_local(voxel_ind_list, global_block_size)

[sub_1, sub_2, sub_3] = ind2sub(global_block_size, voxel_ind_list);
min_1 = min(sub_1);
min_2 = min(sub_2);
min_3 = min(sub_3);
max_1 = max(sub_1);
max_2 = max(sub_2);
max_3 = max(sub_3);
local_mask_size = max([max_1 - min_1 + 1, max_2 - min_2 + 1, max_3 - min_3 + 1], 2);
sub_1 = sub_1 - min_1 + 1;
sub_2 = sub_2 - min_2 + 1;
sub_3 = sub_3 - min_3 + 1;
voxel_ind_list = sub2ind(local_mask_size, sub_1, sub_2, sub_3);
mask = false(local_mask_size);
mask(voxel_ind_list) = true;
if nargout == 1
    varargout = mask;
else
    volumeViewer(mask)
end

end


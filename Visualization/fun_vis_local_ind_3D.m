function fun_vis_local_ind_3D(voxel_ind, array_size);

vis_array = false(array_size);
vis_array(voxel_ind) = true;
volumeViewer(vis_array);
end

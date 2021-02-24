function r_vxl_sub= fun_radius_estimation_voxel_local_sub_ds2r(ds_voxel_sub, ds_local_bbox_mm, ds_vxl_size, r_vxl_size, r_local_bbox_mm)

assert(size(ds_voxel_sub, 2) == 3);
assert(isrow(ds_local_bbox_mm) && isrow(r_local_bbox_mm));
if isempty(ds_voxel_sub)
    r_vxl_sub = ds_voxel_sub;
else
    ds_vxl_sub_global = ds_voxel_sub + ds_local_bbox_mm - 1;
    ds_vxl_sub_global_um = ds_vxl_sub_global .* ds_vxl_size;
    r_vxl_sub_global = max(1, round(ds_vxl_sub_global_um ./ r_vxl_size));
    r_vxl_sub = r_vxl_sub_global - r_local_bbox_mm + 1;
end
end
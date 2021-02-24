function str = fun_graph_get_link_feature_bg(link_sub, vessel_image, vessel_mask, local_bbox_expansion)


image_size = size(vessel_image);
tmp_min = max([1,1,1],  min(link_sub - local_bbox_expansion, [], 1));
tmp_max = min(image_size, max(link_sub + local_bbox_expansion, [], 1));
tmp_mask = vessel_mask(tmp_min(1):tmp_max(1), tmp_min(2):tmp_max(2), tmp_min(3):tmp_max(3));
tmp_mask = tmp_mask(:);
tmp_bg_mask = ~tmp_mask;
tmp_image = vessel_image(tmp_min(1):tmp_max(1), tmp_min(2):tmp_max(2), tmp_min(3):tmp_max(3));
tmp_image = single(tmp_image(:));
tmp_fg_voxel = nnz(tmp_mask);
tmp_bg_voxel = nnz(tmp_bg_mask);

str.bg_mean = sum(tmp_image(:) .* tmp_bg_mask(:)) /tmp_bg_voxel;
tmp_var = (sum(tmp_image(:).^2 .* tmp_bg_mask(:)) /tmp_bg_voxel) - str.bg_mean^2;
if tmp_var < 0
    tmp_var = var(tmp_image);
    if tmp_var < 0
        error('Negative variance');
    end
end
str.bg_std = sqrt(tmp_var);
str.mask_mean = sum(tmp_image .* tmp_mask) /tmp_fg_voxel;
str.mask_std = sqrt(sum(tmp_image.^2 .* tmp_mask) /tmp_fg_voxel - str.bg_mean^2);

end
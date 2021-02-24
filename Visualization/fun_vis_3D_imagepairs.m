function fig_hdl = fun_vis_3D_imagepairs(im1, im2, vis_prctile_list)

assert(all(vis_prctile_list <= 1));
num_vis_dir = 3;
num_vis_per_dir = numel(vis_prctile_list);
fig_hdl = figure;
t_hdl = tiledlayout(num_vis_dir, num_vis_per_dir);
for iter_dir = 1 : num_vis_dir
    
    switch iter_dir
        case 3 
            permute_order = [1, 2, 3];        
        case 2
            permute_order = [1, 3, 2];
        case 1
            permute_order = [3, 2, 1];
    end    
    tmp_vis_im_1 = permute(im1, permute_order);
    tmp_vis_im_2 = permute(im2, permute_order);
    tmp_vis_range = size(tmp_vis_im_1, 3);
    tmp_vis_sec_list = max(1, unique(round(tmp_vis_range * vis_prctile_list), 'stable'));
    for iter_plane = 1 : num_vis_per_dir
        tmp_ax_hdl = nexttile;        
        imshowpair(tmp_vis_im_1(:, :, tmp_vis_sec_list(iter_plane)),...
            tmp_vis_im_2(:, :, tmp_vis_sec_list(iter_plane)));
    end    
end

end
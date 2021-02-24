function link_global_bbox = fun_vis_link_surrounding_by_cc_ind(vessel_image, vessel_mask, vessel_skl, link_ind, link_bbox_expansion_length)

if nargin < 5
    link_bbox_expansion_length = 10;
end

fig_hdl = figure('Name', 'Link connected component visualization', 'Visible', 'off', 'Units', 'normalized');
mask_size = size(vessel_mask);
link_global_sub = zeros(numel(link_ind), 3);
[link_global_sub(:,1), link_global_sub(:,2), link_global_sub(:,3)] = ind2sub(mask_size, link_ind);

link_global_bbox = zeros(1,6);
link_global_bbox(1:3) = max(1, min(link_global_sub, [], 1) - link_bbox_expansion_length);
link_global_bbox_max = min(mask_size, max(link_global_sub, [], 1) + link_bbox_expansion_length) ;
link_global_bbox(4:6) = link_global_bbox_max - link_global_bbox(1:3) + 1;
% Select the link sub witin the block
link_local_sub = bsxfun(@minus, link_global_sub, link_global_bbox(1:3)) + 1;

vis_local_skl = crop_bbox3(vessel_skl, link_global_bbox, 'default');
vis_skl_ind = find(vis_local_skl);
vis_skl_ind = setdiff(vis_skl_ind, sub2ind(link_global_bbox(4:6), link_local_sub(:,1), link_local_sub(:,2), link_local_sub(:,3) ));
vis_skl_sub = zeros(numel(vis_skl_ind), 3);
[vis_skl_sub(:,1), vis_skl_sub(:,2), vis_skl_sub(:,3)] = ind2sub(link_global_bbox(4:6), vis_skl_ind);
vis_local_image = crop_bbox3(vessel_image, link_global_bbox, 'default');
vis_local_mask = crop_bbox3(vessel_mask, link_global_bbox, 'default');
% Plot
subplot(1,5,1:2)
plot3(link_local_sub(:,2),link_local_sub(:,1), link_local_sub(:,3),'-g', 'LineWidth', 3);
hold on
scatter3(vis_skl_sub(:,2),vis_skl_sub(:,1), vis_skl_sub(:,3), 'k', 'filled');
tmp_fv = isosurface(vis_local_mask,false);
patch(tmp_fv, 'FaceColor', 'red', 'EdgeColor', 'none', 'FaceAlpha', 0.25);
xlabel('X');
ylabel('Y');
zlabel('Z');
hold on
tmp_max_proj_xy = max(vis_local_image, [], 3);
tmp_max_proj_yz = squeeze(max(vis_local_image, [], 1));
tmp_max_proj_xz = squeeze( max(vis_local_image, [], 2));
vis_local_image(:,:,1) = tmp_max_proj_xy;
vis_local_image(:,end,:) = tmp_max_proj_xz;
vis_local_image(end,:,:) = tmp_max_proj_yz;
slice(single(vis_local_image), link_global_bbox(5), link_global_bbox(4), 1)
axis equal

subplot(1,5,3)
imagesc(tmp_max_proj_xy);
hold on
plot(link_local_sub(:,2), link_local_sub(:,1), '-r', 'LineWidth', 3)
hold on
scatter(vis_skl_sub(:,2), vis_skl_sub(:,1), 'k', 'filled');
xlabel('X');
ylabel('Y');
daspect([1,1,1]);
subplot(1,5,4)
imagesc(tmp_max_proj_yz);
hold on
plot(link_local_sub(:,3), link_local_sub(:,2), '-r', 'LineWidth', 3)
scatter(vis_skl_sub(:,3), vis_skl_sub(:,2), 'k', 'filled');
xlabel('Z');
ylabel('X');
daspect([1,1,1]);
subplot(1,5,5)
imagesc(tmp_max_proj_xz);
hold on
plot(link_local_sub(:,3), link_local_sub(:,1), '-r', 'LineWidth', 3)
scatter(vis_skl_sub(:,3), vis_skl_sub(:,1), 'k', 'filled');
xlabel('Z');
ylabel('Y');
daspect([1,1,1]);

fig_hdl.Visible = 'on';
end
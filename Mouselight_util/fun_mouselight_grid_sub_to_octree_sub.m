function [octree_sub_str] = fun_mouselight_grid_sub_to_octree_sub(grid_sub, grid_voxel_size, octree_str)
% Use this function to compute the position of a point in the 240-cube
% space to the position in the rendered image octree. 
% Input: 
%   mmxx: 3-by-1 numerical vector. Position of the point in the 240-cube
%   space
%   
grid_sub_um = grid_sub .* grid_voxel_size;
grid_sub_r_pxl = grid_sub_um ./ octree_str.voxel_size;
grid_sub_r_pxl = max(1, min(octree_str.data_size, round(grid_sub_r_pxl)));
% Octree subscript
octree_grid_sub = ceil(grid_sub_r_pxl ./ octree_str.block_size);
assert(all(octree_grid_sub <= octree_str.grid_size));
% Pixel subscript in the octree image tile
octree_bbox_mm = (octree_grid_sub - 1) .* octree_str.block_size + 1;
tile_sub = grid_sub_r_pxl - octree_bbox_mm + 1;
assert(all(tile_sub > 0) & all(tile_sub <= octree_str.block_size));
octree_bbox_xx = octree_bbox_mm + octree_str.block_size - 1;

octree_sub_str.render_grid_sub = octree_grid_sub;
octree_sub_str.render_space_sub = grid_sub_r_pxl;
octree_sub_str.render_in_tile_sub = tile_sub;
octree_sub_str.render_space_bbox_mm = octree_bbox_mm;
octree_sub_str.render_space_bbox_xx = octree_bbox_xx;
% 
octree_sub_str.grid_space_sub = grid_sub;
octree_sub_str.grid_space_bbox_mm = max(1, floor(octree_bbox_mm .* octree_str.voxel_size ./ grid_voxel_size));
octree_sub_str.grid_space_bbox_xx = ceil(octree_bbox_xx.* octree_str.voxel_size ./ grid_voxel_size);
end
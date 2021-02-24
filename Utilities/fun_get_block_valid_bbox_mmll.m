function bbox_mmll = fun_get_block_valid_bbox_mmll(grid_sub, grid_size, block_size, overlap)

bbox_mmxx = fun_get_block_valid_bbox_mmxx(grid_sub, grid_size, block_size, overlap);
bbox_mmll = bbox_mmxx;
bbox_mmll(4:6) = bbox_mmxx(4:6) - bbox_mmxx(1:3) + 1;
end
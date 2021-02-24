function output = fun_generate_grid(block_size, block_overlap_l, data_size)
% CREATE_GRID creates grid of length block_size(squared grid for 2d, cube
% grid for 3 d) and overlap block_overlap_l. data_size specify the size of
% the whole grid. 

% block_half_size = block_size / 2 ;
data_size = double(data_size);
data_dim = length(data_size);
% Left - up points for each box
upper_left_pos = cell(size(data_size));
grid_bbox_length = cell(size(upper_left_pos));
if isscalar(block_size)
    block_size = block_size * ones(1, data_dim);
end
if isscalar(block_overlap_l)
    block_overlap_l = block_overlap_l * ones(1, data_dim);
end
for idim = 1 : data_dim
    upper_left_pos{idim} = 1 : (block_size(idim) - block_overlap_l(idim)) : data_size(idim);
end
switch data_dim
    case 3
        [grid_ul_1, grid_ul_2, grid_ul_3] = ndgrid(upper_left_pos{1}, upper_left_pos{2}, upper_left_pos{3});
        grid_ul_pos = [grid_ul_1(:) grid_ul_2(:) grid_ul_3(:)];
    case 2
        [grid_ul_1, grid_ul_2] = ndgrid(upper_left_pos{1}, upper_left_pos{2});
        grid_ul_pos = [grid_ul_1(:) grid_ul_2(:)];
    case 1
        [grid_ul_1] = ndgrid(upper_left_pos{1});
        grid_ul_pos = [grid_ul_1(:)];
end
% grid_bbox_center = grid_ul_pos + block_half_size -1;

output.ul = grid_ul_pos;
output.mmxx = [grid_ul_pos, min(grid_ul_pos + block_size - 1, data_size)];% [row_min, col_min, row_max, col_max]
output.mmll = [grid_ul_pos, output.mmxx(:,(data_dim+1):(2*data_dim)) - grid_ul_pos + 1 ];% [row, col, l1, l2]
output.center_pos = round(output.mmxx(:,1:data_dim) + output.mmxx(:, (data_dim + 1): end)/2); 
output.mm_idx_pos =  bsxfun(@rdivide, (grid_ul_pos - 1), (block_size - block_overlap_l)) + 1;
output.ll = [output.mmxx(:,(data_dim+1):(2*data_dim)) - grid_ul_pos + 1];
output.grid_size = cellfun(@length, upper_left_pos);
output.num_block = prod(output.grid_size);
output.size = max(output.mm_idx_pos, [], 1);

output.mmll_array = reshape(permute(output.mmll, [2,1]), [size(output.mmll, 2), output.grid_size]);
output.mmxx_array = reshape(permute(output.mmxx, [2,1]), [size(output.mmxx, 2), output.grid_size]);
output.ll_array = reshape(permute(output.ll, [2,1]), [size(output.ll, 2), output.grid_size]);

end
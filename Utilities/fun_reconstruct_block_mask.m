function block_mask = fun_reconstruct_block_mask(mask_str)
% This function convert the mask structure to the mask. 
% Input: 
%   mask_str: struct with fileds:
%       block_size: 3-by-1 double precision array, specify the size of the
%       block 
%       idx1, idx2, idx3: num_mask_pixel-by-1 uint8 arrays, each one
%       specify the position of the mask voxel in the block. The size of
%       the block is chosen to be 240, which fits the range of uint8. 
%       Other fields are not used in this function. 
block_mask = false(mask_str.block_size);
if isfield(mask_str, 'idx1')
%     block_ind_list = sub2ind(mask_str.block_size, mask_str.idx1, mask_str.idx2, mask_str.idx3);
    block_mask(mask_str.idx1, mask_str.idx2, mask_str.idx3) = true;
elseif isfield(mask_str, 'ind')
    block_mask(mask_str.ind) = true;
end
end

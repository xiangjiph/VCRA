function [recon_voxel_ind, recon_voxel_label ]= fun_skeleton_reconstruction_ind_w_label(pos_ind, pos_label, block_size, recon_strel)
% 
% Input: 
%   pos_ind: N-by-1 numerical array, indice of the voxel in the block 
%   pos_label: N-by-1 nuemricla array, label of the voxel( for example, if
%   a voxel is a part of a link, its label can be the label of the link)
%   block_size: 3-by-1 numerical array, the size of the block 
%   strel_array: 3 dimensional logical array. The structure element for
%   reconstruction. 
% Output:
%   recon_mask: 3 dimensional logical array. 
% 
% Implemented by Xiang Ji on 02/20/2019

% Compute the position of the voxel in the structure element
strel_size = size(recon_strel);
% num_block_voxel = prod(block_size);
[pos1, pos2, pos3] = ind2sub(strel_size, find(recon_strel));
% Pad the array
strel_r = max(floor(size(recon_strel)/2));
block_size_pad = block_size + strel_r * 2;
recon_strel_add_pad = sub2ind(block_size_pad, pos1 + strel_r, pos2 + strel_r, pos3 + strel_r);
recon_strel_add_pad = recon_strel_add_pad - recon_strel_add_pad(ceil(length(recon_strel_add_pad)/2));
% Subscript of the centerline voxel in the block 
[pos1, pos2, pos3] = ind2sub(block_size, pos_ind);
% Indices of the centerline voxel in the padded block
pos_ind_pad = sub2ind(block_size_pad, pos1 + strel_r, pos2 + strel_r, pos3 + strel_r);
% Indices of the reconstructed mask in the padded block 
recon_voxel_ind_pad = bsxfun(@plus, pos_ind_pad, recon_strel_add_pad');
% Subscripts of the reconstructed mask in the padded block 
[pos1, pos2, pos3] = ind2sub(block_size_pad, recon_voxel_ind_pad(:));
% Subscripts of the reconstructed mask in the bounding box of the original
% block in the padded block 
select_Q = (pos1 > strel_r & pos1 <= (block_size(1) + strel_r)) & ...
    (pos2 > strel_r & pos2 <= (block_size(2) + strel_r)) & ...
    (pos3 > strel_r & pos3 <= (block_size(3) + strel_r));
% Subscript of the 
pos1 = pos1(select_Q) - strel_r;
pos2 = pos2(select_Q) - strel_r;
pos3 = pos3(select_Q) - strel_r;
recon_voxel_ind = sub2ind(block_size, pos1, pos2, pos3);

recon_voxel_label = repelem(pos_label, 1, numel(recon_strel_add_pad));
recon_voxel_label = recon_voxel_label(select_Q);
end

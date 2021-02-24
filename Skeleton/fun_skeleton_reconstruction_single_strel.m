function recon_mask = fun_skeleton_reconstruction_single_strel(pos_ind, block_size, recon_strel)
% 
% Input: 
%   pos_ind: N-by-1 numerical array, indice of the voxel in the block 
%   block_size: 3-by-1 numerical array, the size of the block 
%   strel_array: 3 dimensional logical array. The structure element for
%   reconstruction. 
% Output:
%   recon_mask: 3 dimensional logical array. 
% 
% This function is equalivent to imdilate(skeleton, strel). This function
% is faster than imdilate if the structure element is large( say more than
% 5x5x5)

% Compute the position of the voxel in the structure element

if nargin < 3
    recon_mask = false(block_size);
    recon_mask(pos_ind) = true;
    return;
else
    strel_size = size(recon_strel);
    [pos1, pos2, pos3] = ind2sub(strel_size, find(recon_strel));
    pad_size = floor(strel_size/2);
    % Initialize the mask (padded to avoid exceeding the boundary)
    mask_size_pad = block_size + 2 * pad_size;
    recon_mask = false(mask_size_pad);
    % Compute the relative index difference for each voxels in the strel with
    % respect to the center
    recon_strel_add = sub2ind(mask_size_pad, pos1, pos2, pos3);
    recon_strel_add = recon_strel_add - recon_strel_add(ceil(length(recon_strel_add)/2));
    % Indices of voxel in the padded array
    [pos1, pos2, pos3] = ind2sub(block_size, pos_ind);
    voxel_ind_pad = sub2ind(mask_size_pad, pos1 + pad_size(1), ...
        pos2 + pad_size(2), pos3 + pad_size(3));
    % Indices of reconstructed voxels in the padded array
    recon_voxel_ind = bsxfun(@plus, voxel_ind_pad, recon_strel_add');
    recon_mask(recon_voxel_ind(:)) = true;
    recon_mask = recon_mask(pad_size(1)+1:end-pad_size(1), ...
        pad_size(2)+1:end-pad_size(2), pad_size(3)+1:end-pad_size(3));
end

end

function mask_label = fun_skeleton_reconstruction_label(skeleton_voxel, skeleton_radius, skeleton_label, block_size)
% fun_skeleton_reconstruction_label reconstruct the vessel labeled array. 
% Input: 
%   skeleton_voxel: numerical vector, linear indices location of the vessel
%   centerline
%   skeleton_radius: numerical vector, radius of the vessel center line
%   voxle 
%   skeleton_lable: lable of the centerline voxel, can be the label of the
%   link and node
%   block_size: 3-by-1 numerical vector, size of the mask for
%   reconstruction 
% Output: 
%   mask_label: reconstruction of the vessel mask, with each voxel labeled
%   with the label of the skeleton voxel from which it is reconstructed. 
%
% Implemented by Xiang Ji on 02/20/2019   

% Classify skeleton voxel
if isempty(skeleton_label)
    mask_label = zeros(block_size, 'int8');
    return
end
label_max = max(abs(skeleton_label));
if label_max < 127
    mask_label = zeros(block_size, 'int8');
elseif label_max < 32767
    mask_label = zeros(block_size, 'int16');
elseif label_max < 2147483647
    mask_label = zeros(block_size, 'int32');
else
    mask_label = zeros(block_size, 'int64');
end

% Classify skeleton voxel
max_r_edge = ceil(max(skeleton_radius));
radius_bin_val = [1, sqrt(2), 2, sqrt(5), 3 : max_r_edge];
radius_bin_edge = [0, movmean(radius_bin_val, 2, 'Endpoints', 'discard'), max_r_edge];

r_class_ind_cell = fun_bin_data_to_idx_list_by_edges(skeleton_radius, radius_bin_edge, true);
assert(sum(cellfun(@numel, r_class_ind_cell)) == numel(skeleton_voxel), 'Number of skeleton voxel in the cell array does not equal the number of input skeleton voxel')

for tmp_idx = 1 : numel(r_class_ind_cell)
    tmp_ind = r_class_ind_cell{tmp_idx};
    if ~isempty(tmp_ind)
        tmp_radius = radius_bin_val(tmp_idx);   
        switch tmp_radius
            case radius_bin_val(1)
                tmp_label = skeleton_label(tmp_ind);
                tmp_ind = skeleton_voxel(tmp_ind);
            case radius_bin_val(2)
                [tmp_ind, tmp_label] = fun_skeleton_reconstruction_ind_w_label(skeleton_voxel(tmp_ind),...
                    skeleton_label(tmp_ind), block_size, strel('sphere',1).Neighborhood);                
            case radius_bin_val(3)
                [tmp_ind, tmp_label] = fun_skeleton_reconstruction_ind_w_label(skeleton_voxel(tmp_ind),...
                    skeleton_label(tmp_ind), block_size, strel('cube', 3).Neighborhood);
            case radius_bin_val(4)
                [tmp_ind, tmp_label] = fun_skeleton_reconstruction_ind_w_label(skeleton_voxel(tmp_ind),...
                    skeleton_label(tmp_ind), block_size, strel('sphere', 2).Neighborhood);
            otherwise
                [tmp_ind, tmp_label] = fun_skeleton_reconstruction_ind_w_label(skeleton_voxel(tmp_ind),...
                    skeleton_label(tmp_ind), block_size, strel('sphere', double(tmp_radius)).Neighborhood);
        end
        mask_label(tmp_ind) = tmp_label;
    end
end    
end


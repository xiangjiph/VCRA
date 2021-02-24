function mask_value = fun_skeleton_reconstruction_with_value(skeleton_voxel, skeleton_radius, skeleton_value, block_size)
% fun_skeleton_reconstruction_with_value reconstruct the vessel labeled array. 
% Input: 
%   skeleton_voxel: numerical vector, linear indices location of the vessel
%   centerline
%   skeleton_radius: numerical vector, radius of the vessel center line
%   voxle 
%   skeleton_value: value of the centerline voxel, can be blood flow speed
%   block_size: 3-by-1 numerical vector, size of the mask for
%   reconstruction 
% Output: 
%   mask_label: reconstruction of the vessel mask, with each voxel labeled
%   with the label of the skeleton voxel from which it is reconstructed. 
%
% Implemented by Xiang Ji on 02/20/2019   

% Classify skeleton voxel
if isempty(skeleton_value)
    mask_value = zeros(block_size);
    return
end
skeleton_radius = round(skeleton_radius);
label_max = max(abs(skeleton_value));
if isinteger(skeleton_value(1))
    if label_max < intmax('int8')
        mask_value = zeros(block_size, 'int8');
    elseif label_max < intmax('int16')
        mask_value = zeros(block_size, 'int16');
    elseif label_max < intmax('int32')
        mask_value = zeros(block_size, 'int32');
    else
        mask_value = zeros(block_size, 'int64');
    end
else
    mask_value = zeros(block_size, 'like', skeleton_value);
end
radius_list = unique(skeleton_radius);
% Binary reconstruction
for tmp_idx = 1 : numel(radius_list)
    tmp_radius = radius_list(tmp_idx);
    if tmp_radius >= 1
        tmp_Q = (skeleton_radius == tmp_radius);
        [tmp_ind, tmp_label]= fun_skeleton_reconstruction_ind_w_label(skeleton_voxel(tmp_Q),  ...
            skeleton_value(tmp_Q), block_size, strel('sphere', double(tmp_radius)).Neighborhood);
        mask_value(tmp_ind) = tmp_label;
    end
end
    
end


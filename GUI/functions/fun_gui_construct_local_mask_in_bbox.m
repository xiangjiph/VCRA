function [local_mask, varargout] = fun_gui_construct_local_mask_in_bbox(mask_global_sub, bbox_mmxx)
% fun_gui_construct_local_mask_in_bbox takes the global mask subscript and
% the boudning box to construct the local mask inside the boundign box. 
% Input: 
%   mask_global_sub: N-by-3 numerical array, 3D subscript of the voxel in
%   the mask 
%   bbox_mmxx: 6-element vector, min and max coordinate of the bounding
%   box. 
% Output:
%   local_mask_str: MATLAB structure
% 
% Implemented by Xiang Ji on Jul 18, 2019

bbox_size = bbox_mmxx(4:6) - bbox_mmxx(1:3) + 1;
in_bbox_Q = fun_voxel_sub_in_bbox_mmxx_Q(mask_global_sub, bbox_mmxx);
in_bbox_global_sub = mask_global_sub(in_bbox_Q, :);

in_bbox_local_sub = in_bbox_global_sub - bbox_mmxx(1:3) + 1;
in_bbox_local_ind = sub2ind(bbox_size, in_bbox_local_sub(:, 1), ...
    in_bbox_local_sub(:, 2), in_bbox_local_sub(:, 3));

info_str = struct;
info_str.global_bbox_mmxx = bbox_mmxx;
info_str.bbox_size = bbox_size;
info_str.local_ind = in_bbox_local_ind;
info_str.global_sub = in_bbox_global_sub;

local_mask = false(bbox_size);
local_mask (in_bbox_local_ind) = true;

if nargout > 1
    varargout{1} = info_str;
end

end
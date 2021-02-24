function [local_linear_ind_t , vargout]= fun_linear_index_coordiante_transform3D(local_linear_ind, bbox_cur_in_target_mmll, bbox_target_mmll)
% This function convert the array indices in one block to the array indices
% in the other block. 
% Input: 
%     local_linear_ind: N-by-1 double precision array

%     bbox_cur_in_target_mmll: the bounding box from which
%     local_linear_ind is subtracted, [min_pos1, min_pos2, min_pos3, l1,
%     l2, l3]. min_pos1 is the minimum index in the first dimensiton (row)
%     in the coordinate of the target

%     bbox_target_mmll: normally equals [1,1,1, size_of_bbox]

if any(bbox_cur_in_target_mmll(1:3) < bbox_target_mmll(1:3) | ...
        (bbox_cur_in_target_mmll(1:3) + bbox_cur_in_target_mmll(4:6)) > ...
        (bbox_target_mmll(1:3) + bbox_target_mmll(4:6)))
    error('The current bouding box is out of the range of the target bounding box');
end

[pos1, pos2, pos3] = ind2sub(bbox_cur_in_target_mmll(4:6), local_linear_ind);
pos1_t = pos1 + bbox_cur_in_target_mmll(1) - bbox_target_mmll(1);
pos2_t = pos2 + bbox_cur_in_target_mmll(2) - bbox_target_mmll(2);
pos3_t = pos3 + bbox_cur_in_target_mmll(3) - bbox_target_mmll(3);
local_linear_ind_t = sub2ind(bbox_target_mmll(4:6), pos1_t, pos2_t, pos3_t);

if nargout > 1
    vargout = cat(2, pos1_t, pos2_t, pos3_t);
end

end

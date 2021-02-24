function surface_area_list = fun_cc_compute_surface_area(cc)
% Function for computing the surface area from the connected component
% (cc). 
% Input: 
%   cc: connected components, output by bwconcomp
% Output:
%   surface_area_list: list of surface area for each connected component in
%       cc
% 
%
%
labeled_matrix = labelmatrix(cc);
% Convert uint type to int type or double for the following computation
num_cc = cc.NumObjects;
if num_cc < intmax('int8')
    labeled_matrix = int8(labeled_matrix);
elseif num_cc < intmax('int16')
    labeled_matrix = int16(labeled_matrix);
elseif num_cc < intmax('int32')
    labeled_matrix = int32(labeled_matrix);
elseif num_cc < intmax('int64')
    labeled_matrix = int64('int64');
else
    labeled_matrix = double(labeled_matrix);
end
% Pad the array for computing the gradient of the boundary voxels. 
mask_labeled_pad = padarray(labeled_matrix, [1,1,1], 0, 'both');
% Compute the gradient in three dimension, correspond to the surface facing
% three direction respectively. 
[lm_d1, lm_d2, lm_d3] = fun_gradient3D(mask_labeled_pad, [1,2,3], 'intermediate');
% Combine both positive and negative gradient value as they correspond to
% the same connected component. 
lm_d1_abs = abs(lm_d1);
lm_d2_abs = abs(lm_d2);
lm_d3_abs = abs(lm_d3);
% Count the histogram in three gradient array. Each bin correspond to the
% surface area of each connected components. 
surface_area_list = histcounts(lm_d1_abs(lm_d1_abs>0), num_cc) + ...
    histcounts(lm_d2_abs(lm_d2_abs>0), num_cc) + ...
    histcounts(lm_d3_abs(lm_d3_abs>0), num_cc);

end
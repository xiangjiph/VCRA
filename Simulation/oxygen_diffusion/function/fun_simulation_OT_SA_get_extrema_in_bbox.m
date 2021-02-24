function extrema_str = fun_simulation_OT_SA_get_extrema_in_bbox(extrema_input_str, bbox_mmxx, upsample_rate)
% Get the local extrema in the bounding box and convert the indices and
% subscripts of the local extrama to the local coordinate inside the
% boudning box. 
assert(upsample_rate >= 1 && mod(upsample_rate, 1) == 0);

extrema_str = struct;
num_cell = numel(extrema_input_str.ind);
[extrema_str.ind, extrema_str.sub, extrema_str.dt, extrema_str.pO2] = ...
    deal(cell(num_cell, 1));

bbox_ll = bbox_mmxx(4:6) - bbox_mmxx(1:3) + 1;
for iter_cell = 1 : num_cell
   tmp_sub = extrema_input_str.sub{iter_cell};
   in_bbox_Q = fun_voxel_sub_in_bbox_mmxx_Q(tmp_sub, bbox_mmxx);
   extrema_str.dt{iter_cell} = extrema_input_str.dt{iter_cell}(in_bbox_Q);
   extrema_str.pO2{iter_cell} = extrema_input_str.pO2{iter_cell}(in_bbox_Q);
   tmp_sub = tmp_sub(in_bbox_Q, :);
   tmp_sub_local = ceil((tmp_sub - bbox_mmxx(1:3) + 1) .* upsample_rate);
   tmp_ind = sub2ind(bbox_ll * upsample_rate, tmp_sub_local(:, 1), tmp_sub_local(:, 2), ...
       tmp_sub_local(:, 3));
   extrema_str.ind{iter_cell} = tmp_ind;
   extrema_str.sub{iter_cell} = tmp_sub_local;    
end
end
function volume_fraction = fun_analysis_get_bbox_in_mask_vol_ratio(mask, bbox_mmxx)


num_bbox = size(bbox_mmxx, 1);
volume_fraction = zeros(num_bbox, 1);

mask_ind = find(mask);
if isempty(mask_ind)
    return;
end
mask_size = size(mask);
mask_sub = fun_ind2sub(mask_size, mask_ind);
mask_sub_min = min(mask_sub, [], 1);
mask_sub_max = max(mask_sub, [], 1);

bbox_mmxx(:, 1:3) = bsxfun(@max, round(bbox_mmxx(:, 1:3)), [1,1,1]);
bbox_mmxx(:, 4:6) = bsxfun(@min, round(bbox_mmxx(:, 4:6)), mask_size);
bbox_mmll(:, 1:3) = bbox_mmxx(:, 1:3);
bbox_mmll(:, 4:6) = bbox_mmxx(:, 4:6) - bbox_mmxx(:, 1:3) + 1;

sub_pair_1 = [1, 4];
sub_pair_2 = [2, 5];
sub_pair_3 = [3, 6];
in_bbox_Q = false(num_bbox, 1);
for iter_1 = 1 : 2
    for iter_2 = 1 : 2
        for iter_3 = 1 : 2
            tmp_is_in_bbox_Q = fun_voxel_sub_in_bbox_mmxx_Q(bbox_mmxx(:, ...
                [sub_pair_1(iter_1), sub_pair_2(iter_2), sub_pair_3(iter_3)]),...
                [mask_sub_min, mask_sub_max]);
            in_bbox_Q = in_bbox_Q | tmp_is_in_bbox_Q;
        end
    end
end
in_mask_bbox_idx = find(in_bbox_Q);

bbox_mmll = bbox_mmll.';
for iter_bbox = 1 : numel(in_mask_bbox_idx)
    tmp_bbox_idx =  in_mask_bbox_idx(iter_bbox);
    tmp_mask = crop_bbox3(mask, bbox_mmll(:, tmp_bbox_idx));
    volume_fraction(tmp_bbox_idx) = nnz(tmp_mask) ./ numel(tmp_mask);    
end
end
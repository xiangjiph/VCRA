function [output, varargout] = crop_center_box(data, center_pos, radius)
% crop_center_box returns the subarray of the input data with center at
% CENTER_POS and size radius*2+1;
center_pos = round(center_pos);
index_max_limits = size(data);

if isvector(data)
    bbox_mmxx = [max(1,center_pos - radius), min(length(data), center_pos + radius)];
    output = data(bbox_mmxx(1):bbox_mmxx(2));
else
    data_dim = length(index_max_limits);
    sub_min = max(1, center_pos - radius);
    sub_max = min(index_max_limits, center_pos + radius);
    bbox_mmxx = [sub_min, sub_max];
    switch data_dim
        case 2
            output = data(sub_min(1):sub_max(1), sub_min(2):sub_max(2));
        case 3
            output = data(sub_min(1):sub_max(1), sub_min(2):sub_max(2), sub_min(3):sub_max(3));
    end
    if nargout > 1
        varargout{1} = bbox_mmxx;
    end
end
end

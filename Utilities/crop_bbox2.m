function output = crop_bbox2(data, bbox_parameters, bbox_order)
% [data] is a 2D array
% default bbox_parameters = [ul1, ul2, l1, l2]
% matlab's regionpros3 output bbox = [ul2, ul1, l2, l1]
% Can be modified to crop the image while pad 0 to left part. 
if nargin < 3
    bbox_order = 'default';
end
if ~iscell(bbox_order)
    bbox_parameters = num2cell(round(bbox_parameters));
end
switch bbox_order
    case {'default'}
        if length(bbox_parameters) <=3
            [row_min, col_min, l_row] = bbox_parameters{:};
            l_col = l_row;
        else
            [row_min, col_min, l_row, l_col] = bbox_parameters{:};
        end
    case {'regionprop'}
        if length(bbox_parameters) <=3
            [col_min, row_min, l_row] = bbox_parameters{:};
            l_col = l_row;
        else
            [col_min, row_min, l_col, l_row] = bbox_parameters{:};
        end
end
data_size = size(data);
max_idx = min([(row_min + l_row - 1) , (col_min + l_col - 1)], data_size); % The space in this line is very tricky. [ a, b -1] = [a, b, -1]!!!!
output = data(row_min: max_idx(1), col_min:max_idx(2));

end
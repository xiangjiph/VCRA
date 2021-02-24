function output = crop_bbox3(data, bbox_mmll, bbox_order)
% CROP_BBOX3 crops part of the array DATA according to the given bounding
% box parameters. 
% default bbox_parameters = [ul1, ul2, ul3, l1, l2, l3]
% matlab's regionpros3 output bbox = [ul2, ul1, ul3, l2, l1, l3]
if nargin < 3
    bbox_order = 'default';
end
if ~iscell(bbox_order)
    bbox_mmll = num2cell(round(bbox_mmll));
end
switch bbox_order
    case {'default'}
        if length(bbox_mmll) <=4
            [ul1, ul2, ul3, l1] = bbox_mmll{:};
            l2 = l1;
            l3 = l1;
        else
            [ul1, ul2, ul3, l1, l2, l3] = bbox_mmll{:};
        end
    case {'regionprop'}
        if length(bbox_mmll) <=4
            [ul2, ul1, ul3, l1] = bbox_mmll{:};
            l2 = l1;
            l3 = l1;
        else
            [ul2, ul1, ul3, l2, l1, l3] = bbox_mmll{:};
        end
end
output = data(ul1:ul1+l1-1, ul2:ul2+l2-1, ul3:ul3+l3-1);
% data_size = size(data);
% max_idx = min([(row_min + l_row - 1) , (col_min + l_col - 1)], data_size); % The space in this line is very tricky. [ a, b -1] = [a, b, -1]!!!!
% output = data(row_min: max_idx(1), col_min:max_idx(2));
end
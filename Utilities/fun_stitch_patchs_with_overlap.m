function [stitched_block, varargout] = fun_stitch_patchs_with_overlap(input_cell_array, overlap)
% Input:
%   block_cell_array: a cell array that organize the block to be stitched
%   overlpap: scalar non-negative value
% Output:
%   block_stitched: numerical array.
%   stitch_info: structure with fields: 
%       stitched_block_bbox: the minimum *3 and maximum *3 coordinates of the
%           voxels that are used for stitching.
%       stitch_block_size: cell array of the size of the blocks that are used
%           for stitching
%       Position of the origin of each cropped block in the combined array
%       Position of the origin of each block in the combined array
%       Bounding box for each cropped block in the block
%       Bounding box for each block in the combined array
%       Bounding box for each cropped block in the combined array
% Implemented by Xiang Ji on 03/18/2019
% Modified by Xiang Ji on 03/21/2019
% 1. Fix a bug when all the cells in a row / column are empty.
if isscalar(overlap)
    overlap = [overlap, overlap];
end
stitch_info(size(input_cell_array)) = struct;
half_overlap_size = round(overlap/2);
[num_row, num_col] = size(input_cell_array);
is_empty_Q = cellfun(@isempty, input_cell_array);
not_empty_idx_1 = find(~is_empty_Q, 1);
data_type = class(input_cell_array{not_empty_idx_1});
im_size_1 = cellfun(@(x) size(x, 1), input_cell_array);
im_size_2 = cellfun(@(x) size(x, 2), input_cell_array);
im_size_1_max = max(im_size_1, [], 2);
im_size_2_max = max(im_size_2, [], 1);
%% Infer the size of the tiles in the empty cell
empty_idx_1 = find(im_size_1_max == 0);
if ~isempty(empty_idx_1)
    num_empty_idx = numel(empty_idx_1);
    for iter_idx = 1 : num_empty_idx
        tmp_idx = empty_idx_1(iter_idx);
        if tmp_idx == 1 || tmp_idx == num_row
            warning('The first or the last row is completely empty.');
        end
            warning('Exist empty row. Infer size automatically');
        non_empty_idx = find(im_size_1_max);
        non_empty_idx = setdiff(non_empty_idx, [1, num_row]);
        % Use the maximum size of the closest nonemtpy patch size
        [~, min_idx] = min(abs(non_empty_idx - tmp_idx));
        im_size_1_max(tmp_idx) = im_size_1_max(non_empty_idx(min_idx));
    end
end

empty_idx_2 = find(im_size_2_max == 0);
if ~isempty(empty_idx_2)
    num_empty_idx = numel(empty_idx_2);
    for iter_idx = 1 : num_empty_idx
        tmp_idx = empty_idx_2(iter_idx);
        if tmp_idx == 1 || tmp_idx == num_row
            warning('The first or the last column is completely empty.');
        end
            warning('Exist empty column. Infer size automatically');
        non_empty_idx = find(im_size_2_max);
        non_empty_idx = setdiff(non_empty_idx, [1, num_col]);
        % Use the size of the closest nonemtpy patch to replace the empty
        % one
        [~, min_idx] = min(abs(non_empty_idx - tmp_idx));
        im_size_2_max(tmp_idx) = im_size_2_max(non_empty_idx(min_idx));
    end
end
%%
im_size_1_fill = bsxfun(@max, im_size_1, im_size_1_max);
im_size_2_fill = bsxfun(@max, im_size_2, im_size_2_max);

for col_idx = 1 : num_col
    for row_idx = 1 : num_row
        if is_empty_Q(row_idx, col_idx)
            % Fill the empty cell with zeros
            block_data = zeros(im_size_1_fill(row_idx, col_idx), ...
                im_size_2_fill(row_idx, col_idx), data_type);
        else
            block_data = input_cell_array{row_idx, col_idx};
        end
        
        block_size_before_crop = size(block_data);
        bbox_in_block_mmll = [1,1, size(block_data)];
        
        block_size = zeros(1,2);
        bbox_in_array_mmll = zeros(1,4);
        
        bbox_before_crop_in_array_mmll = zeros(1,4);
        % Bounding box for each block in the combined array
        bbox_before_crop_in_array_mmll(1) = (row_idx - 1) * (block_size_before_crop(1) - overlap(1)) + 1;
        bbox_before_crop_in_array_mmll(2) = (col_idx - 1) * (block_size_before_crop(2) - overlap(2)) + 1;
        bbox_before_crop_in_array_mmll(3:4) = block_size_before_crop;
        bbox_before_crop_in_array_mmxx = bbox_before_crop_in_array_mmll;
        bbox_before_crop_in_array_mmxx(3:4) = bbox_before_crop_in_array_mmll(1:2) + bbox_before_crop_in_array_mmll(3:4) - 1;
        
        if (row_idx ~= 1) && (row_idx ~= num_row)
            block_size(1) = block_size_before_crop(1) - overlap(1);
            % Bounding box for each cropped block in the block
            bbox_in_block_mmll(1) = half_overlap_size(1) + 1;
            % Bounding box for each cropped block in the combined array:
            bbox_in_array_mmll(1) = (row_idx - 1) * (block_size(1)) + half_overlap_size(1) + 1;
        elseif (row_idx == 1) && (num_row == 1)
            block_size(1) = block_size_before_crop(1);
            bbox_in_block_mmll(1) = 1;
            bbox_in_array_mmll(1) = 1;
        elseif row_idx == 1
            block_size(1) = block_size_before_crop(1) - half_overlap_size(1);
            % Bounding box for each cropped block in the block
            bbox_in_block_mmll(1) = 1;
            % Bounding box for each cropped block in the combined array:
            bbox_in_array_mmll(1) = 1;
        else
            % Bounding box for each cropped block in the block
            block_size(1) = block_size_before_crop(1) - half_overlap_size(1);
            bbox_in_block_mmll(1) = half_overlap_size(1) + 1;
            % Bounding box for each cropped block in the combined array:
            bbox_in_array_mmll(1) = (row_idx - 1) * (block_size(1)) + half_overlap_size(1) + 1;
        end
        
        
        if (col_idx ~= 1) && (col_idx ~= num_col)
            block_size(2) = block_size_before_crop(2) - overlap(2);
            % Bounding box for each cropped block in the block
            bbox_in_block_mmll(2) = half_overlap_size(2) + 1;
            % Bounding box for each cropped block in the combined array:
            bbox_in_array_mmll(2) = (col_idx - 1) * (block_size(2)) + half_overlap_size(2) + 1;
        elseif (col_idx == 1) && (num_col == 1)
            block_size(2) = block_size_before_crop(2);
            bbox_in_block_mmll(2) = 1;
            bbox_in_array_mmll(2) = 1;
        elseif col_idx == 1
            block_size(2) = block_size_before_crop(2) - half_overlap_size(2);
            % Bounding box for each cropped block in the block
            bbox_in_block_mmll(2) = 1;
            % Bounding box for each cropped block in the combined array:
            bbox_in_array_mmll(2) = 1;
        else
            % Bounding box for each cropped block in the block
            block_size(2) = block_size_before_crop(2) - half_overlap_size(2);
            bbox_in_block_mmll(2) = half_overlap_size(2) + 1;
            % Bounding box for each cropped block in the combined array:
            bbox_in_array_mmll(2) = (row_idx - 1) * (block_size(2)) + half_overlap_size(2) + 1;
        end
        bbox_in_block_mmll(3:4) = block_size;
        bbox_in_array_mmll(3:4) = block_size;
        bbox_in_array_mmxx = bbox_in_array_mmll;
        bbox_in_array_mmxx(3:4) = bbox_in_array_mmll(1:2) + bbox_in_array_mmll(3:4) - 1;
        bbox_in_block_mmxx = bbox_in_block_mmll;
        bbox_in_block_mmxx(3:4) = bbox_in_block_mmll(1:2) + bbox_in_block_mmll(3:4) - 1;
        
        block_size = round(block_size);
        bbox_in_block_mmll = round(bbox_in_block_mmll);
        bbox_in_array_mmxx = round(bbox_in_array_mmxx);
        bbox_in_block_mmxx = round(bbox_in_block_mmxx);
        bbox_in_array_mmll = round(bbox_in_array_mmll);
        
        block_data = crop_bbox2(block_data, bbox_in_block_mmll, 'default');
        
        input_cell_array{row_idx, col_idx} = block_data;
        
        stitch_info(row_idx, col_idx).bbox_in_array_mmxx = bbox_in_array_mmxx;
        stitch_info(row_idx, col_idx).bbox_in_array_mmll = bbox_in_array_mmll;
        stitch_info(row_idx, col_idx).bbox_in_block_mmll = bbox_in_block_mmll;
        stitch_info(row_idx, col_idx).bbox_in_block_mmxx = bbox_in_block_mmxx;
        stitch_info(row_idx, col_idx).block_size = block_size;
        stitch_info(row_idx, col_idx).block_size_before_crop = block_size_before_crop;
        stitch_info(row_idx, col_idx).bbox_before_crop_in_array_mmll = bbox_before_crop_in_array_mmll;
        stitch_info(row_idx, col_idx).bbox_before_crop_in_array_mmxx = bbox_before_crop_in_array_mmxx;
        
    end
end

stitched_block = cell2mat(input_cell_array);

if nargout == 2
    varargout = cell(1,1);
    varargout{1} = stitch_info;
end

end
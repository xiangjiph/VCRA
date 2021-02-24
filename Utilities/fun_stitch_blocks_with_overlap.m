function [stitched_block, varargout ]= fun_stitch_blocks_with_overlap(input_cell_array, overlap, z_overlap_correction)

% Input:
%   block_cell_array: a cell array that organize the block to be stitched 
%   overlpap: scalar non-negative value
%   z_overlap_correction: scaling factor in z direction. The block are
%   generated with the original voxel size
%   um). During processing, the voxel size is corrected to be isotropic.
%   Here we are stitching the image stack with [0.3, 0.3, 0.3] um voxel
%   size. The overlap in z direction therefore mush be scale by
%   overlap/0.3;
% Output:
%   block_stitched: numerical array. 
%   stitched_block_bbox: the minimum *3 and maximum *3 coordinates of the
%       voxels that are used for stitching. 
%   stitch_block_size: cell array of the size of the blocks that are used
%       for stitching
%   Position of the origin of each cropped block in the combined array
%   Position of the origin of each block in the combined array
%   Bounding box for each cropped block in the block
%   Bounding box for each block in the combined array
%   Bounding box for each cropped block in the combined array

% input_cell_array = image_cell_array;
% overlap = grid_info.block_overlap;
if nargin < 3
    z_overlap_correction = 1;
end
stitch_info(size(input_cell_array)) = struct;
half_overlap_size = round(overlap/2);
half_overlap_size_z_corrected = round(overlap / (2 * z_overlap_correction));
overlap_z_corrected = 2 * half_overlap_size;
% overlap_z_corrected = round(overlap / (z_overlap_correction));
[num_row, num_col, num_layer] = size(input_cell_array);
for layer_idx = 1 : num_layer
    for col_idx = 1 : num_col
        for row_idx = 1 : num_row
            block_data = input_cell_array{row_idx, col_idx, layer_idx};
            block_size_before_crop = size(block_data);
            bbox_in_block_mmll = [1,1,1, size(block_data)];
            
            block_size = zeros(1,3);
            bbox_in_array_mmll = zeros(1,6);
            
            bbox_before_crop_in_array_mmll = zeros(1,6);
            % Bounding box for each block in the combined array
            bbox_before_crop_in_array_mmll(1) = (row_idx - 1) * (block_size_before_crop(1) - overlap) + 1;
            bbox_before_crop_in_array_mmll(2) = (col_idx - 1) * (block_size_before_crop(2) - overlap) + 1;
            bbox_before_crop_in_array_mmll(3) = (layer_idx - 1) * (block_size_before_crop(3) - overlap_z_corrected) + 1;
            bbox_before_crop_in_array_mmll(4:6) = block_size_before_crop;
            bbox_before_crop_in_array_mmxx = bbox_before_crop_in_array_mmll;
            bbox_before_crop_in_array_mmxx(4:6) = bbox_before_crop_in_array_mmll(1:3) + bbox_before_crop_in_array_mmll(4:6) - 1;
            
            if (row_idx ~= 1) && (row_idx ~= num_row)
                block_size(1) = block_size_before_crop(1) - overlap;
                % Bounding box for each cropped block in the block
                bbox_in_block_mmll(1) = half_overlap_size + 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(1) = (row_idx - 1) * (block_size(1)) + half_overlap_size + 1;
            elseif (row_idx == 1) && (num_row == 1)
                block_size(1) = block_size_before_crop(1);
                bbox_in_block_mmll(1) = 1;
                bbox_in_array_mmll(1) = 1;                
            elseif row_idx == 1
                block_size(1) = block_size_before_crop(1) - half_overlap_size;
                % Bounding box for each cropped block in the block
                bbox_in_block_mmll(1) = 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(1) = 1;
            else
                % Bounding box for each cropped block in the block
                block_size(1) = block_size_before_crop(1) - half_overlap_size;
                bbox_in_block_mmll(1) = half_overlap_size + 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(1) = (row_idx - 1) * (block_size(1)) + half_overlap_size + 1;                
            end
            
            
            if (col_idx ~= 1) && (col_idx ~= num_col)
                block_size(2) = block_size_before_crop(2) - overlap;
                % Bounding box for each cropped block in the block
                bbox_in_block_mmll(2) = half_overlap_size + 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(2) = (col_idx - 1) * (block_size(2)) + half_overlap_size + 1;
            elseif (col_idx == 1) && (num_col == 1)
                block_size(2) = block_size_before_crop(2);
                bbox_in_block_mmll(2) = 1;
                bbox_in_array_mmll(2) = 1;
            elseif col_idx == 1
                block_size(2) = block_size_before_crop(2) - half_overlap_size;
                % Bounding box for each cropped block in the block
                bbox_in_block_mmll(2) = 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(2) = 1;
            else
                % Bounding box for each cropped block in the block
                block_size(2) = block_size_before_crop(2) - half_overlap_size;
                bbox_in_block_mmll(2) = half_overlap_size + 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(2) = (row_idx - 1) * (block_size(2)) + half_overlap_size + 1;                
            end            
            % Problematic part...
            if (layer_idx ~= 1) && (layer_idx ~= num_layer)
                block_size(3) = block_size_before_crop(3) - overlap_z_corrected;
                % Bounding box for each cropped block in the block
                bbox_in_block_mmll(3) = half_overlap_size_z_corrected + 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(3) = (layer_idx - 1) * (block_size(3)) + half_overlap_size_z_corrected + 1;
            elseif (layer_idx == 1) && (num_layer == 1)
                block_size(3) = block_size_before_crop(3);
                bbox_in_block_mmll(3) = 1;
                bbox_in_array_mmll(3) = 1;
            elseif layer_idx == 1
                block_size(3) = block_size_before_crop(3) - half_overlap_size_z_corrected;
                % Bounding box for each cropped block in the block
                bbox_in_block_mmll(3) = 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(3) = 1;
            else
                % Bounding box for each cropped block in the block
                block_size(3) = block_size_before_crop(3) - half_overlap_size_z_corrected;
                bbox_in_block_mmll(3) = half_overlap_size_z_corrected + 1;
                % Bounding box for each cropped block in the combined array:
                bbox_in_array_mmll(3) = (layer_idx - 1) * (block_size(3)) + half_overlap_size_z_corrected + 1;                
            end            
            bbox_in_block_mmll(4:6) = block_size;
            bbox_in_array_mmll(4:6) = block_size;
            bbox_in_array_mmxx = bbox_in_array_mmll;
            bbox_in_array_mmxx(4:6) = bbox_in_array_mmll(1:3) + bbox_in_array_mmll(4:6) - 1;
            bbox_in_block_mmxx = bbox_in_block_mmll;
            bbox_in_block_mmxx(4:6) = bbox_in_block_mmll(1:3) + bbox_in_block_mmll(4:6) - 1;
            
            block_size = round(block_size);
            bbox_in_block_mmll = round(bbox_in_block_mmll);
            bbox_in_array_mmxx = round(bbox_in_array_mmxx);
            bbox_in_block_mmxx = round(bbox_in_block_mmxx);
            bbox_in_array_mmll = round(bbox_in_array_mmll);
            
            block_data = crop_bbox3(block_data, bbox_in_block_mmll, 'default');
            
            input_cell_array{row_idx, col_idx, layer_idx} = block_data;
            
            stitch_info(row_idx, col_idx, layer_idx).bbox_in_array_mmxx = bbox_in_array_mmxx;
            stitch_info(row_idx, col_idx, layer_idx).bbox_in_array_mmll = bbox_in_array_mmll;
            stitch_info(row_idx, col_idx, layer_idx).bbox_in_block_mmll = bbox_in_block_mmll;
            stitch_info(row_idx, col_idx, layer_idx).bbox_in_block_mmxx = bbox_in_block_mmxx;
            stitch_info(row_idx, col_idx, layer_idx).block_size = block_size;
            stitch_info(row_idx, col_idx, layer_idx).block_size_before_crop = block_size_before_crop;
            stitch_info(row_idx, col_idx, layer_idx).bbox_before_crop_in_array_mmll = bbox_before_crop_in_array_mmll;
            stitch_info(row_idx, col_idx, layer_idx).bbox_before_crop_in_array_mmxx = bbox_before_crop_in_array_mmxx;
            
        end
    end
end

stitched_block = cell2mat(input_cell_array);

if nargout == 2
    varargout = cell(1,1);
    varargout{1} = stitch_info;
end

end
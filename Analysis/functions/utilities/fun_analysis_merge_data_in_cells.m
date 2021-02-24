function merge_data = fun_analysis_merge_data_in_cells(cell_array)
% fun_analysis_merge_data_in_cells merges float numerical vectors in cell
% arrays into a matrix, in which each column is from one cell. 
% Input:
%   cell_array: cell array of numerical vector. Since we fill the empty
%   space by NaN, the input cell array should contains double or single
%   values. 
% Output: 
%   merge_data: N-by-M float matrix, where M is the number of cell in the
%   cell_array
% 
% Implemented by Xiang Ji on 02/26/2019

% Find the longest cell
[num_bin, ~] = max(cellfun(@numel, cell_array));
num_cell = numel(cell_array);
merge_data = nan(num_bin, num_cell);
for iter_cell = 1 : num_cell
    tmp_data = cell_array{iter_cell}';
    merge_data(1:numel(tmp_data), iter_cell) = tmp_data;
end
end
function binned_data = fun_bin_data_to_cells_by_ind(data, bin_ind)
% fun_bin_data_to_cells_by_ind reorganize the data into cell arraies
% accoriding to the bin_ind
% Input: 
%   data: numerical vector
%   bin_ind: cell array, indices of the data that need to be put inside
%   each cell
% Output: 
%   binned_data: cell array
%
num_cell = numel(bin_ind);
binned_data = cell(num_cell, 1);
for iter_cell = 1 : num_cell
    binned_data{iter_cell} = data(bin_ind{iter_cell});
end
end
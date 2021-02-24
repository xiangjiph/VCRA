function new_cell_array = fun_concatenate_cell_array_in_arbitrary_dims(cell_array, cat_array_dim, cat_data_dim)

if nargin < 3
    cat_data_dim = 1;
end

if isempty(cat_array_dim)
    new_cell_array = cell_array;
    return;
end
% num_cat_dims = numel(cat_array_dim);
array_size = size(cell_array);
array_dim = ndims(cell_array);
kept_dim_idx = setdiff(1 : array_dim, cat_array_dim, 'stable');

new_array_size = array_size;
new_array_size(cat_array_dim) = 1;

num_elem_cat = prod(array_size(cat_array_dim));
num_elem_uncat = array_size(kept_dim_idx);

cell_array = permute(cell_array, [cat_array_dim, kept_dim_idx]);
num_tmp_col = prod(num_elem_uncat);
cell_array = reshape(cell_array, num_elem_cat, num_tmp_col);
new_cell_array = cell(num_tmp_col, 1);
for iter_cell = 1 : num_tmp_col
    new_cell_array{iter_cell} = cat(cat_data_dim, cell_array{:, iter_cell});
end
if ~isempty(kept_dim_idx)
    new_cell_array = reshape(new_cell_array, new_array_size);
end
end
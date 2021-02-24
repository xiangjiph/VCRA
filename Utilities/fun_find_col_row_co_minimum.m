function [sub_1, sub_2, min_value] = fun_find_col_row_co_minimum(mat, remove_diagonal_Q)
% fun_find_mutual_minimum finds the subscript of the element of a matrix
% who is both the minimum value in row and column.
% Input: 
%   mat: 2D numerical matrix
%   remove_diagonal_Q: if true, do not consider the diagonal term
% Output:
%   sub_1: the first subscript of the element
%   sub_2: the second subscript of the element
%   min_value: the value at the found position
%
if isempty(mat)
    [sub_1, sub_2, min_value] = deal([]);
    return;
end
if nargin < 2
    remove_diagonal_Q = false;
end

if remove_diagonal_Q
    mat_max = max(mat(:)) + 1;
    % Shift the value to be larger than the greatest value in the orignial
    % matrix
    mat = mat + eye(size(mat)) .* mat_max;
end
[~, e2_nearest_e1_idx] = sort(mat, 1, 'ascend');
e2_nearest_e1_idx = e2_nearest_e1_idx(1, :);
[e1_nearest_e2_dist, e1_nearest_e2_idx] = sort(mat(e2_nearest_e1_idx',:), 2, 'ascend');
e1_nearest_e2_idx = e1_nearest_e2_idx(:, 1);
e1_nearest_e2_dist = e1_nearest_e2_dist(:, 1);

sub_2 = find([1 : numel(e1_nearest_e2_idx)]'== e1_nearest_e2_idx);
sub_1 = e2_nearest_e1_idx(sub_2)';
min_value = e1_nearest_e2_dist(sub_2);

end

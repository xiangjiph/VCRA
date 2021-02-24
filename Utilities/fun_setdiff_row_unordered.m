function [valid_list, valid_idx]= fun_setdiff_row_unordered(A, B)
% fun_graph_remove_invalid_ind_pair assumes A and B are unordered list of
% data. A and B are sorted first, then all the rows in A that also appear
% in B are removed. Noticest that the origianl setdiff treate row of same
% elements but different order differently. 
% Input:
%   A: N-by-M numerical array
%   B: N'-by-M numerical array
% Output: 
%   valid_list: N''-by-M numerical array, with all the rows that also
%   appear in B removed. 
%   valid_idx: valid_list = sort(A)(valid_idx, :)


num_pair = size(A, 1);
if isempty(B)
    valid_list = A;
    valid_idx = 1 : num_pair;
    return;
end
if iscolumn(B)
    B = B';
end
sort_A = sort(A, 2);
sort_B = sort(B, 2);
[valid_list, valid_idx] = setdiff(sort_A, sort_B, 'rows', 'stable');
end
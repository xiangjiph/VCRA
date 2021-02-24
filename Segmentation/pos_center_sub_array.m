function output = pos_center_sub_array(array_size, sub_array_size, return_option)
% POS_CENTER_SUB_ARRAY returns the indice position of the subarray at the
% center of the array.
% Input:
%   array_size: the size of the array;
%   sub_array_size: the size of the sub-array
%   return_option: 
%       mmm: the starting position of the array(smallest indice value);
%       mmmxxx: the starting and end position of the array;
%       mmmlll: the starting position and the lengh of the array
%       mxmxmx: [min_i, max_i, min_j, max_j, min_k, max_k]


if nargin < 3
    return_option = 'mmmxxx';
end
start_ind = ceil((array_size - sub_array_size)/2) + 1;
end_ind = start_ind + sub_array_size - 1;

switch return_option
    case 'mmm'
        output = start_ind;
    case 'mmmxxx'
        output = [start_ind, end_ind];
    case 'mmmlll'
        output = [start_ind, sub_array_size];
    case 'mxmxmx'
        output = [start_ind;end_ind];
        output = output(:);
end


end

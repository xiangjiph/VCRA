function in_bin_stat = fun_analysis_get_basic_statistics_in_bins(data, bin_ind)
% fun_analysis_get_basic_statistics_in_bins compute the statistics of the
% data specified by the bin indices 
% Input: 
%   data: numerical vector
%   bin_ind: cell array, each cell specify the data indices to be used
%
% Output: 
%   in_bin_stat: strcutre array, each structure is defined in
%   fun_analysis_get_basic_statistics
%
% Implemented by Xiang Ji on 07/26/2019

num_val_in_bin = cellfun(@numel, bin_ind);
for iter_bin = numel(num_val_in_bin) : -1 : 1
    tmp_ind = bin_ind{iter_bin};
    tmp_data = data(tmp_ind);
    in_bin_stat(iter_bin) = fun_analysis_get_basic_statistics(tmp_data);
end
end
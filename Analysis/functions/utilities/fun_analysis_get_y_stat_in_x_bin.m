function [y_str, varargout]= fun_analysis_get_y_stat_in_x_bin(x, y, x_edge)

if isrow(x_edge)
    x_edge = x_edge';
end
is_valid_Q = isfinite(x) & isfinite(y);
x = x(is_valid_Q);
y = y(is_valid_Q);

bin_ind = fun_bin_data_to_idx_list_by_edges(x, x_edge);
num_bin = numel(bin_ind);
% Initialization
prctile_edge = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95];
mum_prctile = numel(prctile_edge);

[y_str.y_median, y_str.y_mean, y_str.y_std, y_str.y_max, y_str.y_min, y_str.y_count, ...
    ] = deal(nan(num_bin, 1));
y_str.prctile_edge = prctile_edge;
y_str.y_prctile = nan(mum_prctile, num_bin);
for iter_bin = 1 : num_bin
    tmp_list_ind = bin_ind{iter_bin}; 
    if ~isempty(tmp_list_ind)
        tmp_data = sort(y(tmp_list_ind), 'ascend');
        tmp_num_data = numel(tmp_data);
        tmp_prctile_ind = ceil(tmp_num_data .* prctile_edge);
        y_str.y_count(iter_bin) = tmp_num_data;
        y_str.y_prctile(:, iter_bin) = tmp_data(tmp_prctile_ind);
        y_str.y_mean(iter_bin) = mean(tmp_data);
        y_str.y_std(iter_bin) = std(tmp_data);
        y_str.y_min(iter_bin) = tmp_data(1);
        y_str.y_max(iter_bin) = tmp_data(end);
        if tmp_num_data == 1
            y_str.y_median(iter_bin) = tmp_data;
        elseif mod(tmp_num_data, 2) == 1
            y_str.y_median(iter_bin) = tmp_data((tmp_num_data+1)/2);
        else
            y_str.y_median(iter_bin) = (tmp_data((tmp_num_data/2) + 1) + ...
                tmp_data((tmp_num_data)/2)) / 2;
        end
    end
end
y_str.y_prctile = y_str.y_prctile.';
y_str.x_bin_val = movmean(x_edge, 2, 'Endpoints', 'discard');
y_str.x_bin_edge = x_edge;
if nargout > 1
    varargout{1} = x; 
    varargout{2} = y;
end
end
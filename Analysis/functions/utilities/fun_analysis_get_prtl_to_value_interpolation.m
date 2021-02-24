function prctile_itp = fun_analysis_get_prtl_to_value_interpolation(basic_stat_str)

int_x = basic_stat_str.hist_cdf;
int_y = basic_stat_str.hist_edge(2:end);
is_valid_pair_Q = isfinite(int_x) & isfinite(int_y);
int_x = int_x(is_valid_pair_Q);
int_y = int_y(is_valid_pair_Q);
[int_x, tmp_unique_idx] = unique(int_x, 'sorted');
int_y = int_y(tmp_unique_idx);
if ~isempty(int_y)
    prctile_itp = griddedInterpolant(int_x, int_y, 'linear', 'nearest');
else
    warning('Empty interpolation value. Return empty');
    prctile_itp = [];
end
end
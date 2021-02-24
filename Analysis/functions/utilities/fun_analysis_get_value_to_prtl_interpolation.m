function prctile_itp = fun_analysis_get_value_to_prtl_interpolation(basic_stat_str)

int_x = basic_stat_str.hist_edge(2:end);
int_y = basic_stat_str.hist_cdf;

[int_x, tmp_unique_idx] = unique(int_x, 'sorted');
int_y = int_y(tmp_unique_idx);
prctile_itp = griddedInterpolant(int_x, int_y, 'linear', 'nearest');
end
function basic_stat_str = fun_analysis_sc_basic_stat_str(basic_stat_str, scale_factor)

basic_stat_str = fun_analysis_sc_structure_fields(basic_stat_str, ...
     {'min', 'max', 'mean', 'median', 'sum', 'std', 'range', 'hist_bin_val', 'hist_edge', 'hist_bin_size'}, ...
     scale_factor);

if isfield(basic_stat_str, 'prtl2val_itp')
    basic_stat_str.prtl2val_itp.Values = basic_stat_str.prtl2val_itp.Values .* scale_factor;
end

if isfield(basic_stat_str, 'val2ptrl_itp')
    basic_stat_str.val2ptrl_itp.GridVectors = {basic_stat_str.val2ptrl_itp.GridVectors{:} .* scale_factor};
end

end
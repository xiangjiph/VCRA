function txt = fun_analysis_basic_stat_str_to_string(basic_stat_str)


txt = sprintf('Data size:%d\nMin:%.2e\nMax:%.2e\nMean:%.2e\nMedian:%.2e\nSTD:%.2e', ...
    basic_stat_str.num_data, basic_stat_str.min, basic_stat_str.max, ...
    basic_stat_str.mean, basic_stat_str.median, basic_stat_str.std);
end
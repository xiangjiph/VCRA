function hist_str = fun_get_stat_from_histogram_hdl(hist_hdl, field_name)

if nargin < 2
    field_name = {'Values', 'BinEdges', 'Normalization'};
end

hist_str = struct;
for iter_field = 1 : numel(field_name)
    hist_str.(field_name{iter_field}) = hist_hdl.(field_name{iter_field});
end
hist_str.NumData = numel(hist_hdl.Data);
end
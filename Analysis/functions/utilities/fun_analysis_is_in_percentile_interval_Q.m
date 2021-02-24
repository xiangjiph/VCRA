function is_in_interval_Q = fun_analysis_is_in_percentile_interval_Q(data, ptl_low, ptl_high, preselected_Q)

if nargin < 4
    preselected_Q = true(size(data));
end
assert(ptl_low <= ptl_high)

is_finite_Q = isfinite(data);
if ~all(is_finite_Q)
    warning('Input data contains non-finite number. Remvoe non-finite elements before computing percentile');
end

prctile_range = prctile(data(preselected_Q & is_finite_Q), [ptl_low, ptl_high]);
is_in_interval_Q = (data >= prctile_range(1)) & (data <= prctile_range(2)) & ...
    is_finite_Q & preselected_Q;
end
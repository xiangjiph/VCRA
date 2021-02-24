function normalized_diff = fun_normalized_difference(data_1, data_2, abs_Q)

if nargin < 3
    abs_Q = false;
end
normalized_diff = 2 .* (data_1 - data_2) ./ (data_1 + data_2);

if abs_Q
    normalized_diff = abs(normalized_diff);
end
end
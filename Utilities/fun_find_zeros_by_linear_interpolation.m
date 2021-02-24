function x = fun_find_zeros_by_linear_interpolation(x, y)

if ~issorted(x, 'ascend')
    [x, si] = sort(x, 'ascend');
    y = y(si);
end

num_x = numel(x);
assert(num_x == numel(y));

i = 1;
while i < num_x
    if y(i) * y(i+1) < 0
        x(i) = x(i) - (x(i+1) - x(i)) / (y(i+1) - y(i)) * y(i);
    elseif y(i) ~= 0
        x(i) = nan;
    end
    i = i + 1;
end
if y(i) ~= 0
    x(i) = nan;
end
x = x(~isnan(x));
end
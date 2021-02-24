function [bin_count, new_edge] = fun_analysis_rebin_histcount_edge(count, edge, new_edge, spacing_method)


if nargin < 4
    spacing_method = 'linear';
end

if iscolumnvector(count)
    count = count .';
end
count_n = count .* diff(edge);
ori_cdf = griddedInterpolant(edge, [0, cumsum(count_n)], 'linear', 'linear');
% [~, edge] = fun_analysis_select_histcount_edge_by_percentile(count, edge, [1e-3, 1 - 1e-3]);
if isscalar(new_edge)
    edge_min = edge(1);
    edge_max = edge(end);
    switch spacing_method
        case 'linear'
            new_edge = linspace(edge_min, edge_max, new_edge + 1);
        case 'log'
            if edge_min == 0
                edge_exp_space = log10(edge(end) / edge(2)) / new_edge;
                edge_min = edge(2) / 10^(edge_exp_space);
            end
            new_edge = 10 .^ (linspace(log10(edge_min), log10(edge_max), new_edge + 1));
        otherwise
            error('Unknown spacing method');
    end
end
% Interpolate the counts and resample the histogram
cdf_value = ori_cdf(new_edge);
bin_count = diff(cdf_value);
bin_count = bin_count ./ diff(new_edge);
if any(bin_count < 0)
    warning('Exist bin count interpolation value less than 0. Set the minimum bin count to 0');
    bin_count = max(bin_count, 0);
end
end
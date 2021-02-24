function [count, edge] = fun_analysis_select_histcount_edge_by_percentile(count, edge, ptl_range, inputType)

if nargin < 3
    ptl_range = [1e-3, 1 - 1e-3];
    inputType = 'count';
elseif nargin < 4
    inputType = 'count';
end
assert(numel(ptl_range) == 2 && ptl_range(2) > ptl_range(1));

switch inputType
    case {'count', 'pdf'}
        count_n = count ./ sum(count, 'omitnan');
        count_ncs = cumsum(count_n);
        count_ncs = [0, count_ncs];
        
        ptl_low_ind = find(count_ncs <= ptl_range(1), 1, 'last');
        ptl_high_ind = find(count_ncs >= ptl_range(2), 1, 'first');
    case 'cdf'
        count_ncs = [0, count];
        ptl_low_ind = find(count_ncs <= ptl_range(1), 1, 'last');
        ptl_high_ind = find(count_ncs >= ptl_range(2), 1, 'first');
    otherwise 
        error('To be implemented');    
end
if isempty(ptl_high_ind)
    ptl_high_ind = inf;
end
count = count(ptl_low_ind : min(numel(count), ptl_high_ind));

edge = edge(ptl_low_ind : min(ptl_high_ind + 1, numel(edge)));
end
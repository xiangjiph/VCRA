function [is_outlier_Q, varargout]= fun_analysis_is_outlier_by_percentile(data, whisker)

if nargin < 2
    whisker = 1.5;
end
if isempty(data)
    return;
end
data_ptrl = prctile(data, [25, 50, 75]);
half_width = data_ptrl(3) - data_ptrl(1);
data_up_lim = data_ptrl(3) + whisker * half_width;
data_low_lim = data_ptrl(1) - whisker * half_width;
is_outlier_Q = (data < data_low_lim) | (data > data_up_lim);

if nargout > 1
    assert(nargout == 3)
    varargout{1} = data_low_lim;
    varargout{2} = data_up_lim;
end
% Note that the limit estimated by percentile could be out of range of the
% original data.
end
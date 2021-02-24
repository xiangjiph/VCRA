function th = fun_select_threshold_by_max_avg_cc_volume(data, option)


% Default values
if nargin < 2
    option = struct;
end

if ~isfield(option, 'min_th')
    option.min_th = min(data(:));
end

if ~isfield(option, 'max_th')
    option.max_th = max(data(:));
end

if ~isfield(option, 'num_scan_step')
    option.num_scan_step = 10;
end

if ~isfield(option, 'th_list')
    option.th_list = option.min_th: (option.max_th - option.min_th)/(option.num_scan_step - 1) : option.max_th;
end

if ~isfield(option, 'return_first_max')
    option.return_first_max = true;
end
if ~isfield(option, 'show_score_plot')
    option.show_score_plot = true;
end

% bwconncomp does not apply to gpuArray
if isa(data, 'gpuArray')
    data = gather(data);
end
% Scan
avg_cc_vol = zeros(size(option.th_list));
for th_idx = 1 : length(option.th_list)
    cc = bwconncomp(data > option.th_list(th_idx));
    avg_cc_vol(th_idx) = mean(cellfun(@length, cc.PixelIdxList));
end

if option.show_score_plot
    plot(option.th_list, avg_cc_vol);
    xlabel('Threshold');
    ylabel('Average Size of Connected Components');
end

th = option.th_list(avg_cc_vol == max(avg_cc_vol));
if ~isscalar(th)
    disp('Multiple threshold value, return the first one by default')
    disp(th)
    if option.return_first_max
        th = th(1);
    end
end

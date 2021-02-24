function [th, th_info] = fun_select_threshold_by_histogram_local_min(inputArray, options)
% This function only works for integer array
if nargin < 2
    options = struct;
end

if ~isfield(options, 'brain_volume_ratio')
    options.brain_volume_ratio = 1;
end


if ~isfield(options, 'th_max')
    options.th_max = max(inputArray(:));
end
if ~isfield(options, 'th_min')
    options.th_min = min(inputArray(:));
end

if ~isfield(options, 'min_peak_distance')
    options.min_peak_distance = 8;
end

options.num_bins = options.th_max - options.th_min + 2;
[options.counts, options.edges] = histcounts(inputArray,options.num_bins, 'BinLimits', [options.th_min-0.5 ,options.th_max+0.5], 'BinMethod', 'integers');
if isa(inputArray, 'gpuArray')
    [options.counts, options.edges] = gather(options.counts, options.edges);
    block_max_value = single(intmax(classUnderlying(inputArray)));
%     inputArray_class = classUnderlying(inputArray);
else
    block_max_value = single(intmax(class(inputArray)));
%     inputArray_class = class(inputArray);
end




block_length = numel(inputArray);


options.th_list = movmean(options.edges,2, 'Endpoints', 'discard');
options.num_th = length(options.th_list);
options.min_peak_height = 0.1 * options.brain_volume_ratio/options.num_th;
options.vol_ratio = options.counts./block_length;
% Find the largest peak
[options.peak.val, options.peak.pos, options.peak.width, options.peak.prominences] = ...
    findpeaks(options.vol_ratio, 'SortStr', 'descend', 'MinPeakHeight',options.min_peak_height, ...
    'MinPeakDistance',options.min_peak_distance, 'WidthReference','halfheight');
% disp(options)
if isempty(options.peak.val)
    disp('No peak found!');
    options.failQ = true;
    options.search_th_min = options.th_min;
    options.th = options.th_min;
else
    if all(options.peak.pos(1) <= options.peak.pos)
        % If the largest peak is on the left of all the rest of the peak,
        % use it
%         disp('The highest peak is on the left of the rest of the peak. Use it to estimate the background histogram');
        options.peak.bg_peak_pos = options.peak.pos(1);
        options.peak.bg_peak_width = options.peak.width(1);
    else
        disp('The highest peak is not on the left of the rest of the peak. Choose the peak that locates closest to the center(Normalized by peak width)');
        tmp_peak_th_values = options.th_list(options.peak.pos);
        tmp_relative_dis_to_center = abs(tmp_peak_th_values - block_max_value/2)./options.peak.width;
        [~, tmp_peak_idx] = min(tmp_relative_dis_to_center);
        options.peak.bg_peak_pos = options.peak.pos(tmp_peak_idx);
        options.peak.bg_peak_width = options.peak.width(tmp_peak_idx);
    end
    
    
    options.peak.bg_peak_tail = options.peak.bg_peak_pos + round(options.peak.bg_peak_width);
    if options.peak.bg_peak_tail > length(options.th_list)
        options.peak.bg_peak_tail = options.peak.bg_peak_pos;
    end

    options.search_th_min = options.th_list(options.peak.bg_peak_tail);
    options.inv_th_list = options.th_list(options.peak.bg_peak_tail: end);
    options.inv_vol_ratio = - options.vol_ratio(options.peak.bg_peak_tail: end);
    % Find the local peak of the inverted volume ratio that is cloest
    % to the threshold
    [options.inv_peak.val, options.inv_peak.pos, options.inv_peak.width, ...
        options.inv_peak.prominences] = findpeaks(options.inv_vol_ratio);
    if ~isempty(options.inv_peak.val)
        tmp_relative_dis_to_center = options.inv_peak.pos ./ options.inv_peak.width;
        [~, tmp_peak_idx] = min(tmp_relative_dis_to_center);
        options.th = options.inv_th_list(options.inv_peak.pos(tmp_peak_idx));
    else
        disp('No peak for the inversed histogram. Estimated as the tail of the peak');
        tmp_est_tail_pos = options.peak.bg_peak_pos + round(2 * options.peak.bg_peak_width);
        if tmp_est_tail_pos >= options.num_th
            options.th = options.th_list(options.peak.bg_peak_pos + round(options.peak.bg_peak_width));
        else
            options.th = options.th_list(tmp_est_tail_pos);
        end
    end
end
th = options.th;
if nargout > 1
    th_info = options;
end

    

end

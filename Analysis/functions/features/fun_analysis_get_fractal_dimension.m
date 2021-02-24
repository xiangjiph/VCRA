function stat_str = fun_analysis_get_fractal_dimension(mask, cutoff_length, visQ)
% fun_analysis_get_fractal_dimension computes the box-counting dimension of
% the mask 
% Input: 
%   mask: 3D logical array
%   cutoff_length: Fit the data greater or smaller than this length
%   seperately. If is empty, use all the data for fitting
%   visQ: visualize the data and fitting. 
% Output:
%   stat_str: structure with fields
%
% Implemented by Xiang Ji on 02/26/2019
if nargin < 2
    cutoff_length = [];
    visQ = false;
elseif nargin < 3
    visQ = false;
end
mask_size = size(mask);
mask_size_min = min(mask_size);

exp_step = 0.2;
exp_start = 0;
exp_end = floor(log2(mask_size_min));
bbox_size_list = unique(round(2 .^ (exp_start:exp_step:exp_end)));

num_bbox = numel(bbox_size_list);
bbox_count = zeros(num_bbox,1); 
vol_ratio = zeros(num_bbox,1);
for iter_bbox = 1 : num_bbox
    if bbox_size_list(iter_bbox) == 1
        bbox_count(iter_bbox) = nnz(mask); 
        vol_ratio(iter_bbox) =  bbox_count(iter_bbox)/ numel(mask);
    else
        tmp_bbox_size = ones(1,3) .* bbox_size_list(iter_bbox);
        tmp_im = fun_downsample_by_block_operation(mask, @max, tmp_bbox_size, 'true');
        bbox_count(iter_bbox) = nnz(tmp_im); 
        vol_ratio(iter_bbox) = bbox_count(iter_bbox) / numel(tmp_im);
    end
end
cutoff_length = min(min(bbox_size_list), cutoff_length);

stat_str = struct;
stat_str.bbox_size_list = bbox_size_list;
stat_str.bbox_count_list = bbox_count;
stat_str.fit_cutoff_length = cutoff_length;
if ~isempty(cutoff_length)
    small_scall_Q = bbox_size_list < cutoff_length;
    % Use linear regression to extract the box-counting dimension of the vessel
    % network 
    vessel_dim_linear_fit_1 = fitlm(log(1./bbox_size_list(small_scall_Q)), log(bbox_count(small_scall_Q)), 'linear');
    fit_y_1 = exp(vessel_dim_linear_fit_1.Coefficients.Estimate(1)) * (bbox_size_list(small_scall_Q)).^(- vessel_dim_linear_fit_1.Coefficients.Estimate(2));
    vessel_dim_linear_fit_2 = fitlm(log(1./bbox_size_list(~small_scall_Q)), log(bbox_count(~small_scall_Q)), 'linear');
    fit_y_2 = exp(vessel_dim_linear_fit_2.Coefficients.Estimate(1)) * (bbox_size_list(~small_scall_Q)).^(- vessel_dim_linear_fit_2.Coefficients.Estimate(2));
    if visQ
        figure;
        scatter(bbox_size_list, bbox_count, 'o');
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
        hold on
        plot(bbox_size_list(small_scall_Q), fit_y_1, 'LineWidth', 2)
        hold on
        plot(bbox_size_list(~small_scall_Q), fit_y_2, 'LineWidth', 2)
        grid on
        xlabel('Box size/\mum');
        ylabel('Number of covering boxes');
        legend('Data', sprintf('Box-counting dimension = %0.2f+/-%0.2f', vessel_dim_linear_fit_1.Coefficients.Estimate(2), vessel_dim_linear_fit_1.Coefficients.SE(2)), ...
            sprintf('Box-counting dimension = %0.2f+/-%0.2f', vessel_dim_linear_fit_2.Coefficients.Estimate(2), vessel_dim_linear_fit_2.Coefficients.SE(2)));
    end
    stat_str.lmfit_1 = vessel_dim_linear_fit_1;
    stat_str.fit_y_1 = fit_y_1;
    stat_str.dimension_1 = vessel_dim_linear_fit_1.Coefficients.Estimate(2);
    stat_str.lmfit_2 = vessel_dim_linear_fit_2;
    stat_str.fit_y_2 = fit_y_2;
    stat_str.dimension_2 = vessel_dim_linear_fit_2.Coefficients.Estimate(2);
else
    vessel_dim_linear_fit_1 = fitlm(log(1./bbox_size_list), log(bbox_count), 'linear');
    fit_y_1 = exp(vessel_dim_linear_fit_1.Coefficients.Estimate(1)) * (bbox_size_list).^(- vessel_dim_linear_fit_1.Coefficients.Estimate(2));
    if visQ
        figure;
        scatter(bbox_size_list, bbox_count, 'o');
        set(gca, 'XScale', 'log');
        set(gca, 'YScale', 'log');
        hold on
        plot(bbox_size_list(small_scall_Q), fit_y_1, 'LineWidth', 2)
        grid on
        xlabel('Box size/\mum');
        ylabel('Number of covering boxes');
        legend('Data', sprintf('Box-counting dimension = %0.2f+/-%0.2f', vessel_dim_linear_fit_1.Coefficients.Estimate(2), vessel_dim_linear_fit_1.Coefficients.SE(2)));
    end
    stat_str.lmfit = vessel_dim_linear_fit_1;
    stat_str.fit_y = fit_y_1;
    stat_str.dimension = vessel_dim_linear_fit_1.Coefficients.Estimate(2);
end


end
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML20200201';
psf_data = DataManager.load_data(DataManager.fp_metadata_file(dataset_name, ...
    stack, 'PSF_estimation'));
%%
wb_psf_est_result = psf_data.result;
wb_psf_est_result = cat(2, wb_psf_est_result{:});
[wb_psf_stat_xy.mean, wb_psf_stat_xy.median, wb_psf_stat_xy.std] = deal(nan(size(wb_psf_est_result)));
[wb_psf_stat_z.mean, wb_psf_stat_z.median, wb_psf_stat_z.std] = deal(nan(size(wb_psf_est_result)));
[wb_psf_xy_stat_cell, wb_psf_z_stat_cell] = deal(cell(size(wb_psf_est_result)));
for iter_region = 1 : numel(wb_psf_est_result)
    tmp_data = wb_psf_est_result{iter_region};
    tmp_vis_Q = (tmp_data.refine_r_est <= 3 & tmp_data.best_fit_corr > 0.99 & ...
        abs(tmp_data.link_ori_vec(:, 3)) < sqrt(2)/2);
    tmp_xy_stat = fun_analysis_get_basic_statistics(tmp_data.best_fit_PSF_FWHM(tmp_vis_Q, 1));
    wb_psf_stat_xy.mean(iter_region) = tmp_xy_stat.mean;
    wb_psf_stat_xy.median(iter_region) = tmp_xy_stat.median;
    wb_psf_stat_xy.std(iter_region) = tmp_xy_stat.std;
    wb_psf_xy_stat_cell{iter_region} = tmp_xy_stat;
    
    tmp_z_stat = fun_analysis_get_basic_statistics(tmp_data.best_fit_PSF_FWHM(tmp_vis_Q, 3));
    wb_psf_stat_z.mean(iter_region) = tmp_z_stat.mean;
    wb_psf_stat_z.median(iter_region) = tmp_z_stat.median;
    wb_psf_stat_z.std(iter_region) = tmp_z_stat.std;
    wb_psf_z_stat_cell{iter_region} = wb_psf_stat_z;
end
fprintf('FWHM xy: %.2f +/- %.2f\n', mean(wb_psf_stat_xy.mean), std(wb_psf_stat_xy.mean));
fprintf('FWHM z: %.2f +/- %.2f\n', mean(wb_psf_stat_z.mean), std(wb_psf_stat_z.mean));
%% For each anaotimical region, plot one global hist 
num_regions = numel(psf_data.result);
region_name = psf_data.psf_est_setting(1, :);

fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 1];
ax_hdl_1 = subplot(1,2,1);
ax_hdl_2 = subplot(1,2,2);
tmp_num_valid_points = nan(num_regions, 1);
leg_string_cell = cell(num_regions, 1);
for iter_region = 1 : num_regions
    tmp_data = psf_data.result{iter_region};
    tmp_PSF = cellfun(@(x) x.best_fit_PSF_FWHM, tmp_data, 'UniformOutput', false);
    tmp_PSF = cat(1, tmp_PSF{:});
    
    tmp_est_r = cellfun(@(x) x.refine_r_est, tmp_data, 'UniformOutput', false);
    tmp_est_r = cat(1, tmp_est_r{:});
    
    tmp_fit_corr = cellfun(@(x) x.best_fit_corr, tmp_data, 'UniformOutput', false);
    tmp_fit_corr = cat(1, tmp_fit_corr{:});
    
    tmp_segment_ori = cellfun(@(x) x.link_ori_vec(:, 3), tmp_data, 'UniformOutput', false);
    tmp_segment_ori = cat(1, tmp_segment_ori{:});
    
    tmp_vis_Q = (tmp_est_r <= 3 & tmp_fit_corr > 0.99 & abs(tmp_segment_ori) < sqrt(2)/2);
    tmp_num_valid_points(iter_region) = nnz(tmp_vis_Q);
    tmp_xy_fwhm = tmp_PSF(tmp_vis_Q, 1);    
    tmp_z_fwhm = tmp_PSF(tmp_vis_Q, 3);
    
    tmp_x_bin_val = 0.00 : 0.11 : 2;
    tmp_pdf = histcounts(tmp_xy_fwhm, tmp_x_bin_val, 'Normalization', 'pdf');
    plot(ax_hdl_1, movmean(tmp_x_bin_val, 2, 'Endpoints', 'discard'), ...
        tmp_pdf, 'LineWidth', 2);
    hold(ax_hdl_1, 'on');  
    
    tmp_z_bin_val = 0 : 0.55 : 10;
    tmp_pdf = histcounts(tmp_z_fwhm, tmp_z_bin_val, 'Normalization', 'pdf');
    plot(ax_hdl_2, movmean(tmp_z_bin_val, 2, 'Endpoints', 'discard'), ...
        tmp_pdf, 'LineWidth', 2);
    hold(ax_hdl_2, 'on');         
    
    leg_string_cell{iter_region} = sprintf('%s: %d', region_name{iter_region}, ...
        tmp_num_valid_points(iter_region));    
end

leg_hdl_1 = legend(ax_hdl_1, leg_string_cell);
leg_hdl_2 = legend(ax_hdl_2, region_name, 'Location', 'northwest');
ax_hdl_1.FontSize = 14;
ax_hdl_2.FontSize = 14;
ax_hdl_1.Box = 'off';
ax_hdl_2.Box = 'off';
ax_hdl_1.XLabel.String = 'FWHM_{xy} (\mum)';
ax_hdl_2.XLabel.String = 'FWHM_z (\mum)';
ax_hdl_1.YLabel.String = 'PDF';
ax_hdl_2.YLabel.String = 'PDF';
leg_hdl_1.FontSize = 10;
fig_fp = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), ...
    'Paper', 'PSF_estimation_regional_overall_pdf.png');
fun_print_image_in_several_formats(fig_hdl, fig_fp);
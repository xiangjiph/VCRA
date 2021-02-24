%% Random orientation isotropy 
% num_trial = 50000;
% num_sample = 200;
% vec_sum = zeros(3, num_trial);
% clearvars('isotropy_str');
% for iter_trial = 1 : num_trial
%     rand_ori = randn(num_sample, 3);
%     rand_ori = rand_ori ./ vecnorm(rand_ori')';
%     isotropy_str(iter_trial) = fun_analysis_get_link_anisotropy(rand_ori, true);
% end
% figure;
% subplot(1,2,1)
% scatter3(rand_ori(:, 2), rand_ori(:, 1), rand_ori(:, 3))
% daspect([1,1,1]);
% subplot(1,2,2)
% histogram([isotropy_str.svd_min2max], 'Normalization', 'pdf');
% title(sprintf('Isotropy = %.2f \\pm %.2f', mean([isotropy_str.svd_min2max]), std([isotropy_str.svd_min2max])), 'FontSize', 14);
% xlabel('Isotropy', 'FontSize', 14);
% ylabel('Probability', 'FontSize', 14);
%%
DataManager = FileManager;
num_sample_list = [3 : 9, 10 : 5 : 100, 110 : 10 : 1000, 1020 : 20 : 1500];
num_trial = 5e4;
[min2max_list, svd_list, FA_list] = deal(cell(size(num_sample_list)));
parfor iter_sample = 1 : numel(num_sample_list)
    tmp_tic = tic;
    tmp_num_sample = num_sample_list(iter_sample);
    tmp_min2max_list = nan(num_trial, 1);
    tmp_svd = nan(3, num_trial);
    tmp_fa = nan(num_trial, 1);
    for iter_trial = 1 : num_trial
        rand_ori = randn(tmp_num_sample, 3);
        rand_ori = rand_ori ./ vecnorm(rand_ori')';
        tmp_str = fun_analysis_get_link_anisotropy(rand_ori, true);
        tmp_min2max_list(iter_trial) = tmp_str.svd_min2max;
        tmp_svd(:, iter_trial) = tmp_str.svd_value_ratio';
        tmp_fa(iter_trial) = tmp_str.fractional_anisotropy;
    end
    min2max_list{iter_sample} = tmp_min2max_list;
    FA_list{iter_sample} = tmp_fa;
    svd_list{iter_sample} = tmp_svd';
    fprintf('Finish computation for sample with %d linkes. Elapse time is %d seconds.\n', tmp_num_sample, toc(tmp_tic));
end
%%
min2max_mean = cellfun(@mean, min2max_list);
min2max_std = cellfun(@std, min2max_list);

fa_mean = cellfun(@mean, FA_list);
fa_std = cellfun(@std, FA_list);

svd_mean = cellfun(@mean, svd_list, 'UniformOutput', false);
svd_mean = cat(1, svd_mean{:});

svd_std = cellfun(@std, svd_list, 'UniformOutput', false);
svd_std = cat(1, svd_std{:});
% Generate the interpolation for average isotropy and its standard
% deviation 
interpolation_method = 'spline';
extropolation_method = 'linear';

isotropy_str = struct;
isotropy_str.min2max.mean = [];
isotropy_str.min2max.mean = griddedInterpolant(num_sample_list, min2max_mean, interpolation_method, extropolation_method);
isotropy_str.min2max.std = [];
isotropy_str.min2max.std = griddedInterpolant(num_sample_list, min2max_std, interpolation_method, extropolation_method);

isotropy_str.svd_1.mean = [];
isotropy_str.svd_1.mean = griddedInterpolant(num_sample_list, svd_mean(:, 1), interpolation_method, extropolation_method);
isotropy_str.svd_1.std = [];
isotropy_str.svd_1.std = griddedInterpolant(num_sample_list, svd_std(:, 1), interpolation_method, extropolation_method);

isotropy_str.svd_2.mean = [];
isotropy_str.svd_2.mean = griddedInterpolant(num_sample_list, svd_mean(:, 2), interpolation_method, extropolation_method);
isotropy_str.svd_2.std = [];
isotropy_str.svd_2.std = griddedInterpolant(num_sample_list, svd_std(:, 2), interpolation_method, extropolation_method);

isotropy_str.svd_3.mean = [];
isotropy_str.svd_3.mean = griddedInterpolant(num_sample_list, svd_mean(:, 3), interpolation_method, extropolation_method);
isotropy_str.svd_3.std = [];
isotropy_str.svd_3.std = griddedInterpolant(num_sample_list, svd_std(:, 3), interpolation_method, extropolation_method);

isotropy_str.fa.mean = [];
isotropy_str.fa.mean = griddedInterpolant(num_sample_list, fa_mean, interpolation_method, extropolation_method);
isotropy_str.fa.std = [];
isotropy_str.fa.std = griddedInterpolant(num_sample_list, fa_std, interpolation_method, extropolation_method);

save('./Metadata/uni_ori_isotropy.mat', '-struct', 'isotropy_str');
%% Visualization
fig_hdl = figure;
fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [3, 6];
ax1 = subplot(3,1,1);
errorbar(ax1, num_sample_list, min2max_mean, min2max_std, 'LineWidth', 2);
grid(ax1, 'on');
ax1.XLabel.String = 'Number of points';
ax1.YLabel.String = 'Average min-to-max singluar value ratio';
ax1.FontSize = 14;
ax1.FontWeight = 'bold';
set(ax1, 'XScale', 'log');

ax2 = subplot(3,1,3);
errorbar(ax2, num_sample_list, fa_mean, fa_std, 'LineWidth', 2);
grid(ax2, 'on');
ax2.XLabel.String = 'Number of points';
ax2.YLabel.String = 'Average fractional anisotropy';
ax2.FontSize = 14;
ax2.FontWeight = 'bold';
set(ax2, 'XScale', 'log');

ax3 = subplot(3,1,2);
errorbar(ax3, num_sample_list, svd_mean(:, 1), svd_std(:, 1), 'LineWidth', 2);
grid(ax3, 'on');
ax3.XLabel.String = 'Number of points';
ax3.YLabel.String = 'Average largest singular value fraction';
ax3.FontSize = 14;
ax3.FontWeight = 'bold';
set(ax3, 'XScale', 'log');

fig_fp = fullfile(DataManager.AnalysisResultRootPath, 'Unweighted_isotropy_statistics_vs_num_links_log.png');
fun_print_image_in_several_formats(fig_hdl, fig_fp);
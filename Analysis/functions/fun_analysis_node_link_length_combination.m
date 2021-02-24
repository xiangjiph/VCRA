function result_str = fun_analysis_node_link_length_combination(node_l_min, node_l_med, node_l_max, link_length_dist, visibleQ)

if nargin < 5
    visibleQ = false;
end
result_str = struct;

valid_node_data_Q = ~isnan(node_l_min) & ~isnan(node_l_med) & ~isnan(node_l_max);
node_l_min = node_l_min(valid_node_data_Q);
node_l_med = node_l_med(valid_node_data_Q);
node_l_max = node_l_max(valid_node_data_Q);
node_link_comb = cat(2, node_l_min, node_l_med, node_l_max);
link_length_dist = link_length_dist(~isnan(link_length_dist));

fprintf('Start sampling link combination from the link length distribution with replacement\n');
tic
simu_link_stat = struct;

random_sample_comb = randsample(link_length_dist, numel(node_l_min) * 3, true);
random_sample_comb = reshape(random_sample_comb, [], 3);
[result_str.simu.RS_R, result_str.simu.RS_P] = corrcoef(random_sample_comb);
random_sample_sorted = sort(random_sample_comb, 2, 'ascend');
[result_str.simu.R, result_str.simu.P] = corrcoef(random_sample_sorted);

simu_link_stat.length_max = random_sample_sorted(:, 3);
simu_link_stat.length_min = random_sample_sorted(:, 1);
simu_link_stat.length_median = random_sample_sorted(:, 2);
toc

result_str.simu.stat.min = fun_analysis_get_basic_statistics(simu_link_stat.length_min);
result_str.simu.stat.median = fun_analysis_get_basic_statistics(simu_link_stat.length_median);
result_str.simu.stat.max = fun_analysis_get_basic_statistics(simu_link_stat.length_max);

warning('Node link with length equal to 1 are exlucded. A temporary fix to a bug.');
node_l_min = node_l_min(node_l_min > 1);
node_l_max = node_l_max(node_l_max > 1);
node_l_med = node_l_med(node_l_med > 1);

result_str.data.stat.min = fun_analysis_get_basic_statistics(node_l_min);
result_str.data.stat.median = fun_analysis_get_basic_statistics(node_l_med);
result_str.data.stat.max = fun_analysis_get_basic_statistics(node_l_max);
%% Visualize the histogram and compute KStest2
if visibleQ
    fprintf('Visualizing the histogram\n');
    fig_hdl = figure;
    fig_hdl.Position(3:4) = fig_hdl.Position(3:4) .* [2, 3];
    ax_1 = subplot(3,1,1);
    histogram(ax_1, node_l_min, 'Normalization', 'pdf');
    hold(ax_1, 'on');
    histogram(ax_1, simu_link_stat.length_min, 'Normalization', 'pdf');
    legend(ax_1, 'Data', 'Simulation');
    [~, p1, ~] = kstest2(node_l_min, simu_link_stat.length_min);
    ax_1.Title.String = sprintf('KSTest2 p-Value = %.2e', p1);
    ax_1.XLabel.String = 'Shortest vessel segment length (\mum)';
    ax_1.YLabel.String = 'PDF';
    ax_1.FontWeight = 'bold';
    grid(ax_1, 'on');
    legend(ax_1, sprintf('Observation\n%s', fun_analysis_basic_stat_str_to_string(result_str.data.stat.min)), ...
        sprintf('Simulation\n%s', fun_analysis_basic_stat_str_to_string(result_str.simu.stat.min)));
    
    ax_2 = subplot(3,1,2);
    histogram(ax_2, node_l_med, 'Normalization', 'pdf');
    hold(ax_2, 'on');
    histogram(ax_2, simu_link_stat.length_median, 'Normalization', 'pdf');
    legend(ax_2, 'Data', 'Simulation');
    [~, p2, ~] = kstest2(node_l_med, simu_link_stat.length_median);
    ax_2.Title.String = sprintf('KSTest2 p-Value = %.2e', p2);
    ax_2.XLabel.String = 'Median Vessel segment length (\mum)';
    ax_2.YLabel.String = 'PDF';
    ax_2.FontWeight = 'bold';
    grid(ax_2, 'on');
    legend(ax_2, sprintf('Observation\n%s', fun_analysis_basic_stat_str_to_string(result_str.data.stat.median)), ...
        sprintf('Simulation\n%s', fun_analysis_basic_stat_str_to_string(result_str.simu.stat.median)));
    
    ax_3 = subplot(3,1,3);
    histogram(ax_3, node_l_max, 'Normalization', 'pdf');
    hold(ax_3, 'on');
    histogram(ax_3, simu_link_stat.length_max, 'Normalization', 'pdf');
    legend(ax_3, 'Data', 'Simulation');
    [~, p3, ~] = kstest2(node_l_max, simu_link_stat.length_max);
    ax_3.Title.String = sprintf('KSTest2 p-Value = %.2e', p3);
    ax_3.XLabel.String = 'Longest vessel segment length (\mum)';
    ax_3.YLabel.String = 'PDF';
    ax_3.FontWeight = 'bold';
    grid(ax_3, 'on');
    legend(ax_3, sprintf('Observation\n%s', fun_analysis_basic_stat_str_to_string(result_str.data.stat.max)), ...
        sprintf('Simulation\n%s', fun_analysis_basic_stat_str_to_string(result_str.simu.stat.max)));
    
    result_str.KSTest2.min = p1;
    result_str.KSTest2.median = p2;
    result_str.KSTest2.max = p3;
    result_str.fig_hdl = fig_hdl;
    result_str.CData = getframe(fig_hdl);
end
%% Random shuffle the order of node length combination and compute the correlation matrix
randShuffled_comb = node_link_comb';
for iter_node = 1 : size(node_link_comb, 1)
    randShuffled_comb(:, iter_node) = randShuffled_comb(randperm(3), iter_node);
end
randShuffled_comb = randShuffled_comb.';
[result_str.data.RS_R, result_str.data.RS_P] = corrcoef(randShuffled_comb);

[result_str.data.R, result_str.data.P] = corrcoef(node_link_comb);
end


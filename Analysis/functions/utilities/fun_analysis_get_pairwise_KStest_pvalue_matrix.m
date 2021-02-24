function pair_pvalue_mat = fun_analysis_get_pairwise_KStest_pvalue_matrix(data_cell)

num_data = numel(data_cell);
pair_pvalue_mat = nan(num_data, num_data);

for iter_cell_1 = 1 : num_data
    tmp_data_1 = data_cell{iter_cell_1};
    for iter_cell_2 = iter_cell_1 + 1 : num_data
        tmp_data_2 = data_cell{iter_cell_2};
        [~, tmp_pvalue, ~] = kstest2(tmp_data_1, tmp_data_2);
        pair_pvalue_mat(iter_cell_1, iter_cell_2) = tmp_pvalue;
        pair_pvalue_mat(iter_cell_2, iter_cell_1) = tmp_pvalue;
    end
end
end
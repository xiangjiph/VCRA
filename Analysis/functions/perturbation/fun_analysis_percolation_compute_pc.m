function pc = fun_analysis_percolation_compute_pc(occupancy_p, gaint_cc_fraction)

% Compute the threshold by finding the occupancy with largest derivative
diff_p = occupancy_p(3 : end) - occupancy_p(1 : end-2);
diff_s = gaint_cc_fraction(3 : end) - gaint_cc_fraction(1 : end - 2);
ds_dp = diff_s ./ diff_p;
int_x = occupancy_p(2 : end - 1);
[~, max_ind] = max(ds_dp);
pc = int_x(max_ind);
end
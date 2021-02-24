function cc_dist_to_surf = fun_analysis_get_cc_distance_to_brain_surface(cc_ind, cropped_brain_mask_dt, sampleFun)

if nargin < 3
    sampleFun = @mean;
end

num_cc = numel(cc_ind);
cc_dist_to_surf = nan(num_cc, 1);
if num_cc == 0
    return;
else
   for iter_cc = 1 : num_cc
       tmp_dt = cropped_brain_mask_dt(cc_ind{iter_cc});
       cc_dist_to_surf(iter_cc) = sampleFun(tmp_dt);
   end
end
end
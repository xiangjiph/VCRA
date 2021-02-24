function pO2_cube_data = fun_analysis_sc_pO2_result(pO2_cube_data, linear_scaling_coeff, radius_scaling_coeff)

% Krogh fitting - scale the radius of capillary together. Otherwise, no way
% to rescale here actually. And therefore, the correction coefficient does
% not change. 
if isfield(pO2_cube_data, 'fit_Krogh')
    pO2_cube_data.fit_Krogh.d_max = pO2_cube_data.fit_Krogh.d_max .* linear_scaling_coeff;
    pO2_cube_data.fit_Krogh.cap_r = pO2_cube_data.fit_Krogh.cap_r .* radius_scaling_coeff;
    pO2_cube_data.fit_fun_hdl = @(x) (-1/2) * ((pO2_cube_data.fit_Krogh.d_max + pO2_cube_data.fit_Krogh.cap_r) .^ 2 .* ...
        log((x + pO2_cube_data.fit_Krogh.cap_r) ./ pO2_cube_data.fit_Krogh.cap_r) - ...
        ((x + pO2_cube_data.fit_Krogh.cap_r) .^ 2 - pO2_cube_data.fit_Krogh.cap_r ^2)/2) * linear_scaling_coeff^2;
end
% local voxel statistics
pO2_cube_data.local_dt_stat = fun_analysis_sc_basic_stat_str(pO2_cube_data.local_dt_stat, linear_scaling_coeff);

pO2_cube_data.local_pO2_stat = fun_analysis_sc_basic_stat_str(pO2_cube_data.local_pO2_stat, linear_scaling_coeff ^ 2);
% pO2_dt_hist_2: 
pO2_cube_data.pO2_dt_hist2.dt_edge = pO2_cube_data.pO2_dt_hist2.dt_edge .* linear_scaling_coeff;
pO2_cube_data.pO2_dt_hist2.pO2_edge = pO2_cube_data.pO2_dt_hist2.pO2_edge .* linear_scaling_coeff ^ 2;
% Local DT maximum
lm_field = {'pO2_lm', 'dt_lm', 'pO2_lm_stat', 'dt_lm_stat'};
for iter_field = 1 : numel(lm_field)
   tmp_field_name = lm_field{iter_field};
   tmp_subfield_list = fieldnames(pO2_cube_data.(tmp_field_name));
   for iter_subfield = 1 : numel(tmp_subfield_list)
       tmp_subfield_name = tmp_subfield_list{iter_subfield};
       if contains(tmp_subfield_name, 'dt')
           tmp_cor_coeff = linear_scaling_coeff;
       elseif contains(tmp_subfield_name, 'pO2')
           tmp_cor_coeff = linear_scaling_coeff ^ 2;
       else
           continue;
       end       
       tmp_data = pO2_cube_data.(tmp_field_name).(tmp_subfield_name);
       if iscell(tmp_data)
           pO2_cube_data.(tmp_field_name).(tmp_subfield_name) = cellfun(@(x) x .* tmp_cor_coeff, ...
               tmp_data, 'UniformOutput', false);
       elseif isnumeric(tmp_data)
           pO2_cube_data.(tmp_field_name).(tmp_subfield_name) = tmp_data .* tmp_cor_coeff;
       end
   end    
end
% pO2_cube_data.dt_lm.v = pO2_cube_data.dt_lm.v .* linear_scaling_coeff;
% pO2_cube_data.dt_lm.pO2_v = pO2_cube_data.dt_lm.pO2_v .* (linear_scaling_coeff ^ 2);
% Local pO2 minimum
% pO2_cube_data.pO2_lm.v = pO2_cube_data.pO2_lm.v .* (linear_scaling_coeff ^ 2);
% pO2_cube_data.pO2_lm.dt_v = pO2_cube_data.pO2_lm.dt_v .* linear_scaling_coeff;
% pO2 vs dt bin
pO2_cube_data.pO2_stat_in_dt_bin = fun_analysis_sc_structure_fields(pO2_cube_data.pO2_stat_in_dt_bin, ...
    {'y_median', 'y_mean', 'y_std', 'y_max', 'y_min', 'y_prctile'}, linear_scaling_coeff ^ 2);
pO2_cube_data.pO2_stat_in_dt_bin = fun_analysis_sc_structure_fields(pO2_cube_data.pO2_stat_in_dt_bin, ...
    {'x_bin_val', 'x_bin_edge'}, linear_scaling_coeff);
% paired_extrema
% pO2_cube_data.paired_extrema.dist = pO2_cube_data.paired_extrema.dist .* linear_scaling_coeff;

end
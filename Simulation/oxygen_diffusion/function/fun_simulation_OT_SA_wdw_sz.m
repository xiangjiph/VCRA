function [dt_le, pO2_le]= fun_simulation_OT_SA_wdw_sz(dt_array, pO2_array, ...
    wz_list_pxl, valid_mask)

if nargin < 4
    valid_mask = true(size(dt_array));
end

num_wz = numel(wz_list_pxl);
[dt_le, pO2_le] = deal(struct);
[dt_le.ind, dt_le.sub, dt_le.dt, dt_le.pO2, ...
    pO2_le.ind, pO2_le.sub, pO2_le.dt, pO2_le.pO2] = deal(cell(num_wz, 1));
for iter_wz = 1 : num_wz
    tmp_wz = wz_list_pxl(iter_wz);
    % Search for local maximum of DT
    tmp_extrema_info = fun_analysis_get_local_extrema_info(dt_array, ...
        tmp_wz, 'max');
    tmp_extrema_info = fun_structure_field_indexing(tmp_extrema_info, ...
        valid_mask(tmp_extrema_info.ind));
    dt_le.ind{iter_wz} = tmp_extrema_info.ind;    
    dt_le.sub{iter_wz} = tmp_extrema_info.sub;    
    dt_le.dt{iter_wz} = tmp_extrema_info.v;    
    dt_le.pO2{iter_wz} = pO2_array(tmp_extrema_info.ind);
    
    % Search for local minimum of pO2
    tmp_extrema_info = fun_analysis_get_local_extrema_info(pO2_array, ...
        tmp_wz, 'min');
    tmp_extrema_info = fun_structure_field_indexing(tmp_extrema_info, ...
        valid_mask(tmp_extrema_info.ind));
    pO2_le.ind{iter_wz} = tmp_extrema_info.ind;    
    pO2_le.sub{iter_wz} = tmp_extrema_info.sub;    
    pO2_le.pO2{iter_wz} = tmp_extrema_info.v;    
    pO2_le.dt{iter_wz} = dt_array(tmp_extrema_info.ind);
end

end
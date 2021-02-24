function extrema_info = fun_analysis_get_local_extrema_info(data_array, window_size, extrema_type)

switch extrema_type
    case {'max', 'maximum', 'maxima'}
        is_local_extrema_Q = fun_array_local_maximum(data_array, window_size);
    case {'min', 'minimum'}
        is_local_extrema_Q = fun_array_local_maximum(-data_array, window_size);
    otherwise
        error('Unknown extrema type');
end

extrema_info.ind = find(is_local_extrema_Q);
extrema_info.sub = fun_ind2sub(size(is_local_extrema_Q), extrema_info.ind);
extrema_info.v = data_array(extrema_info.ind);
end
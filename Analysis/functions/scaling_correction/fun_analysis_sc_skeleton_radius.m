function skel_str = fun_analysis_sc_skeleton_radius(skel_str, radius_scaling_factor)

if isnumeric(radius_scaling_factor)
    assert(isscalar(radius_scaling_factor) && isfinite(radius_scaling_factor), 'The radius scaling correction factor should be finite');
    skel_str.r = skel_str.r * radius_scaling_factor;
elseif isa(radius_scaling_factor, 'griddedInterpolant')
    input_r_type = class(skel_str.r);
    min_convert_r = min(radius_scaling_factor.GridVectors{1});
    convert_Q = (skel_str.r >= min_convert_r);
    skel_str.r(convert_Q) = cast(radius_scaling_factor(skel_str.r(convert_Q)), input_r_type);    
end
end
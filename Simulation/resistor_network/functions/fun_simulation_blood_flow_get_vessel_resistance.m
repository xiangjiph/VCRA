function r = fun_simulation_blood_flow_get_vessel_resistance(vessel_r, vessel_l, input_unit)
% fun_simulation_get_vessel_resistance computes the vessel resistance using
% the formula from the Online Method of Blinder et al 2013 Nature neuroscience
% Input: 
%   vessel_r: numerical array, radius of vessel in micron 
%   vessel_l: numerical array of the same size as vessel_r, length of
%   vessels in micron. 
%   r: numerical array, resistance of the vessels in SI unit. 
% 
if nargin < 3
    input_unit = 'm';
end
viscosity_water = 0.6913e-3; % Pa s. Source: https://wiki.anton-paar.com/en/water/
% viscosity_water = 1e-3;
exp_coeff_1 = 14.285714 * 10^-6;
exp_coeff_2 = 0.3511235955056 * 10^-6;
% Convert the radius and length to meter
switch input_unit
    case 'm'
        % Do nothing
    case {'um', 'micron'}
        vessel_r = vessel_r ./ 10^6;
        vessel_l = vessel_l ./ 10^6;
    otherwise
        error('Unrecognized unit');
end
r = (32 * viscosity_water / pi) .* (vessel_l ./ vessel_r .^ 4) .* ...
    (1 - 0.8625 .* exp( - vessel_r ./ exp_coeff_1) + 27.5 .* exp(- vessel_r ./ exp_coeff_2)); 
end
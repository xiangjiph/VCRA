function vesselness = fun_gpu_Frangi_vesselness(A11,A12,A13,A22,A23,A33)
%% Compute eigenvalues
c = 2 * pi/3;
precision_error = eps('double');
p = A12 * A12 + A13 * A13 + A23 * A23;
eig1 = A11*0;
eig2 = A11*0;
eig3 = A11*0;
eig2_t = (A11 + A22 + A33)/3;
p = (A11 - eig2_t)^2 + (A22 - eig2_t)^2 + (A33 - eig2_t)^2 + 2 * p + precision_error;
p = sqrt(p/6);
A11 = (A11 - eig2_t)/p;
A22 = (A22 - eig2_t)/p;
A33 = (A33 - eig2_t)/p;
A12 = A12 / p;
A13 = A13 / p;
A23 = A23 / p;
eig3_t = 0.5 * ( A11 * A22 * A33 + A12 * A23 * A13 * 2 - ...
    (A13*A13) * A22  - (A12*A12) * A33 - (A23*A23) * A11);
eig3_t = acos( min(1, max(-1, eig3_t)) )/3;
eig1_t = eig2_t + 2 * p * cos(eig3_t);
eig3_t = eig2_t + 2 * p * cos(eig3_t + c);
eig2_t = 3 * eig2_t - eig1_t - eig3_t;
% Sort eigenvalues from small to large according to their absolute
abs_e1 = abs(eig1_t);
abs_e2 = abs(eig2_t);
abs_e3 = abs(eig3_t);

if abs_e1 < abs_e2
    if abs_e2 < abs_e3 % 1,2,3
        eig1 = eig1_t;
        eig2 = eig2_t;
        eig3 = eig3_t;
    elseif abs_e1 < abs_e3 % 1,3,2
        eig1 = eig1_t;
        eig2 = eig3_t;
        eig3 = eig2_t;
    else
        eig1 = eig3_t;% 3,1,2
        eig2 = eig1_t;
        eig3 = eig2_t;
    end
else% abs_e2 < abs_e1
    if abs_e1 < abs_e3
        eig1 = eig2_t;
        eig2 = eig1_t;
        eig3 = eig3_t;
    elseif abs_e2 < abs_e3   % abs_e2 < abs_e1 & abs_e3 < abs_e1
        eig1 = eig2_t;
        eig2 = eig3_t;
        eig3 = eig1_t;
    else
        eig1 = eig3_t;
        eig2 = eig2_t;
        eig3 = eig1_t;
    end
end
%% Compute Frangi vesselness
alpha = 0.5;
beta = 0.5;
noise_level = 0.002;
vesselness = eig1*0;
if eig2 >= 0 
    return;
elseif eig3 >= 0
    return;
else
    signal_level = eig1 * eig1 + eig2 * eig2 + eig3 * eig3;
    signal_level = 1 - exp(- signal_level/noise_level);
    vesselness =  (1 - exp(- ((eig2/eig3)^2 ) / alpha)) * ...
        exp(-( eig1^2 /abs(eig2 * eig3) ) / beta) * signal_level;
    if isnan(vesselness)
        vesselness = eig1*0;
    end
end
end

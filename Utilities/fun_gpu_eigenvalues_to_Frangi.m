function vesselness = fun_gpu_eigenvalues_to_Frangi(eig1, eig2, eig3, noise_level)

alpha = 0.5;
beta = 0.5;
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
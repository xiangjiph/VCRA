function pdf_solu = fun_simulation_OT_solve_ct_diff_itr_n(vessel_mask, inhomogeneous_term, ic, useGPUQ)

if nargin < 4
    useGPUQ = true;
end
vis_Q = false;
%% Parameters
% Computation parameters
num_max_iter = 8e3;
target_accuracy = 0.01;
%%
bc_at_mask = 0;

vessel_ind = find(vessel_mask);
int_tissue_mask = ~vessel_mask(2:(end-1), 2:(end-1), 2:(end-1));
sys_size = size(vessel_mask);
%% Initialization 
inhomogeneous_term = single(inhomogeneous_term);
if isscalar(ic)
    simu_array = ones(sys_size, 'single') * ic;
else
    simu_array = ic;
end
simu_array(vessel_ind) = bc_at_mask;
if useGPUQ
    simu_array = gpuArray(simu_array);
end
check_error_period = 100;
max_err = 1;
iter = 0;
simu_tic = tic;
if vis_Q
    fig_hdl = figure;
    ax_hdl = axes(fig_hdl);
    ax_hdl.DataAspectRatio = [1,1,1];
    ax_hdl.CLim = [-600, 0];
end
vis_sec = round(sys_size(3)/2);
while iter < num_max_iter && max_err > target_accuracy
    iter = iter + 1;
    % Update non-boundary voxels    
%     tmp_new_array = fun_add_central_difference_terms_in_dir(simu_array, 1);
%     tmp_new_array = tmp_new_array + fun_add_central_difference_terms_in_dir(simu_array, 2);
%     tmp_new_array = tmp_new_array + fun_add_central_difference_terms_in_dir(simu_array, 3);
    
    tmp_new_array = simu_array(1 : (end - 2), 2 : (end - 1), 2 : (end - 1));
    
    tmp_new_array = tmp_new_array + simu_array(3 : (end    ), 2 : (end - 1), 2 : (end - 1));
    tmp_new_array = tmp_new_array + simu_array(2 : (end - 1), 1 : (end - 2), 2 : (end - 1));
    tmp_new_array = tmp_new_array + simu_array(2 : (end - 1), 3 : (end    ), 2 : (end - 1));
    tmp_new_array = tmp_new_array + simu_array(2 : (end - 1), 2 : (end - 1), 1 : (end - 2));
    tmp_new_array = tmp_new_array + simu_array(2 : (end - 1), 2 : (end - 1), 3 : (end    ));
    
%     tmp_new_array = bsxfun(@plus, tmp_new_array, simu_array(3 : (end    ), 2 : (end - 1), 2 : (end - 1)));
%     tmp_new_array = bsxfun(@plus, tmp_new_array, simu_array(2 : (end - 1), 1 : (end - 2), 2 : (end - 1)));
%     tmp_new_array = bsxfun(@plus, tmp_new_array, simu_array(2 : (end - 1), 3 : (end    ), 2 : (end - 1)));
%     tmp_new_array = bsxfun(@plus, tmp_new_array, simu_array(2 : (end - 1), 2 : (end - 1), 1 : (end - 2)));
%     tmp_new_array = bsxfun(@plus, tmp_new_array, simu_array(2 : (end - 1), 2 : (end - 1), 3 : (end    )));
%     
    tmp_new_array = tmp_new_array ./ 6 - inhomogeneous_term / 6;
    if mod(iter, check_error_period) == 0
        max_err = gather(bsxfun(@minus, simu_array(2 : (end - 1), 2 : (end - 1), 2 : (end - 1)), tmp_new_array));
        max_err = max_err(int_tissue_mask);
        max_err = max(abs(max_err));
        fprintf('Finish computing iteration %d. Maximum update is %f.\n', ...
            iter, max_err);
        if vis_Q
            hold(ax_hdl, 'on');
            imagesc(ax_hdl, simu_array(:, :, vis_sec));
            hold(ax_hdl, 'off');
            colorbar(ax_hdl);
            drawnow
        end
    end
    simu_array(2 : (end - 1), 2 : (end - 1), 2 : (end - 1)) = tmp_new_array;
    % Update boundary voxels
    simu_array(1, :, :) = simu_array(3, :, :);
    simu_array(end, :, :) = simu_array(end - 2, :, :);
    simu_array(:, 1, :) = simu_array(:, 3, :);
    simu_array(:, end, :) = simu_array(:, end - 2, :);
    simu_array(:, :, 1) = simu_array(:, :, 3);
    simu_array(:, :, end) = simu_array(:, :, end - 2);
    % Fix the oxygen partial pressure on the boundary
    simu_array(vessel_ind) = bc_at_mask;
end
fprintf('Finish solving Poisson equation. Elapse time is %f seconds. Final maximum voxel value updated was %f\n', ...
    toc(simu_tic), max_err);
if useGPUQ
    simu_array = gather(simu_array);
end
%% Result
pdf_solu.vessel_mask = vessel_mask;
pdf_solu.bc_at_mask = bc_at_mask;
pdf_solu.pO2_array = simu_array;
pdf_solu.inhomogeneous_term = inhomogeneous_term;
pdf_solu.final_maximum_update = max_err;
pdf_solu.final_iteration = iter;
end
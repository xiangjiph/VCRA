function [eig1, eig2, eig3, v1_1, v1_2, v1_3, v2_1, v2_2, v2_3, v3_1, v3_2, v3_3]...
    = fun_vectorized_3x3_real_sym_mat_eigenvectors(A11, A12, A13, A22, A23, A33, eig1, eig2, eig3, eig_order)
% fun_vectorized_3x3_real_sym_mat_eigenvectors computes the eigenvectors
% and (corrected) eigenvalues, return [eig1, eig2, eig3, v1_1, v1_2, v1_3,
% v2_1, v2_2, v2_3, v3_1, v3_2, v3_3] in vectors
% Input: 
%   Aij: N * 1 vector of the (i,j) element of real symmetric matrix
%   eig1: the largest eigenvale of the input matrix
% Output: 
%   eig1: (corrected) eigenvalues
%   v1_1: the first component of eigenvector 1( correspond to eigenvalue 1)
% Note
%   Debug file: debug_eigenvectors_vectorized_simplified.m
if nargin < 10
    eig_order = 'abs_eig1_smallest';
end
coeff_error_thresh = 50;
coeff_error_small_angle = 256;
expected_max_error = (1e-10)^2;
% expected_max_error = 1;
precision_epsilon = eps(class(eig1));
switch eig_order
    case 'abs_eig1_smallest'
        thresh_1 = eps(eig3);
        thresh = ( coeff_error_thresh * precision_epsilon .* eig3) .^2;
    case 'eig1_largest'
        lambda_abs_max = max(abs([eig1'; eig3']), [], 1)';
        thresh_1 = eps(lambda_abs_max);
        thresh = ( coeff_error_thresh * precision_epsilon .* lambda_abs_max) .^2;
        clear lambda_abs_max;
end
% Test if the matrix is degenerated
degenerated_1_2 = abs(eig1 - eig2) < thresh_1;
three_fold_degenerated_Q = degenerated_1_2 & (abs(eig2 - eig3) < thresh_1);
% #############################
% Compute the frist eigenvector
% #############################
disp('Compute eigenvector 1');
processed_eig_list1 = false(length(eig1),1);
E1_11 = A11 - eig1;
E1_22 = A22 - eig1;
v1_1 = A12 .* A23 - A13 .* E1_22;
v1_2 = A13 .* A12 - E1_11 .* A23;
v1_3 = E1_11 .* E1_22 - A12 .^2;
% 3 fold degenerated
% if any(three_fold_degenerated_Q)
    
if ~isempty(find(three_fold_degenerated_Q, 1))
    disp('Eigenvector 1: 3 fold degenerated eigenvalues');
    v1_1(three_fold_degenerated_Q) = 1;
    v1_2(three_fold_degenerated_Q) = 0;
    v1_3(three_fold_degenerated_Q) = 0;
%     v1_1 = v1_1 + (1 - v1_1) .* three_fold_degenerated_Q;
%     v1_2 = v1_2 - v1_2 .* three_fold_degenerated_Q;
%     v1_3 = v1_3 - v1_3 .* three_fold_degenerated_Q;
    processed_eig_list1 = processed_eig_list1 | three_fold_degenerated_Q;
end
% column 1 is [0,0,0]
n1 = E1_11.^2 + A12.^2 + A13.^2;
colume_one_is_zero = (n1 <= thresh) & (~processed_eig_list1);
if any(colume_one_is_zero)
    disp('Eigenvector 1: First column is approximately 0');
    processed_eig_list1 = processed_eig_list1 | colume_one_is_zero;
    v1_1(colume_one_is_zero) = 1;
    v1_2(colume_one_is_zero) = 0;
    v1_3(colume_one_is_zero) = 0;
%     v1_1 = v1_1 + (1 - v1_1) .* colume_one_is_zero;
%     v1_2 = v1_2 - v1_2 .* colume_one_is_zero;
%     v1_3 = v1_3 - v1_3 .* colume_one_is_zero;
end
% colume 2 is [0,0,0] - The following implements have some overlap with the
% last part, i.e. colume 1 and colume 2 are both 0
n2 = A12.^2 + E1_22.^2 + A23.^2; 
errorVal = n1 .* n2;
clear n1 colume_one_is_zero;
colume_two_is_zero = (n2 <= thresh) & (~processed_eig_list1);
if any(colume_two_is_zero)
    disp('Eigenvector 1: Second column is approximately 0');
    processed_eig_list1 = processed_eig_list1 | colume_two_is_zero;
%     v1_1 = v1_1 - v1_1 .* colume_two_is_zero;
%     v1_2 = v1_2 + (1 - v1_2) .* colume_two_is_zero;
%     v1_3 = v1_3 - v1_3 .* colume_two_is_zero;
    v1_1(colume_two_is_zero) = 0;
    v1_2(colume_two_is_zero) = 1;
    v1_3(colume_two_is_zero) = 0;
end
clear n2 colume_two_is_zero;
% Colume 1 and 2 are not 0, but parallel 
v1_norm = v1_1 .^2 + v1_2 .^2 + v1_3.^2;
v1_small_angle = (v1_norm < (coeff_error_small_angle * precision_epsilon)^2 .* errorVal) & (~processed_eig_list1);
if any(v1_small_angle)
    disp('Eigenvector 1: First two column vectors are in parallel');
    E1_2 = [A12(v1_small_angle)';E1_22(v1_small_angle)'; A23(v1_small_angle)'];
    [max_V, max_idx] = max(abs([E1_11(v1_small_angle)'; A12(v1_small_angle)'; A13(v1_small_angle)']),[], 1);
    tmp_f = - max_V ./ E1_2(sub2ind(size(E1_2), max_idx, 1:length(max_idx)));
    v1_1(v1_small_angle) = 1;
    v1_2(v1_small_angle) = tmp_f;
    v1_3(v1_small_angle) = 0;
%     v1_norm(v1_small_angle) = 1 + tmp_f.^2;
end
clear E1_2 E1_22 E1_11 max_V max_idx tmp_f processed_eig_list1;
% Normalized:
if ~isempty(find(v1_small_angle, 1))
    v1_norm = sqrt(v1_1 .^2 + v1_2 .^2 + v1_3.^2);
end
v1_1 = v1_1 ./ v1_norm;
v1_2 = v1_2 ./ v1_norm;
v1_3 = v1_3 ./ v1_norm;
clear v1_norm v1_small_angle
% ##############################
% Compute the second eigenvector
% ##############################
disp('Compute eigenvector 2');
processed_eig_list2 = false(length(eig1),1);
E2_11 = A11 - eig2;
E2_22 = A22 - eig2;
E2_33 = A33 - eig2;
v2_1 = A12 .* A23 - A13 .* E2_22;
v2_2 = A13 .* A12 - E2_11 .* A23;
v2_3 = E2_11 .* E2_22 - A12 .^2;
% 3 fold degenerated
if any(three_fold_degenerated_Q)
    disp('Eigenvector 2: 3 fold eigenvalues');
    v2_1 = v2_1 - v2_1 .* three_fold_degenerated_Q;
    v2_2 = v2_2 + (1 - v2_2) .* three_fold_degenerated_Q;
    v2_3 = v2_3 - v2_3 .* three_fold_degenerated_Q;
%     v2_1(three_fold_degenerated_Q) = 0;
%     v2_2(three_fold_degenerated_Q) = 1;
%     v2_3(three_fold_degenerated_Q) = 0;
    processed_eig_list2 = processed_eig_list2 | three_fold_degenerated_Q;
end
% colume 1 is [0,0,0]
n1 = E2_11.^2 + A12.^2 + A13.^2;
colume_one_is_zero = (n1 <= thresh) & (~processed_eig_list2);
if any(colume_one_is_zero)
    disp('Eigenvector 2: The first column vector is 0');
    v2_1 = v2_1 + (1 - v2_1) .* colume_one_is_zero;
    v2_2 = v2_2 - v2_2 .* colume_one_is_zero;
    v2_3 = v2_3 - v2_3 .* colume_one_is_zero;
%     v2_1(colume_one_is_zero) = 1;
%     v2_2(colume_one_is_zero) = 0;
%     v2_3(colume_one_is_zero) = 0;
    processed_eig_list2 = processed_eig_list2 | colume_one_is_zero;
end
clear colume_one_is_zero
% colume 2 is [0,0,0] - The following implements have some overlap with the
% last part, i.e. colume 1 and colume 2 are both 0
n2 = A12.^2 + E2_22.^2 + A23.^2; 
errorVal = n1 .* n2;

colume_two_is_zero = (n2 <= thresh) & (~processed_eig_list2);
if any(colume_two_is_zero)
    disp('Eigenvector 2: The second column vector is 0');
    processed_eig_list2 = processed_eig_list2 | colume_two_is_zero;
    v2_1 = v2_1 - v2_1 .* colume_two_is_zero;
    v2_2 = v2_2 + (1 - v2_2) .* colume_two_is_zero;
    v2_3 = v2_3 - v2_3 .* colume_two_is_zero;
%     v2_1(colume_two_is_zero) = 0;
%     v2_2(colume_two_is_zero) = 1;
%     v2_3(colume_two_is_zero) = 0;
end
clear colume_two_is_zero
% Colume 1 and 2 are not 0, but parallel 
v2_norm = v2_1 .^2 + v2_2 .^2 + v2_3.^2;

v2_small_angle = (v2_norm < (coeff_error_small_angle * precision_epsilon)^2 .* errorVal) & (~processed_eig_list2);
if any(v2_small_angle)
    disp('Eigenvector 2: The first two column vectors are in parallel');
    E2_2 = [A12(v2_small_angle)';E2_22(v2_small_angle)'; A23(v2_small_angle)'];
    [max_V, max_idx] = max(abs([E2_11(v2_small_angle)'; A12(v2_small_angle)'; A13(v2_small_angle)']),[], 1);
    tmp_f = - max_V ./ E2_2(sub2ind(size(E2_2), max_idx, 1:length(max_idx)));
    v2_1(v2_small_angle) = 1;
    v2_2(v2_small_angle) = tmp_f;
    v2_3(v2_small_angle) = 0;
end
clear v2_small_angle E2_2 max_V max_idx tmp_f processed_eig_list2 errorVal
% two fold Degenerated eigenvalues
degenerated_1_2_left = degenerated_1_2 & (~three_fold_degenerated_Q);
clear degenerated_1_2 three_fold_degenerated_Q
if any(degenerated_1_2_left)
    disp('Eigenvector 2: The first two eigenvalues are two fold degenerated');
    % If colume 1 of (A - eig1 * I) is not zero, take the cross product with v1
    % Find the index of eigen matrix with non-zero colume one. 
    col_v1_valid = find(degenerated_1_2_left & (n1 > thresh)); 
    if ~isempty(col_v1_valid)
        disp('Use the first column vector to take the cross product with the first eigenvector');
        % Calculate v2 for degenerated eigenvalue
        deg_v2_1_list = v1_2(col_v1_valid) .* A13(col_v1_valid) - v1_3(col_v1_valid) .* A12(col_v1_valid);
        deg_v2_2_list = v1_3(col_v1_valid) .* E2_11(col_v1_valid) - v1_1(col_v1_valid) .* A13(col_v1_valid);
        deg_v2_3_list = v1_1(col_v1_valid) .* A12(col_v1_valid) - v1_2(col_v1_valid) .* E2_11(col_v1_valid);
        % The resulting eigenvector is valid if it's norm is not zero.
        deg_v2_norm = deg_v2_1_list .^2 + deg_v2_2_list.^2 + deg_v2_3_list.^2;
        deg_v2_valid_Q = deg_v2_norm > ((256 * precision_epsilon)^2 .* n1(col_v1_valid));
        % Remove the index of eigen matrix that has 0 norm eigenvector calculated
        % above
        col_v1_valid(~deg_v2_valid_Q) = [];
        % Remove the processed eigenvalue from the list
        degenerated_1_2_left(col_v1_valid) = false;
        % Update the degenerated eigenvector
        v2_1(col_v1_valid) = deg_v2_1_list(deg_v2_valid_Q);
        v2_2(col_v1_valid) = deg_v2_2_list(deg_v2_valid_Q);
        v2_3(col_v1_valid) = deg_v2_3_list(deg_v2_valid_Q);
        % Update the degenerated eigenvector norm
    end
    clear col_v1_valid E2_11
    
    if any(degenerated_1_2_left)
        % If colume 2 of (A - eig1 * I) is not zero, take the cross product with v1
        col_v2_valid = find(degenerated_1_2_left & (n2 > thresh));
        if ~isempty(col_v2_valid)
            disp('Use the second column vector to take the cross product with the first eigenvector');
            deg_v2_1_list = v1_2(col_v2_valid) .* A23(col_v2_valid) - v1_3(col_v2_valid) .* E2_22(col_v2_valid);
            deg_v2_2_list = v1_3(col_v2_valid) .* A12(col_v2_valid) - v1_1(col_v2_valid) .* A23(col_v2_valid);
            deg_v2_3_list = v1_1(col_v2_valid) .* E2_22(col_v2_valid) - v1_2(col_v2_valid) .* A12(col_v2_valid);
            deg_v2_norm = deg_v2_1_list .^2 + deg_v2_2_list.^2 + deg_v2_3_list.^2;
            deg_v2_valid_Q = deg_v2_norm > ((256 * precision_epsilon)^2 .* n2(col_v2_valid));
            col_v2_valid(~deg_v2_valid_Q) = [];
            degenerated_1_2_left(col_v2_valid) = false;
            v2_1(col_v2_valid) = deg_v2_1_list(deg_v2_valid_Q);
            v2_2(col_v2_valid) = deg_v2_2_list(deg_v2_valid_Q);
            v2_3(col_v2_valid) = deg_v2_3_list(deg_v2_valid_Q);
        end
        clear col_v2_valid E2_22
        if any(degenerated_1_2_left)
            % If colume 3 of (A - eig1 * I) is not zero, take the cross product with v1
            n3 = A13.^2 + A23.^2 + (A33 - eig2).^2;
            col_v3_valid = find(degenerated_1_2_left & (n3 > thresh));
            if ~isempty(col_v3_valid)
                disp('Use the third column vector to take the cross product with the first eigenvector');
                deg_v2_1_list = v1_2(col_v3_valid) .* E2_33(col_v3_valid) - v1_3(col_v3_valid) .* A23(col_v3_valid);
                deg_v2_2_list = v1_3(col_v3_valid) .* A13(col_v3_valid) - v1_1(col_v3_valid) .* E2_33(col_v3_valid);
                deg_v2_3_list = v1_1(col_v3_valid) .* A23(col_v3_valid) - v1_2(col_v3_valid) .* A13(col_v3_valid);
                deg_v2_norm = deg_v2_1_list .^2 + deg_v2_2_list.^2 + deg_v2_3_list.^2;
                deg_v2_valid_Q = deg_v2_norm > ((256 * precision_epsilon)^2 .* n3(col_v3_valid));
                col_v3_valid(~deg_v2_valid_Q) = [];
                degenerated_1_2_left(col_v3_valid) = false;
                v2_1(col_v3_valid) = deg_v2_1_list(deg_v2_valid_Q);
                v2_2(col_v3_valid) = deg_v2_2_list(deg_v2_valid_Q);
                v2_3(col_v3_valid) = deg_v2_3_list(deg_v2_valid_Q);
            end
            clear col_v3_valid E2_33
            % Any vector orthogonal to v1 is eigenvector
            if any(degenerated_1_2_left)
                disp(' Any vector orthogonal to the first eigenvector is an eigenvector');
                degenerated_1_2_left_idx = find(degenerated_1_2_left);
                v2_1(degenerated_1_2_left_idx) = 0;
                v2_2(degenerated_1_2_left_idx) = 0;
                v2_3(degenerated_1_2_left_idx) = 0;
                deg_v1_1_list = v1_1(degenerated_1_2_left_idx);
                deg_v1_2_list = v1_2(degenerated_1_2_left_idx);
                deg_v1_3_list = v1_3(degenerated_1_2_left_idx);
                % If v1_1 is the first non-zero term
                valid_deg_v1_1 = deg_v1_1_list > thresh_1;
                if any(valid_deg_v1_1)
                    disp('The first component of eigenvector 1 is nonzero');
                    v2_1(degenerated_1_2_left_idx(valid_deg_v1_1)) = deg_v1_2_list(valid_deg_v1_1);
                    v2_2(degenerated_1_2_left_idx(valid_deg_v1_1)) = - deg_v1_1_list(valid_deg_v1_1);
                    degenerated_1_2_left_idx(valid_deg_v1_1) = [];
                end
                
                if ~isempty(degenerated_1_2_left_idx)
                    % If v1_2 is the first non-zero term
                    deg_v1_1_list = deg_v1_1_list(~valid_deg_v1_1);
                    deg_v1_2_list = deg_v1_2_list(~valid_deg_v1_1);
                    deg_v1_3_list = deg_v1_3_list(~valid_deg_v1_1);
                    valid_deg_v1_2 = deg_v1_2_list > thresh_1;
                    if any(valid_deg_v1_2)
                        disp('The second component of eigenvector 1 is nonzero');
                        v2_2(degenerated_1_2_left_idx(valid_deg_v1_2)) = deg_v1_3_list(valid_deg_v1_2);
                        v2_3(degenerated_1_2_left_idx(valid_deg_v1_2)) = - deg_v1_2_list(valid_deg_v1_2);
                        degenerated_1_2_left_idx(valid_deg_v1_2) = [];
                    end
                    
                    if ~isempty(degenerated_1_2_left_idx)
                        % If v1_3 is the first non-zero term
                        deg_v1_1_list = deg_v1_1_list(~valid_deg_v1_2);
%                         deg_v1_2_list = deg_v1_2_list(~valid_deg_v1_2);
                        deg_v1_3_list = deg_v1_3_list(~valid_deg_v1_2);
                        valid_deg_v1_3 = deg_v1_3_list > thresh_1;
                        if any(valid_deg_v1_3)
                            disp('The third component of eigenvector 1 is nonzero');
                            v2_3(degenerated_1_2_left_idx(valid_deg_v1_3)) = deg_v1_1_list(valid_deg_v1_3);
                            v2_1(degenerated_1_2_left_idx(valid_deg_v1_3)) = - deg_v1_3_list(valid_deg_v1_3);
                            degenerated_1_2_left_idx(valid_deg_v1_3) = [];
                        end
                        if any(degenerated_1_2_left_idx)
                            warning('Condition not considered! Index:');
                            disp(find(degenerated_1_2_left_idx));                            
                        end
                    end
                end
            end
        end
    end
end
clear degenerated_1_2_left_idx degenerated_1_2_left n2 n1 n3 thresh thresh_1
% Normalized the second eigenvector
v2_norm = sqrt(v2_1 .^2 + v2_2 .^2 + v2_3.^2);
v2_1 = v2_1 ./ v2_norm;
v2_2 = v2_2 ./ v2_norm;
v2_3 = v2_3 ./ v2_norm;
clear v2_norm
% ###################################################################
% Compute the third eigenvector as the cross product of the first two
% ###################################################################
disp('Compute eigenvector 3');
v3_1 = v1_2 .* v2_3 - v1_3 .* v2_2;
v3_2 = v1_3 .* v2_1 - v1_1 .* v2_3;
v3_3 = v1_1 .* v2_2 - v1_2 .* v2_1;
% #########################################
% Self-consistant evaluation and correction
% #########################################
if autoCorrectionQ
    error_check1 = (A11 - eig1) .* v1_1 + A12 .* v1_2 + A13 .* v1_3;
    error_check2 = A12 .* v1_1 + (A22 - eig1) .* v1_2 + A23 .* v1_3;
    error_check3 = A13 .* v1_1 + A23 .* v1_2 + (A33 - eig1) .* v1_3;
    v1_norm_error = (error_check1.^2 + error_check2.^2 + error_check3.^2)';
    error_check1 = (A11 - eig2) .* v2_1 + A12 .* v2_2 + A13 .* v2_3;
    error_check2 = A12 .* v2_1 + (A22 - eig2) .* v2_2 + A23 .* v2_3;
    error_check3 = A13 .* v2_1 + A23 .* v2_2 + (A33 - eig2) .* v2_3;
    v2_norm_error = (error_check1.^2 + error_check2.^2 + error_check3.^2)';
    error_check1 = (A11 - eig3) .* v3_1 + A12 .* v3_2 + A13 .* v3_3;
    error_check2 = A12 .* v3_1 + (A22 - eig3) .* v3_2 + A23 .* v3_3;
    error_check3 = A13 .* v3_1 + A23 .* v3_2 + (A33 - eig3) .* v3_3;
    v3_norm_error = (error_check1.^2 + error_check2.^2 + error_check3.^2)';
    clear error_check1 error_check2 error_check3
    % Apply Correction
    % Fix the probably wrong eigensystems with MATLAB's built-in eigs
    all_need_correction_idx = find((v1_norm_error > expected_max_error) | (v2_norm_error > expected_max_error) | (v3_norm_error > expected_max_error));
    clear v1_norm_error v2_norm_error v3_norm_error
    if ~isempty(all_need_correction_idx)
        fprintf('Correct %d numerical error using built-in eigs function\n', length(all_need_correction_idx));
        num_need_correction = length(all_need_correction_idx);
        fix_matList = zeros(3,3,num_need_correction);
        fix_matList(1,1,:) = A11(all_need_correction_idx);
        fix_matList(1,2,:) = A12(all_need_correction_idx);
        fix_matList(1,3,:) = A13(all_need_correction_idx);
        fix_matList(2,1,:) = A12(all_need_correction_idx);
        fix_matList(2,2,:) = A22(all_need_correction_idx);
        fix_matList(2,3,:) = A23(all_need_correction_idx);
        fix_matList(3,1,:) = A13(all_need_correction_idx);
        fix_matList(3,2,:) = A23(all_need_correction_idx);
        fix_matList(3,3,:) = A33(all_need_correction_idx);
        fix_eig_mat = zeros(3,num_need_correction);
        vec_list = zeros(3,3,num_need_correction);
        for mat_idx = 1 : num_need_correction
            [vec_list(:,:,mat_idx), fix_eig_mat(:,mat_idx)] = eig(fix_matList(:,:,mat_idx), 'vector');
        end
        clear fix_matList
        [fix_eig_list_sort, fix_eig_idx_sort] = sort(fix_eig_mat, 1,'ascend','ComparisonMethod', 'abs');
        clear fix_eig_list fix_eig_mat
        fix_eig_list_sort = fix_eig_list_sort';
        fix_eig_idx_sort = fix_eig_idx_sort';
        eig1(all_need_correction_idx) = fix_eig_list_sort(:,1);
        eig2(all_need_correction_idx) = fix_eig_list_sort(:,2);
        eig3(all_need_correction_idx) = fix_eig_list_sort(:,3);
        clear fix_eig_list_sort
        v1_1(all_need_correction_idx) = vec_list(sub2ind(size(vec_list),     ones(1,num_need_correction), fix_eig_idx_sort(:,1)', 1:num_need_correction));
        v1_2(all_need_correction_idx) = vec_list(sub2ind(size(vec_list), 2 * ones(1,num_need_correction), fix_eig_idx_sort(:,1)', 1:num_need_correction));
        v1_3(all_need_correction_idx) = vec_list(sub2ind(size(vec_list), 3 * ones(1,num_need_correction), fix_eig_idx_sort(:,1)', 1:num_need_correction));
        v2_1(all_need_correction_idx) = vec_list(sub2ind(size(vec_list),     ones(1,num_need_correction), fix_eig_idx_sort(:,2)', 1:num_need_correction));
        v2_2(all_need_correction_idx) = vec_list(sub2ind(size(vec_list), 2 * ones(1,num_need_correction), fix_eig_idx_sort(:,2)', 1:num_need_correction));
        v2_3(all_need_correction_idx) = vec_list(sub2ind(size(vec_list), 3 * ones(1,num_need_correction), fix_eig_idx_sort(:,2)', 1:num_need_correction));
        v3_1(all_need_correction_idx) = vec_list(sub2ind(size(vec_list),     ones(1,num_need_correction), fix_eig_idx_sort(:,3)', 1:num_need_correction));
        v3_2(all_need_correction_idx) = vec_list(sub2ind(size(vec_list), 2 * ones(1,num_need_correction), fix_eig_idx_sort(:,3)', 1:num_need_correction));
        v3_3(all_need_correction_idx) = vec_list(sub2ind(size(vec_list), 3 * ones(1,num_need_correction), fix_eig_idx_sort(:,3)', 1:num_need_correction));
    end
end

end

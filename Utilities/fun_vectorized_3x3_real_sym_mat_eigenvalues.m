function [eig1, eig2, eig3] = fun_vectorized_3x3_real_sym_mat_eigenvalues(A11,A12,A13,A22,A23,A33, eig_order)
% fun_vectorized_3x3_sym_real_mat_diagonalization(A11, A12, A13, A22, A23, A33)
% Assuming there are not much diagonalized hessian matrix
% Input: 
%   A12: vector of (1,2) element of the matrics
%   eig_order: string, specify the order of output eigenvalues
%       eig1_largest: eig1 >= eig2 >= eig3
%       abs_eig1_smallest: |eig1| <= |eig2| <= |eig3|, default option

% switch class(A11)
%     case {'float', 'double', 'uint8', 'int8', 'uint16', 'int16', 'uint32', 'int32'}
%         A11 = single(A11);
%         A12 = single(A12);
%         A13 = single(A13);
%         A22 = single(A22);
%         A23 = single(A23);
%         A33 = single(A33);
% end
if nargin < 7
    eig_order = 'abs_eig1_smallest';
end

c = 2 * pi/3;
precision_error = eps('double');
p = A12.^2 + A13.^2 + A23.^2;

switch eig_order
    case 'eig1_largest'
        eig2 = (A11 + A22 + A33)./3;
        % Add a very small constant to avoid divided by 0 error.
        p = (A11 - eig2).^2 + (A22 - eig2).^2 + (A33 - eig2).^2 + 2 * p + precision_error;
        p = sqrt(p./6);
        A11 = (A11 - eig2)./p;
        A22 = (A22 - eig2)./p;
        A33 = (A33 - eig2)./p;
        A12 = A12 ./ p;
        A13 = A13 ./ p;
        A23 = A23 ./ p;
        eig3 = 0.5 .* ( A11 .* A22 .* A33 + A12 .* A23 .* A13 .* 2 - ...
            (A13.^2) .* A22  - (A12.^2) .* A33 - (A23.^2) .* A11);
        if gpuUsedQ
            clear A11 A12 A13 A22 A23 A33
        end
        eig3 = acos( min(1, max(-1, eig3)) )/3;
        % eig1 >= eig2 >= eig3
        eig1 = eig2 + 2 .* p .* cos(eig3);
        eig3 = eig2 + 2 .* p .* cos(eig3 + c);
        eig2 = 3 .* eig2 - eig1 - eig3;
    case 'abs_eig1_smallest'
        eig2_t = (A11 + A22 + A33)./3;
        p = (A11 - eig2_t).^2 + (A22 - eig2_t).^2 + (A33 - eig2_t).^2 + 2 * p + precision_error;
        p = sqrt(p./6);
        A11 = (A11 - eig2_t)./p;
        A22 = (A22 - eig2_t)./p;
        A33 = (A33 - eig2_t)./p;
        A12 = A12 ./ p;
        A13 = A13 ./ p;
        A23 = A23 ./ p;
        eig3_t = 0.5 .* ( A11 .* A22 .* A33 + A12 .* A23 .* A13 .* 2 - ...
            (A13.^2) .* A22  - (A12.^2) .* A33 - (A23.^2) .* A11);
        clear A11 A12 A13 A22 A23 A33
        eig3_t = acos( min(1, max(-1, eig3_t)) )/3;
        eig1_t = eig2_t + 2 .* p .* cos(eig3_t);
        eig3_t = eig2_t + 2 .* p .* cos(eig3_t + c);
        eig2_t = 3 .* eig2_t - eig1_t - eig3_t;
        % Sort eigenvalues from small to large according to their absolute
        cp_1_2 = abs(eig1_t) >= abs(eig2_t);
        cp_1_3 = abs(eig1_t) >= abs(eig3_t);
        cp_2_3 = abs(eig2_t) >= abs(eig3_t);
        od123 = (~cp_1_2) & (~cp_1_3) & (~cp_2_3);
        od132 = (~cp_1_2) & (~cp_1_3) & (cp_2_3);
        od213 = (cp_1_2) & (~cp_1_3) & (~cp_2_3);
        od231 = (cp_1_2) & (cp_1_3) & (~cp_2_3);
        od312 = (cp_1_3) & (~cp_1_2) & (cp_2_3);
        od321 = (cp_1_3) & (cp_1_2) & (cp_2_3);
        eig1 = eig1_t .* od123 + eig1_t .* od132 + eig2_t .* od213 + eig2_t .* od231 + eig3_t .* od312 + eig3_t .* od321;
        eig2 = eig2_t .* od123 + eig3_t .* od132 + eig1_t .* od213 + eig3_t .* od231 + eig1_t .* od312 + eig2_t .* od321;
        eig3 = eig3_t .* od123 + eig2_t .* od132 + eig3_t .* od213 + eig1_t .* od231 + eig2_t .* od312 + eig1_t .* od321;
end
end

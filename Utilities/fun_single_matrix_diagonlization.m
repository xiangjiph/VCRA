function eigs = fun_single_matrix_diagonlization(A)

c1 = 2 * pi/3;
p1 = A(1,2)^2 + A(1,3)^2 + A(2,3)^2;
q = trace(A)/3;
p2 = ( A(1,1) - q )^2 + ( A(2,2) - q )^2 + ( A(3,3) - q)^2 + 2 * p1;
p = sqrt(p2/6);
B = (1/p) * (A - q * eye([3,3]));

% r = 0.5* (B(1,1) * B(2,2) * B(3,3) + B(1,2) * B(2,3) * B(3,1) + B(2,1) * B(3,2) * B(1,3) - ...
%     B(1,3) * B(2,2) * B(3,1) - B(1,2) * B(2,1) * B(3,3) - B(3,2) * B(2,3) * B(1,1)); 
r = 0.5* (B(1,1) * B(2,2) * B(3,3) + B(1,2) * B(2,3) * B(1,3) *2 - ...
    B(1,3) * B(2,2) * B(1,3) - B(1,2).^2 * B(3,3) - B(2,3) .^ 2* B(1,1));

phi = acos( min(1, max(-1, r)) )/3;
eig1 = q + 2 * p * cos(phi);
eig3 = q + 2 * p * cos(phi + c1);
eig2 = 3 * q - eig1 - eig3;
eigs = [eig1, eig2, eig3];

end

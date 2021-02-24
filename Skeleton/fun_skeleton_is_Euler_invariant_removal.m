function is_euler_inv =  fun_skeleton_is_Euler_invariant_removal(neighbor_mat)
% fun_vectorize_is_Euler_invariant check if the neighbor after removing the
% center point has the same euler characteristic as the one before removal.
% 
% invariant, which is issential for skeletonization. 
% Input:
%   neighbor_mat: N-by-27 bindary array. Each row of the array is the 26
%   neighbors (including the center point at (:, 14))
% Output:
%   is_euler_inv: N-by-1 logical array. True if the removal operation is
%   euler invariant. 

% Adapt from Philip Kollmannsberger's implementation on MATLAB File
% Exchange.

%% Euler characteristic look up table from Lee's paper. 
LUT = zeros(255,1, 'int8');
LUT(1)  =  1;LUT(3)  = -1;LUT(5)  = -1;LUT(7)  =  1;LUT(9)  = -3;LUT(11) = -1;
LUT(13) = -1;LUT(15) =  1;LUT(17) = -1;LUT(19) =  1;LUT(21) =  1;LUT(23) = -1;
LUT(25) =  3;LUT(27) =  1;LUT(29) =  1;LUT(31) = -1;LUT(33) = -3;LUT(35) = -1;
LUT(37) =  3;LUT(39) =  1;LUT(41) =  1;LUT(43) = -1;LUT(45) =  3;LUT(47) =  1;
LUT(49) = -1;LUT(51) =  1;

LUT(53) =  1;LUT(55) = -1;LUT(57) =  3;LUT(59) =  1;LUT(61) =  1;LUT(63) = -1;
LUT(65) = -3;LUT(67) =  3;LUT(69) = -1;LUT(71) =  1;LUT(73) =  1;LUT(75) =  3;
LUT(77) = -1;LUT(79) =  1;LUT(81) = -1;LUT(83) =  1;LUT(85) =  1;LUT(87) = -1;
LUT(89) =  3;LUT(91) =  1;LUT(93) =  1;LUT(95) = -1;LUT(97) =  1;LUT(99) =  3;
LUT(101) =  3;LUT(103) =  1;

LUT(105) =  5;LUT(107) =  3;LUT(109) =  3;LUT(111) =  1;LUT(113) = -1;
LUT(115) =  1;LUT(117) =  1;LUT(119) = -1;LUT(121) =  3;LUT(123) =  1;
LUT(125) =  1;LUT(127) = -1;LUT(129) = -7;LUT(131) = -1;LUT(133) = -1;
LUT(135) =  1;LUT(137) = -3;LUT(139) = -1;LUT(141) = -1;LUT(143) =  1;
LUT(145) = -1;LUT(147) =  1;LUT(149) =  1;LUT(151) = -1;LUT(153) =  3;
LUT(155) =  1;

LUT(157) =  1;LUT(159) = -1;LUT(161) = -3;LUT(163) = -1;LUT(165) =  3;
LUT(167) =  1;LUT(169) =  1;LUT(171) = -1;LUT(173) =  3;LUT(175) =  1;
LUT(177) = -1;LUT(179) =  1;LUT(181) =  1;LUT(183) = -1;LUT(185) =  3;
LUT(187) =  1;LUT(189) =  1;LUT(191) = -1;LUT(193) = -3;LUT(195) =  3;
LUT(197) = -1;LUT(199) =  1;LUT(201) =  1;LUT(203) =  3;LUT(205) = -1;
LUT(207) =  1;LUT(209) = -1;LUT(211) =  1;LUT(213) =  1;LUT(215) = -1;
LUT(217) =  3;LUT(219) =  1;LUT(221) =  1;LUT(223) = -1;LUT(225) =  1;
LUT(227) =  3;LUT(229) =  3;LUT(231) =  1;LUT(233) =  5;LUT(235) =  3;
LUT(237) =  3;LUT(239) =  1;LUT(241) = -1;LUT(243) =  1;LUT(245) =  1;
LUT(247) = -1;LUT(249) =  3;LUT(251) =  1;LUT(253) =  1;LUT(255) = -1;
%% Processing 
% The Euler characteristic of the 26 neighbors can be computed by summing
% up the Euler characteristic of the eight overlapping subblock ( the 8
% octants). Refer to the original paper for the proof and more information.

% Calculate Euler characteristic for each octant and sum up
is_euler_inv = zeros(size(neighbor_mat,1),1, 'int8');
% Octant SWU
bitorTable = uint8([128; 64; 32; 16; 8; 4; 2]);
n = ones(size(neighbor_mat,1),1, 'uint8');
n(neighbor_mat(:,25)==1) = bitor(n(neighbor_mat(:,25)==1), bitorTable(1));
n(neighbor_mat(:,26)==1) = bitor(n(neighbor_mat(:,26)==1), bitorTable(2));
n(neighbor_mat(:,16)==1) = bitor(n(neighbor_mat(:,16)==1), bitorTable(3));
n(neighbor_mat(:,17)==1) = bitor(n(neighbor_mat(:,17)==1), bitorTable(4));
n(neighbor_mat(:,22)==1) = bitor(n(neighbor_mat(:,22)==1), bitorTable(5));
n(neighbor_mat(:,23)==1) = bitor(n(neighbor_mat(:,23)==1), bitorTable(6));
n(neighbor_mat(:,13)==1) = bitor(n(neighbor_mat(:,13)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant SEU
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:,27)==1) = bitor(n(neighbor_mat(:,27)==1), bitorTable(1));
n(neighbor_mat(:,24)==1) = bitor(n(neighbor_mat(:,24)==1), bitorTable(2));
n(neighbor_mat(:,18)==1) = bitor(n(neighbor_mat(:,18)==1), bitorTable(3));
n(neighbor_mat(:,15)==1) = bitor(n(neighbor_mat(:,15)==1), bitorTable(4));
n(neighbor_mat(:,26)==1) = bitor(n(neighbor_mat(:,26)==1), bitorTable(5));
n(neighbor_mat(:,23)==1) = bitor(n(neighbor_mat(:,23)==1), bitorTable(6));
n(neighbor_mat(:,17)==1) = bitor(n(neighbor_mat(:,17)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant NWU
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:,19)==1) = bitor(n(neighbor_mat(:,19)==1), bitorTable(1));
n(neighbor_mat(:,22)==1) = bitor(n(neighbor_mat(:,22)==1), bitorTable(2));
n(neighbor_mat(:,10)==1) = bitor(n(neighbor_mat(:,10)==1), bitorTable(3));
n(neighbor_mat(:,13)==1) = bitor(n(neighbor_mat(:,13)==1), bitorTable(4));
n(neighbor_mat(:,20)==1) = bitor(n(neighbor_mat(:,20)==1), bitorTable(5));
n(neighbor_mat(:,23)==1) = bitor(n(neighbor_mat(:,23)==1), bitorTable(6));
n(neighbor_mat(:,11)==1) = bitor(n(neighbor_mat(:,11)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant NEU
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:,21)==1) = bitor(n(neighbor_mat(:,21)==1), bitorTable(1));
n(neighbor_mat(:,24)==1) = bitor(n(neighbor_mat(:,24)==1), bitorTable(2));
n(neighbor_mat(:,20)==1) = bitor(n(neighbor_mat(:,20)==1), bitorTable(3));
n(neighbor_mat(:,23)==1) = bitor(n(neighbor_mat(:,23)==1), bitorTable(4));
n(neighbor_mat(:,12)==1) = bitor(n(neighbor_mat(:,12)==1), bitorTable(5));
n(neighbor_mat(:,15)==1) = bitor(n(neighbor_mat(:,15)==1), bitorTable(6));
n(neighbor_mat(:,11)==1) = bitor(n(neighbor_mat(:,11)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant SWB
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:, 7)==1) = bitor(n(neighbor_mat(:, 7)==1), bitorTable(1));
n(neighbor_mat(:,16)==1) = bitor(n(neighbor_mat(:,16)==1), bitorTable(2));
n(neighbor_mat(:, 8)==1) = bitor(n(neighbor_mat(:, 8)==1), bitorTable(3));
n(neighbor_mat(:,17)==1) = bitor(n(neighbor_mat(:,17)==1), bitorTable(4));
n(neighbor_mat(:, 4)==1) = bitor(n(neighbor_mat(:, 4)==1), bitorTable(5));
n(neighbor_mat(:,13)==1) = bitor(n(neighbor_mat(:,13)==1), bitorTable(6));
n(neighbor_mat(:, 5)==1) = bitor(n(neighbor_mat(:, 5)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant SEB
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:, 9)==1) = bitor(n(neighbor_mat(:, 9)==1), bitorTable(1));
n(neighbor_mat(:, 8)==1) = bitor(n(neighbor_mat(:, 8)==1), bitorTable(2));
n(neighbor_mat(:,18)==1) = bitor(n(neighbor_mat(:,18)==1), bitorTable(3));
n(neighbor_mat(:,17)==1) = bitor(n(neighbor_mat(:,17)==1), bitorTable(4));
n(neighbor_mat(:, 6)==1) = bitor(n(neighbor_mat(:, 6)==1), bitorTable(5));
n(neighbor_mat(:, 5)==1) = bitor(n(neighbor_mat(:, 5)==1), bitorTable(6));
n(neighbor_mat(:,15)==1) = bitor(n(neighbor_mat(:,15)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant NWB
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:, 1)==1) = bitor(n(neighbor_mat(:, 1)==1), bitorTable(1));
n(neighbor_mat(:,10)==1) = bitor(n(neighbor_mat(:,10)==1), bitorTable(2));
n(neighbor_mat(:, 4)==1) = bitor(n(neighbor_mat(:, 4)==1), bitorTable(3));
n(neighbor_mat(:,13)==1) = bitor(n(neighbor_mat(:,13)==1), bitorTable(4));
n(neighbor_mat(:, 2)==1) = bitor(n(neighbor_mat(:, 2)==1), bitorTable(5));
n(neighbor_mat(:,11)==1) = bitor(n(neighbor_mat(:,11)==1), bitorTable(6));
n(neighbor_mat(:, 5)==1) = bitor(n(neighbor_mat(:, 5)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);
% Octant NEB
n = ones(size(neighbor_mat,1),1, 'uint8'); 
n(neighbor_mat(:, 3)==1) = bitor(n(neighbor_mat(:, 3)==1), bitorTable(1));
n(neighbor_mat(:, 2)==1) = bitor(n(neighbor_mat(:, 2)==1), bitorTable(2));
n(neighbor_mat(:,12)==1) = bitor(n(neighbor_mat(:,12)==1), bitorTable(3));
n(neighbor_mat(:,11)==1) = bitor(n(neighbor_mat(:,11)==1), bitorTable(4));
n(neighbor_mat(:, 6)==1) = bitor(n(neighbor_mat(:, 6)==1), bitorTable(5));
n(neighbor_mat(:, 5)==1) = bitor(n(neighbor_mat(:, 5)==1), bitorTable(6));
n(neighbor_mat(:,15)==1) = bitor(n(neighbor_mat(:,15)==1), bitorTable(7));
is_euler_inv = is_euler_inv + LUT(n);

is_euler_inv = (is_euler_inv==0);

end

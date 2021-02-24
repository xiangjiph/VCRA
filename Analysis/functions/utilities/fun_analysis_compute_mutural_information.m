function mutural_info = fun_analysis_compute_mutural_information(X, Y, num_edge_vec)
% fun_analysis_compute_mutural_information computes the mutural information
% between two array X and Y. The probability distribution is estimated by
% histcounts and histcount. The edges of bins obtained from the joint
% probabiliy distribution is used for computing the probability
% distribution of X and Y respectively. 
% Input: 
%   X, Y: numerical array of the same size. 
%   num_edge_vec: 2-by-1 numerical vector, specify the number of edges used for the 
%   estimating the probabiliy
%   distribution.
% Output: 
%   mutural_info: numerical scalar, mutural information of the probability
%   distribution of X and Y. 
%
% Implemented by Xiang Ji on May 30, 2019
if nargin < 3
    num_edge_vec = [];
end
if ~isfloat(X)
    X = double(X(:));
    Y = double(Y(:));
end
if isempty(num_edge_vec)
    [prob_xy, x_edge, y_edge] = histcounts2(X, Y, 'Normalization', 'pdf');
else
    [prob_xy, x_edge, y_edge] = histcounts2(X, Y, num_edge_vec, 'Normalization', 'pdf');
end
x_width = diff(x_edge);
y_width = diff(y_edge);
prob_x = y_width * (prob_xy.');
prob_y = x_width * prob_xy;
assert(isrow(prob_x) && isrow(prob_y));
mutural_info = sum(prob_xy .* log2(prob_xy ./ (prob_x' * prob_y)) .* (x_width.' * y_width), 'all', 'omitnan');
end
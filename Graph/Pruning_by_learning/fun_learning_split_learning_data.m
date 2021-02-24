function [training_set, test_set] = fun_learning_split_learning_data(data, label, training_set_ratio, method)
% fun_learning_split_learning_data split the data into training set and
% validation set. If the traning data are unbalanced, create a balanced
% training set and set the rest to be validation set. 
% Input: 
%   data: table or matrix, where each row is the feature vector of the
%   data point. 
%   label: numerical or logical vector, whose length is the same as the
%   number of row in data
%   training_set_ratio: ratio of the data for training, scalar between 0
%   and 1
%
% Output: 
%   training_set: structure with field data and label
%   test_set: structure with field data and label 
% 

if nargin < 3
    training_set_ratio = 0.8;
    method = 'balance';
elseif nargin < 4
    method = 'balance';
end
if isa(data, 'table')
    data = table2array(data);
end
switch method
    case 'balance'
        true_idx = find(label);
        false_idx = find(~label);
        num_true = numel(true_idx);
        num_false = numel(false_idx);
        min_num = min([num_true, num_false]);
        % Balance training set
        num_true_training = round(training_set_ratio * min_num);
        % Random sample
        training_true_idx = randsample(true_idx, num_true_training, false);
        training_false_idx = randsample(false_idx, num_true_training, false);
        % Seperating data
        is_training_set_Q = false(size(label));
        is_training_set_Q([training_true_idx; training_false_idx]) = true;
    case 'random'
        num_data = numel(label);
        num_training = round(training_set_ratio * num_data);
        training_idx = randsample(1:num_data, num_training, false);
        is_training_set_Q = false(num_data, 1);
        is_training_set_Q(training_idx) = true;
    otherwise
        error('Unrecognized method');
end
training_set.data = data(is_training_set_Q, :);
training_set.label = label(is_training_set_Q);
training_set.ind = find(is_training_set_Q);

test_set.data = data(~is_training_set_Q, :);
test_set.label = label(~is_training_set_Q);
test_set.ind = find(~is_training_set_Q);
end
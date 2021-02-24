function stat = fun_learning_get_validation_statistics(predicted_label, ground_true_label)
% fun_learning_get_validation_statistics computes the statistics for the
% validation set by comparing the predicted label to the ground true label.
% Input: 
%   predicted_label: logical vector, list of label prediceted by the
%   classifier
%   ground_true_label: logical vector of the same size as the predicted
%   label. List of the ground truth label
% Output: 
%   stat: structure with fields defined below
% 
% 
stat.true_positive_Q = (ground_true_label == true & predicted_label == true);
stat.true_negative_Q = (ground_true_label == false & predicted_label == false);
stat.false_positive_Q = (ground_true_label == false & predicted_label == true);
stat.false_negative_Q = (ground_true_label == true & predicted_label == false);

stat.num_true_positive = nnz(stat.true_positive_Q);
stat.num_predicted_positive = nnz(predicted_label == true);
stat.num_label_positive = nnz(ground_true_label == true);

stat.num_true_negative = nnz(stat.true_negative_Q);
stat.num_predicted_negative = nnz(predicted_label == false);
stat.num_label_negative = nnz(ground_true_label == false);

stat.num_false_positive = nnz(stat.false_positive_Q );
stat.num_false_negative = nnz(stat.false_negative_Q);

stat.true_positive = stat.num_true_positive / stat.num_predicted_positive;
stat.true_negative = stat.num_true_negative / stat.num_predicted_negative;
stat.false_positive = stat.num_false_positive / stat.num_predicted_positive;
stat.false_negative = stat.num_false_negative / stat.num_predicted_negative;

stat.recall = stat.num_true_positive / ...
    (stat.num_true_positive + stat.num_false_negative);

stat.precision = stat.num_true_positive / ...
    (stat.num_true_positive + stat.num_false_positive);
stat.confusion_matrix_fig = plotconfusion(ground_true_label', predicted_label');
stat.confusion_matrix_cdata = getframe(stat.confusion_matrix_fig);




end
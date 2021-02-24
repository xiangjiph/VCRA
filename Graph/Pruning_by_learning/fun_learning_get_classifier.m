function classifier = fun_learning_get_classifier(data, feature_name, classifier_name)
% This is just a draft for facilitating checking the performance of the
% features. 
if nargin < 3
    classifier_name = 'AdaBoostM1';
end
% Training parameters
training_set_ratio = 0.80;
sample_method = 'random';
%%
classifier = struct;
classifier.data = data;
% Normalize the intensity data
classifier.data.feature_name = classifier.data.features.Properties.VariableNames;
% for iter_feature = 1 : numel(classifier.data.feature_name)
%     fprintf('%d\t%s\n', iter_feature, classifier.data.feature_name{iter_feature});
% end
% Use background and mask intensity information
% classifier.used_feature_idx = [2, 4, 7, 9, 20, 26, 28, 30, 31, 32, 14, 15, 16, 17];
% Without background or intensity information
classifier.used_feature_name = feature_name;

classifier.num_example = size(classifier.data.features, 1);

[classifier.training_set, classifier.validation_set] = fun_learning_split_learning_data(classifier.data.features(:, classifier.used_feature_name), ...
    classifier.data.label, training_set_ratio, sample_method);
% For hyperparameter tuning, use Learning_hyperparameter_tuning
switch classifier_name
    case 'Bag'
        templS = templateTree('MaxNumSplits', 256, 'Surrogate', 'On');
%         templS = templateTree('Surrogate', 'On');
        classifier.classifier = fitcensemble(classifier.training_set.data, classifier.training_set.label, ...
            'Method','Bag','NumLearningCycles',100, 'Learners',templS);        
    case 'AdaBoostM1'
        classifier.classifier = fitcensemble(classifier.training_set.data, classifier.training_set.label, 'Method', 'AdaBoostM1');
    case 'SVM-gaussian'
        classifier.classifier = fitcsvm(classifier.training_set.data, classifier.training_set.label, 'KernelFunction', 'gaussian');
end
% Test classifier - too few negative sample
test_label = classifier.validation_set.label;
predicted_label = classifier.classifier.predict(classifier.validation_set.data);
% Validation
classification_stat = fun_learning_get_validation_statistics(predicted_label, test_label);

classification_stat.false_positive_idx = classifier.validation_set.ind(classification_stat.false_positive_Q);
classification_stat.false_negative_idx = classifier.validation_set.ind(classification_stat.false_negative_Q);

classifier.validation_set.stat = classification_stat;
end
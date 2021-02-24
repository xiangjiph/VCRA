set_env;
DataManager = FileManager;
dataset_name = 'WholeBrain';
stack = 'ML_2018_08_15';
grid_c_version = '240_cube_combined_5_o_1';
grid_c = DataManager.load_grid(dataset_name, stack, grid_c_version);
grid_info = grid_c.grid_ori;
gpuDevice(2);
%% Load training data
% 1. Block 150 and 151 are from the imresize3 downsampled dataset. The rest
% of the training data are from the max-pooling downsampling. 
% 2. 
combined_grid_idx_list = [219 220 306 338 873];
num_training_block = numel(combined_grid_idx_list);
training_data = cell(num_training_block, 1);
for iter_block = 1 : num_training_block
    fn_training_data = DataManager.fp_training_data(dataset_name, stack, ...
        sprintf('%s_%d', grid_c_version, combined_grid_idx_list(iter_block)));
    training_data{iter_block} = load(fn_training_data);
end
used_training_data_idx = 1:5;
num_block_used = numel(used_training_data_idx);
%% Merge data - normal link_ep1
link_ep1 = struct;
link_ep1.features = [];
link_ep1.label = [];
for iter_block = 1 : num_block_used
    tmp_cell = training_data{iter_block}.link_ep1.normal;
    for iter_cell = 1 : numel(tmp_cell)
        tmp_str = tmp_cell{iter_cell};
        link_ep1.features = cat(1, link_ep1.features, tmp_str.features);
        link_ep1.label = cat(1, link_ep1.label, tmp_str.label);
    end
end
link_ep1.features_mean = mean(table2array(link_ep1.features), 1, 'omitnan');
link_ep1.features_std = std(table2array(link_ep1.features), 1, 'omitnan');
%%
fun_learning_print_feature_name_list(classifier_link_ep1.data)
classifier_link_ep1 = fun_learning_get_classifier(link_ep1, [2, 4, 7, 10, 21, 27, 29, 31], 'Bag');
classifier_link_ep1.validation_set.stat
%%
classifier_link_ep1 = struct;
classifier_link_ep1.data = link_ep1;
classifier_link_ep1.used_feature_idx = ;
classifier_link_ep1.used_feature_name = classifier_link_ep1.data.features.Properties.VariableNames(classifier_link_ep1.used_feature_idx);

classifier_link_ep1.num_example = size(classifier_link_ep1.data.features, 1);
% Baised dataset. Find a balance one. 
true_idx = find(classifier_link_ep1.data.label);
false_idx = find(~classifier_link_ep1.data.label);
num_true_training = round(0.75 * nnz(true_idx));

num_trial = 100;
clearvars('stat_record');
% stat_record = struct;
for iter_trial = num_trial : -1 : 1
    training_true_idx = randsample(true_idx, num_true_training, false);
    training_false_idx = randsample(false_idx, num_true_training, false);
    tmp_training_set_Q = false(classifier_link_ep1.num_example, 1);
    tmp_training_set_Q([training_true_idx; training_false_idx]) = true;
    % tmp_training_set_Q = rand(classifier_link_ep1.num_example, 1) < 0.5;
    classifier_link_ep1.training_set.data = table2array(classifier_link_ep1.data.features(tmp_training_set_Q, classifier_link_ep1.used_feature_idx));
    classifier_link_ep1.training_set.label = classifier_link_ep1.data.label(tmp_training_set_Q);

    classifier_link_ep1.validation_set.data = table2array(classifier_link_ep1.data.features(~tmp_training_set_Q, classifier_link_ep1.used_feature_idx));
    classifier_link_ep1.validation_set.label = classifier_link_ep1.data.label(~tmp_training_set_Q);

    classifier_link_ep1.classifier = fitcensemble(classifier_link_ep1.training_set.data, classifier_link_ep1.training_set.label, 'Method', 'AdaBoostM1');
    % classifier_link_ep1.classifier = fitcsvm(classifier_link_ep1.training_set.data, classifier_link_ep1.training_set.label, 'KernelFunction', 'gaussian');
    % plot(kfoldLoss(classifier_link_ep1.classifier, 'Mode', 'cumulative'));
    % Test classifier - too few negative sample
    test_label = classifier_link_ep1.validation_set.label;
    predicted_label = classifier_link_ep1.classifier.predict(classifier_link_ep1.validation_set.data);
    classification_stat.num_true_positive = nnz(test_label == true & predicted_label == true);
    classification_stat.num_positive = nnz(predicted_label == true);
    classification_stat.num_true_negative = nnz(test_label == false & predicted_label == false);
    classification_stat.num_negative = nnz(predicted_label == false);
    classification_stat.num_false_positive = nnz(test_label == false & predicted_label == true);
    classification_stat.num_false_negative = nnz(test_label == true & predicted_label == false);

    classification_stat.true_positive = classification_stat.num_true_positive / classification_stat.num_positive;
    classification_stat.true_negative = classification_stat.num_true_negative / classification_stat.num_negative;
    classification_stat.false_positive = classification_stat.num_false_positive / classification_stat.num_positive;
    classification_stat.false_negative = classification_stat.num_false_negative / classification_stat.num_negative;

    classification_stat.recall = classification_stat.num_true_positive / ...
        (classification_stat.num_true_positive + classification_stat.num_false_negative);

    classification_stat.precision = classification_stat.num_true_positive / ...
        (classification_stat.num_true_positive + classification_stat.num_true_negative);
    stat_record(iter_trial) = classification_stat;
end

classifier_link_ep1.validation_set.stat = classification_stat;

classifier_link_ep1.fp = fullfile('./Graph/Pruning_by_learning/', sprintf('classifier_link_ep1_to_remove_%s.mat', datetime('today')));
save(classifier_link_ep1.fp , '-struct', 'classifier_link_ep1');








%% Merge data - Linker
linker_data = struct;
linker_data.data = {};
linker_data.label.removed = {};
linker_data.label.not_sure = {};
num_round_data = 0;
for idx = 1 : num_block_used
    tmp_data = training_data{used_training_data_idx(idx)}.linker;
    tmp_round_name = fieldnames(tmp_data);
    tmp_num_round = numel(tmp_round_name);
    for iter_round = 1 : tmp_num_round
        num_round_data = num_round_data + 1;
        tmp_round_data = tmp_data.(tmp_round_name{iter_round});
        linker_data.data{num_round_data} = tmp_round_data.data;
        linker_data.label.removed{num_round_data} = tmp_round_data.label.to_removeQ;
        linker_data.label.not_sure{num_round_data} = tmp_round_data.label.not_sureQ;
    end
end
linker_data.data = cat(1, linker_data.data{:});
linker_data.label.removed = cat(1, linker_data.label.removed{:});
linker_data.label.not_sure = cat(1, linker_data.label.not_sure{:});
%% Merge data - erroneous link
error_link_data = struct;
error_link_data.data = {};
error_link_data.label.removed = {};
error_link_data.label.not_sure = {};
num_round_data = 0;
for idx = 1 : num_block_used
    tmp_data = training_data{used_training_data_idx(idx)}.fake_link;
    tmp_round_name = fieldnames(tmp_data);
    tmp_num_round = numel(tmp_round_name);
    for iter_round = 1 : tmp_num_round
        num_round_data = num_round_data + 1;
        tmp_round_data = tmp_data.(tmp_round_name{iter_round});
        error_link_data.data{num_round_data} = tmp_round_data.features;
        error_link_data.label.removed{num_round_data} = tmp_round_data.label.to_removeQ;
        error_link_data.label.not_sure{num_round_data} = tmp_round_data.label.not_sureQ;
    end
end
error_link_data.data = cat(1, error_link_data.data{:});
error_link_data.label.removed = cat(1, error_link_data.label.removed{:});
error_link_data.label.not_sure = cat(1, error_link_data.label.not_sure{:});
%% Merge data - link with 1 endpoint 
link_ep1_data = struct;
link_ep1_data.data = {};
link_ep1_data.label.removed = {};
link_ep1_data.label.not_sure = {};
num_round_data = 0;
for idx = 1 : num_block_used
    tmp_data = training_data{used_training_data_idx(idx)}.link_ep1;
    tmp_round_name = fieldnames(tmp_data);
    tmp_num_round = numel(tmp_round_name);
    for iter_round = 1 : tmp_num_round
        num_round_data = num_round_data + 1;
        tmp_round_data = tmp_data.(tmp_round_name{iter_round});
        if isfield(tmp_round_data, 'data')
            link_ep1_data.data{num_round_data} = struct2table(table2struct(tmp_round_data.data));        
        elseif isfield(tmp_round_data, 'features')
            link_ep1_data.data{num_round_data} = struct2table(table2struct(tmp_round_data.features));        
        end
        link_ep1_data.label.removed{num_round_data} = tmp_round_data.label.to_removeQ;
        link_ep1_data.label.not_sure{num_round_data} = tmp_round_data.label.not_sureQ;
    end
end
link_ep1_data.data = cat(1, link_ep1_data.data{:});
link_ep1_data.label.removed = cat(1, link_ep1_data.label.removed{:});
link_ep1_data.label.not_sure = cat(1, link_ep1_data.label.not_sure{:});
%% Linker - Too few positive training set
% Conclusion: need to further improve the feature desgned for linker. For
% example, need to incorporate the feature from the link from which the
% linker is found. Also, maybe the target endpoint distance transform is
% also important. 
clearvars classifier_linker
classifier_linker = struct;
classifier_linker.data = linker_data;
% Normalize the intensity data
int_max = 65535;
classifier_linker.data.data.int_mean_n = classifier_linker.data.data.int_mean ./ int_max;
classifier_linker.data.data.int_om_mean_n = classifier_linker.data.data.int_om_mean./ int_max;

classifier_linker.data.data.int_std_n = classifier_linker.data.data.int_std ./ int_max;
classifier_linker.data.data.int_om_std_n = classifier_linker.data.data.int_om_std./ int_max;

classifier_linker.data.feature_name = classifier_linker.data.data.Properties.VariableNames;
for iter_feature = 1 : numel(classifier_linker.data.feature_name)
    fprintf('%d\t%s\n', iter_feature, classifier_linker.data.feature_name{iter_feature});
end
classifier_linker.used_feature_idx = [1, 3, 8, 11];
classifier_linker.used_feature_name = classifier_linker.data.feature_name(classifier_linker.used_feature_idx);

classifier_linker.num_example = size(classifier_linker.data.data, 1);
% Baised dataset. Find a balance one. 
true_idx = find(classifier_linker.data.label.removed);
false_idx = find(~classifier_linker.data.label.removed);
num_true_training = round(0.50 * nnz(true_idx));
training_true_idx = randsample(true_idx, num_true_training, false);
training_false_idx = randsample(false_idx, num_true_training, false);
tmp_training_set_Q = false(classifier_linker.num_example, 1);
tmp_training_set_Q([training_true_idx; training_false_idx]) = true;
% tmp_training_set_Q = rand(classifier_linker.num_example, 1) < 0.5;
classifier_linker.training_set.data = table2array(classifier_linker.data.data(tmp_training_set_Q, classifier_linker.used_feature_idx));
classifier_linker.training_set.label = classifier_linker.data.label.removed(tmp_training_set_Q);

classifier_linker.validation_set.data = table2array(classifier_linker.data.data(~tmp_training_set_Q, classifier_linker.used_feature_idx));
classifier_linker.validation_set.label = classifier_linker.data.label.removed(~tmp_training_set_Q);

% classifier_linker.classifier = fitcensemble(classifier_linker.training_set.data, classifier_linker.training_set.label, 'Method', 'LogitBoost');
% classifier_linker.classifier = fitcsvm(classifier_linker.training_set.data, classifier_linker.training_set.label, 'KernelFunction', 'gaussian');
% plot(kfoldLoss(classifier_linker.classifier, 'Mode', 'cumulative'));
% Test classifier - too few negative sample
test_label = classifier_linker.validation_set.label;
predicted_label = classifier_linker.classifier.predict(classifier_linker.validation_set.data);
classification_stat.num_true_positive = nnz(test_label == true & predicted_label == true);
classification_stat.num_positive = nnz(predicted_label == true);
classification_stat.num_true_negative = nnz(test_label == false & predicted_label == false);
classification_stat.num_negative = nnz(predicted_label == false);
classification_stat.num_false_positive = nnz(test_label == false & predicted_label == true);
classification_stat.num_false_negative = nnz(test_label == true & predicted_label == false);

classification_stat.true_positive = classification_stat.num_true_positive / classification_stat.num_positive;
classification_stat.true_negative = classification_stat.num_true_negative / classification_stat.num_negative;
classification_stat.false_positive = classification_stat.num_false_positive / classification_stat.num_positive;
classification_stat.false_negative = classification_stat.num_false_negative / classification_stat.num_negative;
classifier_linker.validation_set.stat = classification_stat;
classifier_linker.fp = fullfile('./Graph/Pruning_by_learning/', sprintf('classifier_linker_to_remove_%s.mat', datetime('today')));
save(classifier_linker.fp , '-struct', 'classifier_linker');
%% Plot
tmp_feature = 'linker_SNR';
figure;
histogram(classifier_linker.data.data.(tmp_feature)(classifier_linker.data.label.removed));
hold on 
histogram(classifier_linker.data.data.(tmp_feature)(~classifier_linker.data.label.removed));
legend('Removed', 'Kept');
hold off
%% Link with one endpoint
clearvars classifier_link_ep1
classifier_link_ep1 = struct;
classifier_link_ep1.data = link_ep1_data;
% Normalize the intensity data
int_max = 65535;
classifier_link_ep1.data.data.int_mean_n = classifier_link_ep1.data.data.int_mean ./ int_max;
classifier_link_ep1.data.data.int_median_n = classifier_link_ep1.data.data.int_median./ int_max;
classifier_link_ep1.data.feature_name = classifier_link_ep1.data.data.Properties.VariableNames;
for iter_feature = 1 : numel(classifier_link_ep1.data.feature_name)
    fprintf('%d\t%s\n', iter_feature, classifier_link_ep1.data.feature_name{iter_feature});
end
classifier_link_ep1.used_feature_idx = [2, 4, 7, 15, 16, 22, 24, 26, 27, 28, 29];
classifier_link_ep1.used_feature_name = classifier_link_ep1.data.feature_name(classifier_link_ep1.used_feature_idx);

classifier_link_ep1.num_example = size(classifier_link_ep1.data.data, 1);
% Baised dataset. Find a balance one. 
true_idx = find(classifier_link_ep1.data.label.removed);
false_idx = find(~classifier_link_ep1.data.label.removed);
num_true_training = round(0.75 * nnz(true_idx));
training_true_idx = randsample(true_idx, num_true_training, false);
training_false_idx = randsample(false_idx, num_true_training, false);
tmp_training_set_Q = false(classifier_link_ep1.num_example, 1);
tmp_training_set_Q([training_true_idx; training_false_idx]) = true;
% tmp_training_set_Q = rand(classifier_link_ep1.num_example, 1) < 0.5;
classifier_link_ep1.training_set.data = table2array(classifier_link_ep1.data.data(tmp_training_set_Q, classifier_link_ep1.used_feature_idx));
classifier_link_ep1.training_set.label = classifier_link_ep1.data.label.removed(tmp_training_set_Q);

classifier_link_ep1.validation_set.data = table2array(classifier_link_ep1.data.data(~tmp_training_set_Q, classifier_link_ep1.used_feature_idx));
classifier_link_ep1.validation_set.label = classifier_link_ep1.data.label.removed(~tmp_training_set_Q);

classifier_link_ep1.classifier = fitcensemble(classifier_link_ep1.training_set.data, classifier_link_ep1.training_set.label, 'Method', 'AdaBoostM1');
% classifier_link_ep1.classifier = fitcsvm(classifier_link_ep1.training_set.data, classifier_link_ep1.training_set.label, 'KernelFunction', 'gaussian');
% plot(kfoldLoss(classifier_link_ep1.classifier, 'Mode', 'cumulative'));
% Test classifier - too few negative sample
test_label = classifier_link_ep1.validation_set.label;
predicted_label = classifier_link_ep1.classifier.predict(classifier_link_ep1.validation_set.data);
classification_stat.num_true_positive = nnz(test_label == true & predicted_label == true);
classification_stat.num_positive = nnz(predicted_label == true);
classification_stat.num_true_negative = nnz(test_label == false & predicted_label == false);
classification_stat.num_negative = nnz(predicted_label == false);
classification_stat.num_false_positive = nnz(test_label == false & predicted_label == true);
classification_stat.num_false_negative = nnz(test_label == true & predicted_label == false);

classification_stat.true_positive = classification_stat.num_true_positive / classification_stat.num_positive;
classification_stat.true_negative = classification_stat.num_true_negative / classification_stat.num_negative;
classification_stat.false_positive = classification_stat.num_false_positive / classification_stat.num_positive;
classification_stat.false_negative = classification_stat.num_false_negative / classification_stat.num_negative;

classification_stat.recall = classification_stat.num_true_positive / ...
    (classification_stat.num_true_positive + classification_stat.num_false_negative);

classification_stat.precision = classification_stat.num_true_positive / ...
    (classification_stat.num_true_positive + classification_stat.num_true_negative);


classifier_link_ep1.validation_set.stat = classification_stat;

classifier_link_ep1.fp = fullfile('./Graph/Pruning_by_learning/', sprintf('classifier_link_ep1_to_remove_%s.mat', datetime('today')));
save(classifier_link_ep1.fp , '-struct', 'classifier_link_ep1');
%% 
clearvars classifier_error_link
classifier_error_link = struct;
classifier_error_link.data = error_link_data;
% Normalize the intensity data
int_max = 65535;
classifier_error_link.data.data.int_mean_n = classifier_error_link.data.data.int_mean ./ int_max;
classifier_error_link.data.data.int_median_n = classifier_error_link.data.data.int_median ./ int_max;
classifier_error_link.data.data.int_middle_point_n = classifier_error_link.data.data.int_middle_point ./ int_max;
classifier_error_link.data.feature_name = classifier_error_link.data.data.Properties.VariableNames;
for iter_feature = 1 : numel(classifier_error_link.data.feature_name)
    fprintf('%d\t%s\n', iter_feature, classifier_error_link.data.feature_name{iter_feature});
end
classifier_error_link.used_feature_idx = [2, 52, 17, 30, 31, 33, 36, 38, 40, 41, 44, 46, 48, 49, 50];
classifier_error_link.used_feature_name = classifier_error_link.data.feature_name(classifier_error_link.used_feature_idx);

classifier_error_link.num_example = size(classifier_error_link.data.data, 1);
% Baised dataset. Find a balance one. 
true_idx = find(classifier_error_link.data.label.removed);
false_idx = find(~classifier_error_link.data.label.removed);
num_true_training = round(0.75 * nnz(true_idx));
training_true_idx = randsample(true_idx, num_true_training, false);
training_false_idx = randsample(false_idx, num_true_training, false);
tmp_training_set_Q = false(classifier_error_link.num_example, 1);
tmp_training_set_Q([training_true_idx; training_false_idx]) = true;
validation_idx = find(~tmp_training_set_Q);
% tmp_training_set_Q = rand(classifier_error_link.num_example, 1) < 0.5;
classifier_error_link.training_set.data = table2array(classifier_error_link.data.data(tmp_training_set_Q, classifier_error_link.used_feature_idx));
classifier_error_link.training_set.label = classifier_error_link.data.label.removed(tmp_training_set_Q);

classifier_error_link.validation_set.data = table2array(classifier_error_link.data.data(~tmp_training_set_Q, classifier_error_link.used_feature_idx));
classifier_error_link.validation_set.label = classifier_error_link.data.label.removed(~tmp_training_set_Q);

classifier_error_link.classifier = fitcensemble(classifier_error_link.training_set.data, classifier_error_link.training_set.label, 'Method', 'AdaBoostM1');
% classifier_error_link.classifier = fitcsvm(classifier_error_link.training_set.data, classifier_error_link.training_set.label, 'KernelFunction', 'gaussian');
% plot(kfoldLoss(classifier_error_link.classifier, 'Mode', 'cumulative'));
% Test classifier - too few negative sample
test_label = classifier_error_link.validation_set.label;
predicted_label = classifier_error_link.classifier.predict(classifier_error_link.validation_set.data);
classification_stat.num_true_positive = nnz(test_label == true & predicted_label == true);
classification_stat.num_positive = nnz(predicted_label == true);
classification_stat.num_true_negative = nnz(test_label == false & predicted_label == false);
classification_stat.num_negative = nnz(predicted_label == false);
classification_stat.num_false_positive = nnz(test_label == false & predicted_label == true);
classification_stat.num_false_negative = nnz(test_label == true & predicted_label == false);

classification_stat.true_positive = classification_stat.num_true_positive / classification_stat.num_positive;
classification_stat.true_negative = classification_stat.num_true_negative / classification_stat.num_negative;
classification_stat.false_positive = classification_stat.num_false_positive / classification_stat.num_positive;
classification_stat.false_negative = classification_stat.num_false_negative / classification_stat.num_negative;
classifier_error_link.validation_set.stat = classification_stat;

classifier_error_link.fp = fullfile('./Graph/Pruning_by_learning/', sprintf('classifier_error_link_to_remove_%s.mat', datetime('today')));
save(classifier_error_link.fp , '-struct', 'classifier_error_link');
%% Look at the false positive case
tmp_false_positive_Q = (test_label == false & predicted_label == true);
tmp_false_positive_idx = validation_idx(tmp_false_positive_Q);
classifier_error_link.data.data(tmp_false_positive_idx, :)
% Conclusion:
% 1. You want the classifier to do too much. One classifier should
% classifier one type of erroneous connections. Focus on the one that is
% due to segmentation error first. Think more about how you make the
% decision on which on to delete. 
% 2. Record the erroneous link that are due to sample preparation later.
% Those links are very hard to be distinguish from normal links. 
% 3. Also, consider the possibility to train a neural network that works
% directly on the local images, just like what I am doing in practice, to
% detect link for delection. 
% 4. Currently the false positive rate is too high. I want to be
% conservative in this part. 
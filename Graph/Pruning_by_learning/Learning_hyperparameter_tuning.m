% Prepare data
training_set_ratio = 0.90;
rng(1234)
selected_feature_name = fun_learning_print_feature_name_list(link_ep1, [2, 4, 7, 10, 12, 21, 23, 25, 27, 29, 31]);
[classifier.training_set, classifier.validation_set] = fun_learning_split_learning_data(link_ep1.features(:, selected_feature_name), ...
    link_ep1.label, training_set_ratio, 'random');
% Training set
X = classifier.training_set.data;
Y = classifier.training_set.label;
%% Train a 3 layers neural network for classification 

    %% Tree + Bag
max_num_splits = 2 .^ (4 : 2 : 8);
method_name = {'Bag', 'LogitBoost'};
num_method = numel(method_name);
model_cell = cell(num_method, numel(max_num_splits));
numTrees = 100;
cross_validation_fold = 10;
for iter_method = 1 : num_method
    for iter_split = 1 : numel(max_num_splits)
        t = templateTree('MaxNumSplits', max_num_splits(iter_split), 'Surrogate','on');
        switch iter_method
            case 1
                model_cell{iter_method, iter_split} = fitcensemble(X, Y, 'Method','Bag',...
                    'NumLearningCycles',numTrees,'Learners',t,...
                    'KFold',cross_validation_fold);
%             case 2
%                 model_cell{iter_method, iter_split} = fitcensemble(X, Y, 'Method','AdaBoostM1',...
%                     'NumLearningCycles',numTrees,'Learners',t,...
%                     'KFold',cross_validation_fold, 'LearnRate', 0.1);
            case 2
                model_cell{iter_method, iter_split} = fitcensemble(X, Y, 'Method','LogitBoost',...
                    'NumLearningCycles',numTrees,'Learners',t,...
                    'KFold',cross_validation_fold, 'LearnRate', 0.1);
        end
    end
end
%% Visualize cross-validation error for model comparison
plot_opt = struct;
plot_opt.FileType = 'png';
plot_opt.Output_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Classifier');

errorCell = cellfun(@(x)kfoldLoss(x, 'Mode', 'cumulative') , model_cell, 'Uniform', false);
for iter_method = 1 : numel(method_name)
    fig_handle = figure;
    for iter_split = 1 : numel(max_num_splits)
        plot(errorCell{iter_method, iter_split}, 'LineWidth', 2);
        hold on;
    end
    xlabel('Number of trees');
    ylabel('Cross-validated MSE');
    legend([cellstr(num2str(max_num_splits','Maximum number of splits= %d'))]);
    title(sprintf('Training method: %s', method_name{iter_method}));
    hold off;
end
% Conclusion: Bag seems to work better
%% Plot cross-validation error
% test_classifier = fitcensemble(X, Y, 'Method','LogitBoost',...
%                     'NumLearningCycles',numTrees,'Learners',t,...
%                     'LearnRate', 0.1);
test_classifier_cv = fitcensemble(X, Y, 'Method','Bag',...
    'NumLearningCycles',numTrees,'Learners',t, 'KFold', cross_validation_fold);
test_classifier = fitcensemble(X, Y, 'Method','Bag',...
    'NumLearningCycles',numTrees,'Learners',t);
%
test_cxval_err = kfoldLoss(test_classifier_cv, 'Mode', 'cumulative');
fig_handle = figure;
ax = axes(fig_handle);
plot(ax, test_cxval_err, 'LineWidth', 2)
xlabel(ax, 'Number of trees');
ylabel(ax, 'Cross-validation MSE');
ax.LineWidth = 2;
ax.FontSize = 14;
grid on
plot_opt = struct;
plot_opt.FileType = 'png';
plot_opt.Output_folder = fullfile(DataManager.fp_visualization_folder(dataset_name, stack), 'Classifier');
plot_opt.FileName = ['Classifier_link_ep1_bag_10_fold_cross_validation_error.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image(fig_handle, plot_opt);
%% Plot confusion matrix
%%
test_label = classifier.validation_set.label;
predicted_label = test_classifier.predict(classifier.validation_set.data);
% Validation
classification_stat = fun_learning_get_validation_statistics(predicted_label, test_label);

plot_opt.FileName = ['Classifier_link_ep1_bag_confusion_matrix_random_sample.', plot_opt.FileType];
plot_opt.FilePath = fullfile(plot_opt.Output_folder, plot_opt.FileName);
fun_print_image(classification_stat.fig_confusion_matrix, plot_opt);
%% Profile performance of Tree-ensemble classifier
MdlDeep = fitrtree(X,Y,'CrossVal','on','MergeLeaves','off',...
    'MinParentSize',1,'Surrogate','on');
MdlStump = fitrtree(X,Y,'MaxNumSplits',1,'CrossVal','on','Surrogate','on');

n = size(X,1);
m = max(8, floor(log2(n - 1)));
lr = [0.05 0.1 0.2 0.3, 0.4, 0.5];
% maxNumSplits = 2.^(0:m);
maxNumSplits = 2:2:16;

Mdl = cell(numel(maxNumSplits),numel(lr));
rng(1); % For reproducibility
for k = 1:numel(lr)
    fprintf('Testing learning rate %f\n', lr(k));
    for j = 1:numel(maxNumSplits)
        
%         Mdl{j,k} = fitrensemble(X,Y,'Method','LSBoost',...
%             'NumLearningCycles',numTrees,'Learners',t,...
%             'KFold',5,'LearnRate',lr(k));

    end
end
%%
kflAll = @(x)kfoldLoss(x,'Mode','cumulative');
errorCell = cellfun(kflAll,Mdl,'Uniform',false);
error = reshape(cell2mat(errorCell),[numTrees numel(maxNumSplits) numel(lr)]);
errorDeep = kfoldLoss(MdlDeep);
errorStump = kfoldLoss(MdlStump);
%%
mnsPlot = [1 round(numel(maxNumSplits)/2) numel(maxNumSplits)];
figure;
for k = 1:3
    subplot(2,2,k);  
    plot(squeeze(error(:,mnsPlot(k),:)),'LineWidth',2);
    set(gca, 'YScale', 'log');
%     axis tight;
    hold on;
    h = gca;
    plot(h.XLim,[errorDeep errorDeep],'-.b','LineWidth',2);
    plot(h.XLim,[errorStump errorStump],'-.r','LineWidth',2);
    plot(h.XLim,min(min(error(:,mnsPlot(k),:))).*[1 1],'--k');
%     h.YLim = [10 50];    
    xlabel 'Number of trees';
    ylabel 'Cross-validated MSE';
    title(sprintf('MaxNumSplits = %0.3g', maxNumSplits(mnsPlot(k))));
    hold off;
end
hL = legend([cellstr(num2str(lr','Learning Rate = %0.2f'));...
        'Deep Tree';'Stump';'Min. MSE']);
hL.Position(1) = 0.6;
%%
tmp_classifier = fitrensemble(X, Y, 'Method','Bag',...
            'NumLearningCycles',numTrees,'Learners',t, 'OptimizeHyperparameters', 'all',...
            'HyperparameterOptimizationOptions', struct('MaxObjectiveEvaluations', {100}, 'Verbose', 2));
%% Train a 3 layer neural network for classification 

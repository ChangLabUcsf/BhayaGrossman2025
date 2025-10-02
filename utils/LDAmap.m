% performs LDA with kfold cross-validation and returns the data matrix
% mapped onto 2-dimensional LD space across runs
function [y_pred, mappedX, MdlLinear, fold_accuracies] = LDAmap(X, y, kfolds) 

    MdlLinear = fitcdiscr(X, y, ...
        'DiscrimType', 'linear', 'KFold', ...
        kfolds, 'Prior','uniform');
    dims = min([3, size(X, 2)]);
    mappedX = nan(kfolds, length(y), dims);
    allW = nan(kfolds, size(X, 2), dims);
    
    % Preallocate array for fold accuracies
    fold_accuracies = zeros(kfolds, 1);

    for k = 1:length(MdlLinear.Trained)   
        % Get training and test indices for this fold
        rng(50);
        test_idx = MdlLinear.Partition.test(k);
        
        % Get predictions for this fold
        fold_predictions = predict(MdlLinear.Trained{k}, X(test_idx,:));
        
        % Calculate accuracy for this fold
        fold_accuracies(k) = sum(fold_predictions == y(test_idx)) / sum(test_idx);
        
        [W, LAMBDA] = eig(MdlLinear.Trained{k}.BetweenSigma, ...
            MdlLinear.Trained{k}.Sigma); 
        lambda = diag(LAMBDA);
        [~, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        allLDs = X*W;
        mappedX(k, :, :) = allLDs(:, 1:dims);
        % return corresponding weights for LDs
        allW(k, :, :) = W(:, 1:dims);
    end

    % find average mapping between trained models
    y_pred = kfoldPredict(MdlLinear);
    
    % alternative community built alternative
    % [mappedX, mapping]=lda(coeff(:, 1:n_comp), y);

    %[y_pred,~, ~] = predict(MdlLinear, coeff(:, 1:n_comp));
    % visualize the accuracy as number of components increases
    % [acc, ~] = accOverComp(A, max_comp, y);
    % figure;
    % plot(acc, 'LineWidth', 2);
    % xlabel('Num Components');
    % ylabel('Accuracy');
    % yline(1/length(unique(y)), 'LineStyle', '--', 'Color', 'k');
    % xline(n_comp, 'LineStyle', '--', 'Color', 'r');
    % set(gca, 'FontSize', 15);
    % ylim([0 1]);
    % acs_lda = y_pred;
end


function [AUC, acc, p] = logistic_permute(X, y, pcflag, ~, tps, nfolds, nperm)
    % This function is a permutation test for the logistic regression
    % classifier. It permutes the labels of the data and calculates the
    % AUC for each permutation. The p-value is calculated as the proportion
    % of permutations that have an AUC greater than the AUC of the true
    % labels.
    if nargin < 7
        nperm = 1000;
    end

    % run pca to reduce the dimensionality of the data, run permutation test faster
    if pcflag
        [~, score, ~, ~, exp, ~] = pca(X);
        n_comp = find(diff(cumsum(exp) > 90));

        if isempty(n_comp)
            n_comp = 1;
        end

        % For reconstruction only used components retained
        pcaX = score(:, 1:n_comp);
    else
        pcaX = X;
    end

    % true AUC
    [~, ~, AUC_true] = logistic(X, y, pcflag, [], tps, nfolds);

    % preallocate
    AUC = nan(nperm, nfolds);
    acc = nan(nperm, nfolds);
    for i = 1:nperm
        % permute the labels
        y = y(randperm(length(y)));
        [~, ~, AUC(i, :), ~, ~, acc(i, :)] = ...
            logistic(pcaX, y, 0, [], tps, nfolds);
    end

    % calculate the p-value
    p = sum(AUC(:) > median(AUC_true)) / nperm;
    disp(['p-value: ' num2str(p)]);
end

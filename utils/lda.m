function [fp, tp, AUC, pcaX, scores, acc, mappedX, y_hat] = lda(X, y, pcflag, ...
    ~, ~, nfolds)
    % run linear discriminant analysis after PCA on design matrix
    % leave comp empty if pcaflag on

    if pcflag
        [~, score, ~, ~, exp] = pca(X); 
        n_comp = find(diff(cumsum(exp)>90));
    
        if isempty(n_comp)
            n_comp = 1;
        end
        
        disp(['N-comp used: ' num2str(n_comp)])
        pcaX = score(:, 1:n_comp);
    else
        pcaX = X;
    end
    
    if nargin<6, nfolds = 15; end
    acc = nan(nfolds, 1); %floor(0.2*size(pcaX,1))
    
    % test on remaining 100/nfold % of the data
    scores = nan(size(pcaX, 1), length(unique(y)));
    cv = cvpartition(size(pcaX,1),'KFold',nfolds);
    AUC = nan(nfolds, 1);
    % for mapping to LDA space
    dims = min([3, size(pcaX, 2)]);
    mappedX = nan(length(y), dims);
    allW = nan(nfolds, size(pcaX, 2), dims);
    if iscell(y)
        y_hat = cell(size(y, 1), size(y, 2));
    else
        y_hat = nan(size(y, 1), size(y, 2));
    end
    for rep = 1:nfolds          
        test_split = cv.test(rep);
        train_split = ~test_split;

        % fit model with fixed random seed
        rng(5);
        Mdl = fitcdiscr(pcaX(train_split, :), y(train_split), ...
            'DiscrimType', 'linear', 'Prior', 'uniform'); 
        % changed on 08.04.2025 from : 'Prior', 'empirical'     

        % predicted accuracy using the optimal threshold calculated for
        % this train/test split
        [ypred, test_scores] = Mdl.predict(pcaX(test_split, :));
        fp = [];
        tp = [];
        if length(unique(y))==2 % only works for binary classification
            posclass = unique(y);
            posclass = posclass(2);
            [fp,tp,~, AUC(rep), ~] = perfcurve(y(test_split),...
                test_scores(:, 2), posclass);
        else % use the rocmetrics muticlass AUC OVR average
            % will error if not all categories in test split
            try
                rocObj = rocmetrics(y(test_split),...
                    test_scores, unique(y));
                AUC(rep) = mean(rocObj.AUC);
            catch
                warning('Skipping AUC because not all categories in fold...')
            end
        end

        scores(test_split) = test_scores(:, 2);  
        if iscell(y)
            acc(rep) = mean(strcmp(y(test_split), ypred'));
        else
            acc(rep) = mean(y(test_split) == ypred');
        end
        y_hat(test_split) = ypred;

        % map                 
        [W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); 
        lambda = diag(LAMBDA);
        [~, sortord] = sort(lambda, 'descend');
        W = W(:, sortord);
        allLDs = pcaX(test_split, :)*W;
        mappedX(test_split, :) = allLDs(:, 1:dims);
        % return corresponding weights for LDs
        allW(rep, :, :) = W(:, 1:dims);
    end
end

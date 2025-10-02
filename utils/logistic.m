function [fp, tp, AUC, pcaX, scores, acc, weights, opthresh, Mdl, comp, mu] = ...
logistic(X, y, pcflag, ~, tps, nfolds)
    % This function performs logistic regression after applying Principal Component 
    % Analysis (PCA) on the design matrix over n-folds cross-validation.
    %
    % Inputs:
    %   X       - A 2D matrix of size (trials x (electrodes x time)), representing the design matrix.
    %   y       - A column vector of size (trials x 1), representing the binary response variable.
    %   n_folds - An integer specifying the number of folds for cross-validation.
    %   timepoints - A vector specifying the timepoints to use for the analysis.
    %   pcflag  - A boolean flag indicating whether to perform PCA on the design matrix.
    %   comp    - A vector specifying the number of principal components to retain. Leave empty if 'pcflag' is enabled.
    %
    % Outputs:
    %   AUC - A nfold vector representing the area under the ROC curve of the model (test set). 
    %   pcaX - The design matrix after PCA. 
    %   scores - A column vector of size (trials x 1) representing the predicted scores of the model.
    %   acc - A nfold vector representing the classification accuracy of the model (test set).
    %   weights - A 3D matrix of size (nfold x electrodes x time), representing the weights of the model.
    %   opthresh - A nfold vector representing the optimal threshold for classification.
    %   Mdl   - A struct containing the trained logistic regression model.
    %   comp  - A 2D matrix of size (electrodes x time) representing the principal components.
    %   mu    - A column vector of size (electrodes x 1) representing the mean of the design matrix.
    %   
    % Author: Ilina Bhaya-Grossman 2024
    
    if pcflag
        [comp, score, ~, ~, exp, mu] = pca(X); 
        n_comp = find(diff(cumsum(exp)>90));
    
        if isempty(n_comp)
            n_comp = 1;
        end
                
        % For reconstruction only used components retained
        comp = comp(:, 1:n_comp);
        pcaX = score(:, 1:n_comp);
    else
        comp = [];
        n_comp = size(X, 2);
        pcaX = X;
    end
    disp(['N-comp used: ' num2str(n_comp)])
    
    if nargin<6, nfolds = 20; end
    acc = nan(nfolds, 1); %floor(0.2*size(pcaX,1))
    
    % test on remaining 100/nfold of the data
    scores = nan(size(X, 1), 1);
    rng(1);
    cv = cvpartition(size(pcaX,1),'KFold',nfolds);
    AUC = nan(nfolds, 1);
    opthresh = nan(nfolds, 1);
    for rep = 1:nfolds          
        test_split = cv.test(rep);
        train_split = ~test_split;

        % fit model with fixed random seed
        Mdl = fitclinear(pcaX(train_split, :), y(train_split), ...
            "Learner","logistic", "Regularization","lasso", ...
            "ScoreTransform", 'logit');

        % find auc curve for all relevant thresholds         
        [~, train_scores] = Mdl.predict(pcaX(train_split, :));
        [fp_train, tp_train,thresh, ~, optrec] = perfcurve(y(train_split),...
            train_scores(:, 2), 1);           

        % find optimal threshold for classification
        opthresh(rep) = thresh(optrec(1) == fp_train & optrec(2) == tp_train);

        % predicted accuracy using the optimal threshold calculated for
        % this train/test split
        % ypred = (scores(test_split)>opthresh);  
        [ypred, test_scores] = Mdl.predict(pcaX(test_split, :));
        [fp,tp,~, AUC(rep), ~] = perfcurve(y(test_split),...
            test_scores(:, 2), 1);

        scores(test_split) = test_scores(:, 2);
        acc(rep) = mean(y(test_split) == ypred);

        % calculate weights per initial element 
        % find all weights for components for which there is a 
        % significantly non-zero model weight
        % to account for intercept
        Beta = Mdl.Beta;    
        cutoff = [prctile(Beta, 3) prctile(Beta, 97)];
        % limit the number of weights you are looking at
        if numel(Beta)>20
            mdlcomp = Beta<cutoff(1) | Beta>cutoff(2);
        else
            mdlcomp = ones(n_comp, 1);
        end

        if ~isempty(comp) %&& ~isempty(mdlcomp)
            % revert from PCA space back to electrode x time series
            % columnwise multiplcation of the beta estimates and component
            % weights
            %tmp = comp(:, mdlcomp).*repmat(beta', size(comp, 1), 1);
            %tmp = comp(:, 1:n_comp).*repmat(Beta', size(comp, 1), 1);
            % if isscalar(Beta)                
            %     tmp = comp(:, mdlcomp)*repmat(Beta(mdlcomp)', size(comp, 1), 1);
            % else
                tmp = comp(:, mdlcomp).*repmat(Beta(mdlcomp)', size(comp, 1), 1);
            % end
            % electrode x timepoints x pca components
            tmp = reshape(tmp, [], length(tps), size(tmp, 2));

            % if only one electrode / feature, add a singleton dimension
            % first
            
            % see how weights look across time window (for first component)
            debug = 0;
            if debug && rep 
                
                figure; 
                imagesc(squeeze(tmp(:, :, 1)));

                % same thing using fitglm insstead
                mdl = fitglm(pcaX(train_split, :), y(train_split), ...
                    'Distribution', 'binomial', 'Link','logit');
                %train_scores = mdl.Fitted.Probability;
                mdlcomp = find(mdl.Coefficients.pValue(2:end)<0.0001);
                beta = mdl.Coefficients.Estimate(mdlcomp+1);

                disp(['Correlation between fitglm and fitclinear: ' ...
                    num2str(corr(mdl.Coefficients.Estimate(2:end), Mdl.Beta))])

                figure; 
                scatter(mdl.Coefficients.Estimate(2:end), Beta); hold on;
                xlabel('mdl Coefficient Estimate');
                ylabel('fitclinear estimate')
                scatter(beta, Beta(mdlcomp), 25, 'filled');
            end
            weights(rep, :, :) = squeeze(sum(tmp, 3))';
        else
            weights(rep, :) = Beta;
        end
    end       
end

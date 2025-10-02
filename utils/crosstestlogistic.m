function [AUC, pcaX, ypred] = crosstestlogistic(X, y, Mdl, comp, mu, nreps)
    % X should be 2d matrix: trial x (electrode x time)
    % run logistic regression after PCA on design matrix over nfolds
    % leave comp param empty if pcaflag on
     
    % Perform the PCA construction from the trained model input
    pcaX = (X - mu) * comp;

    rng(1);
    % COmpute the AUC over fixed cross-validation test sets
    AUC = nan(nreps, 1);
    cv = cvpartition(size(pcaX,1),'KFold',nreps);
    ypred = nan(size(pcaX,1), 1);
    for n = 1:nreps
        test_split = cv.test(n);
        [label, test_scores] = Mdl.predict(pcaX(test_split, :));
        [~, ~, ~, AUC(n), ~] = perfcurve(y(test_split), test_scores(:, 2), 1);  
        ypred(test_split) = label == y(test_split)';

        clear test_scores
    end
end
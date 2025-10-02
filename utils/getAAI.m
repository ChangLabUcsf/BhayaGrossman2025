function [ambiguity] = getAAI(window, Swrd, fieldname, pcaflag, nfolds)

    if nargin<3, fieldname = 'aud'; end
    if nargin<4, pcaflag = 1; end
    if nargin<5, nfolds = 10; end
    
    % perform classification on the acoustic feature specified by fieldname
    tps = 51-round(window/2):51+round(window/2); % around
    
    % aggregate feature
    X = cat(3, Swrd.(fieldname){:});
    featsz = size(X, 1); 

    % vectorize feature 
    X = reshape(X(:, tps, :), featsz*length(tps), []);

    % run classification
    y = Swrd.wordOns;
    [~, ~, ~, ~, scores, ~, ~] = logistic(X', y, pcaflag, [], tps, nfolds);  

    % define acoustic ambiguity index
    % this is the difference between the score and correct label (0 or 1)
    ambiguity = abs(scores-y);
end

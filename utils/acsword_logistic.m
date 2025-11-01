function [fp, tp, AUC, pcaX, scores, acc, weights, Mdl, comp, mu] = ...
acsword_logistic(Swrd, field, trlidx, nreps, timing, pcaflag, corpus_details)

    if nargin<6
        pcaflag = 0;
    end

    % phn is a special case (use prec and ons phonemes)
    if ~contains(field, 'phn')
       data = [Swrd.(field)];
    end

    y = Swrd.wordOns(trlidx)>0;
    phnfeatmat = corpus_details.features.mat;

    disp(['timing: ' num2str(timing)])
    tps = 51-round(timing/2):51+round(timing/2); % around

    switch field
        case 'env'
            X_tmp = cell2mat(data(trlidx));
            X = X_tmp(:, tps);
        case 'formant'
            X_tmp = cat(3, data{trlidx});
            X = reshape(X_tmp(1:2, tps, :), 2*length(tps), [])';
        case 'aud'
            X_tmp = cat(3, data{trlidx});
            X = reshape(X_tmp(:, tps, :), 80*length(tps), [])';
        case 'phn' 
            X = [Swrd.precPhn(trlidx) Swrd.onsPhn(trlidx)]; 

            % only use trials with both preceding and succeeding phonemes
            phnidx = X(:, 1)>0 & X(:, 2)>0 & ~any(isnan(X), 2);

            X = X(phnidx, :);
            y = y(phnidx);
        case 'phnfeat'
            X_tmp = [Swrd.precPhn(trlidx) Swrd.onsPhn(trlidx)]; 

            % only use trials with both preceding and succeeding phonemes
            phnidx = X_tmp(:, 1)>0 & X_tmp(:, 2)>0 & ~any(isnan(X_tmp), 2);

            % translate to phonetic features
            onset = arrayfun(@(x) phnfeatmat(:, x), ...
                X_tmp(phnidx, 1), 'UniformOutput', false);
            offset = arrayfun(@(x) phnfeatmat(:, x), ...
                X_tmp(phnidx, 2), 'UniformOutput', false);

            X = cell2mat([onset offset]')';
            y = y(phnidx);
        case 'consphnfeat'
            X_tmp = [Swrd.precPhn(trlidx) Swrd.onsPhn(trlidx)]; 

            % only use trials with both preceding and succeeding phonemes
            phnidx = X_tmp(:, 1)>0 & X_tmp(:, 2)>0 & ~any(isnan(X_tmp), 2);

            % translate to phonetic features
            consfeat = ~ismember(corpus_details.features.names, ...
                {'high', 'front', 'low', 'back', 'syllabic'});
            onset = arrayfun(@(x) phnfeatmat(consfeat, x), ...
                X_tmp(phnidx, 1), 'UniformOutput', false);
            offset = arrayfun(@(x) phnfeatmat(consfeat, x), ...
                X_tmp(phnidx, 2), 'UniformOutput', false);

            X = cell2mat([onset offset]')';
            y = y(phnidx);
    end

    clear X_tmp
    
    if pcaflag
        [fp, tp, AUC, pcaX, scores, acc, weights, ~, Mdl, comp, mu] = ...
            logistic(X, y, pcaflag, [], tps, nreps);
    else
        [fp, tp, AUC, pcaX, scores, acc, weights, ~, Mdl] = ...
            logistic(X, y, pcaflag, [], tps, nreps);
        comp = [];
        mu = [];
    end

end
    

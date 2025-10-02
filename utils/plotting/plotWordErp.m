%% plot erps over word, syllable, phoneme onsets
function plotWordErp(Dvow, SID, singleSet, ~, f, cols, type, bef, baseline, plotmean)
    addpath brewer
    desmat.names = [];
    if nargin<9
        baseline = 0;
    end

    if nargin<10
        plotmean = 1;
    end

    resp = Dvow.(SID).resp;
    switch type
        case 0 % just word onset ERP
            desmat.names = {'word'};
            desmat.condition =  (Dvow.wordOns);
        case 1 % syllable
            desmat.names = {'syll', 'word'};
            desmat.condition =  (Dvow.syllOns + Dvow.wordOns);
        case 2 % preceding word length
            nanidx = isnan(Dvow.precWordLen) | Dvow.precWordLen>7 ...
                | Dvow.precWordLen<2;
            desmat.condition =  Dvow.precWordLen(~nanidx);
            
            desmat.names = arrayfun(@num2str, unique(desmat.condition),...
                'UniformOutput', 0);
            resp = resp(:, :, ~nanidx);
    end   

    if nargin<6 
        cols = brewermap(5, 'Dark2'); 
    end

    if baseline
        dataf = 100;
        % subtract the mean of the baseline period for each electrode separately
        % resp is of size [nElectrodes x nTimepoints x nTrials]
        for i = 1:size(resp, 1)
            resp(i, :, :) = resp(i, :, :) - mean(resp(i, 1:(bef*dataf), :), 2);
        end
    end

    if plotmean
        % plot the mean ERP
        plotResp(resp, desmat, singleSet, SID, f, cols, [], bef);
    end
end
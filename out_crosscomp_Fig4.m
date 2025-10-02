%% Set up
% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;
SIDs = [sSIDs, eSIDs, {'HS11', 'HS9', 'HS10'}];

% change this for the word onset analysis )
bef=50;
aft=50;

if ~exist('Dwrd', 'var')
    load("data/Figure3/Figure3_DIMEXWrd.mat");
end

if ~exist('TDwrd', 'var')
    load("data/Figure3/Figure3_TIMITWrd.mat");
end

window = 20;
fieldname = 'aud';

Dwrd.ambiguity = getAAI(window, Dwrd, fieldname, 20);
TDwrd.ambiguity = getAAI(window, TDwrd, fieldname, 20);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;
%% A - Neural word-boundary decoding 

figure('Position', [100, 300, 250, 400], 'Renderer', 'painters');

% Specifies the type of decoding
decode_type = 'word'; 
ctr = 1;
accls = cell(2, 1);
for s = {Dwrd, TDwrd}
    Swrd = s{1};
    time_label = '600ms';
    % time_label = '10ms-sliding';
    lss = [1, 2, 4];
    type = '';
    
    switch time_label
        case '10ms-sliding'
            accls{ctr} = nan(4, 13, 20); % Initializes a matrix for storing AUC values
        otherwise
            accls{ctr} = nan(4, 20);
    end
     
    for ls = lss  
        if ls == 3
            subj = 'mandarinmono';
        elseif ls == 4
            subj = 'bilingual';
        else
            subj = 'monolingual';
        end
    
        if startsWith(Swrd.sentid{1}, 's')
            corpus = 'DIMEX';
        else
            corpus = 'TIMIT';
        end
        
        filename =  [corpus '_' decode_type '_decode_' subj '_' time_label type '.mat']; 
        load([datapath 'Figure3/decode/' filename], 'decode_details');
    
        switch time_label
            case '10ms-sliding'
                % Stores AUC values in the matrix
                accls{ctr}(ls, :, :) = squeeze(decode_details.AUC(ls, :, :)); 
            otherwise
                accls{ctr}(ls, :) = squeeze(decode_details.AUC(ls, :, :));
        end 

        % Print out information about the decoding (number of elecs, unique subjects, trials)
        disp([subj ':' num2str(ls) ' ' corpus ' ' time_label ' ' type ' ' decode_type]);
        elecs = decode_details.elecs{ls};
        pcaX = decode_details.pcaX{ls};
        disp(['Number of elecs: ' num2str(length(elecs))]);
        disp(['Unique subjects: ' num2str(length(unique(decode_details.encoding.SID(elecs))))]);
        disp(['Number of components retained: ' num2str(size(pcaX, 2))])
        if isfield(decode_details, 'ys')
            disp(['Syllable trials: ' num2str(sum(decode_details.ys{ls}==0)) ' Word trials: ' ...
                num2str(sum(decode_details.ys{ls}==1))]);
        else
            disp(['Syllable trials: ' num2str(sum(decode_details.y==0)) ' Word trials: ' ...
                num2str(sum(decode_details.y==1))]);
        end
        disp(['AUC: ' num2str(median(decode_details.AUC(ls, :)))]);
        disp('-------------------------------------------------------------------------------');
    end
    
    % If using sliding window, use the maximum timepoint across subject groups
    if strcmp(time_label, '10ms-sliding')
        [~, maxtp] = max(mean(accls, [1 3], 'omitnan'));
        accls{ctr} = squeeze(mean(accls{ctr}(:, maxtp-4: maxtp+4, :), 2));
    end
    ctr = ctr + 1;
end

% Plot native vs. nonnative
native = [accls{1}(1, :), accls{2}(2, :)];
nonnative = [accls{1}(2, :), accls{2}(1, :)];
plts{1} = [native; nonnative];
bil = [accls{1}(4, :), accls{2}(4, :)];
plts{2} = [native; bil];
labels = {{'native', 'foreign'}, {'monolingual', 'bilingual'}};

cols = {'b', 'r'};
for t = 1 % Type
    subplot(1, 1, t);
    plt = plts{t};
    for c = 1:2
        h = boxchart(ones(size(plt, 2), 1)*c, plt(c, :), ...
            'BoxFaceColor', [1 1 1], 'BoxEdgeColor', 'k'); % Creates boxplots
        h.JitterOutliers = 'on';
        h.MarkerStyle = '.';
        h.MarkerColor = 'k';
        hold on;
        scatter(randn(size(plt, 2), 1)*0.1+c-0.05, plt(c, :), 25, 'filled', ...
            cols{c});
        disp(['Median AUC for boxplot: ' num2str(median(plt(c, :))) ...
            ', SD=' num2str(std(plt(c, :)))])

        hold on;
    end

    % Perform a ttest on the native vs. non-native decoding
    [h, p] = ttest2(plt(1, :), plt(2, :), "Tail", "right"); 
    line([1.25 1.75], [.85, .85], 'Color', 'k', 'LineWidth', 1.5); % Draws a line
    text(1.35, .87, getSigStr(p, 2), 'FontSize', 13); % Adds text to the plot
    
    % Formatting
    xlabel('Group Type');
    ylabel('AUC');

    ylim([0.45 0.9]);
    yticks(0.5:0.2:0.9)
    yline(0.5, 'Color', 'k');

    xlim([0.5 2.5])
    xticks([1 2]);  
    xticklabels(labels{t})
    
    set(gca, 'FontSize', 15);
    box off;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% B - Histogram of AAI and example low / high AAI instances

titles = {'English', 'Spanish'};
figure;
subplot(2, 1, 2);
histogram(TDwrd.ambiguity, 'EdgeColor', 'none', 'FaceColor', ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'Normalization', 'probability'); hold on;
ylabel('Probability');
yticks([0.0 0.2 0.4]);
xticks(0:0.5:1)
xlabel('Acoustic Ambiguity Index (AAI)');
xline([0.949, 0.0322]) % AAI of examples
yl = ylim();
text([0.949, 0.0322], [yl(2), yl(2)], {'high AAI example', 'low AAI example'}, ...
    'FontSize', 13, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', 13);
title(titles{1}, 'FontSize', 15, 'FontWeight', 'normal'); 
box off;

subplot(2, 1, 1);
histogram(Dwrd.ambiguity, 'EdgeColor', 'none', 'FaceColor', ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'Normalization', 'probability'); hold on;
ylabel('Probability');
yticks([0.0 0.2 0.4]);
xticks([]);
set(gca, 'FontSize', 13);
xline([0.9883, 0.01214]) % AAI of examples
% 'high AAI example' 0.9883 and 'low AAI example' 0.01214
yl = ylim();
text([0.9883, 0.01214], [yl(2), yl(2)], {'high AAI example', 'low AAI example'}, ...
    'FontSize', 13, 'HorizontalAlignment', 'center');
title(titles{2}, 'FontSize', 15, 'FontWeight', 'normal'); 
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% B - High and low AAI examples in Spanish and English 

% highambig = 'msdh0_si2240'; % try another
% lowambig = 'mrew1_si2130'; % this bother
% % wrdidx = find(strcmp(Dwrd.sentid, highambig));
wrdidx = [1832, 1676];
maxnum = 7;
Swrd = TDwrd;

% % find 5 highest and 5 lowest AAI that are also word onsets
[~, idx] = sort(Swrd.ambiguity, 'descend');
wrdidx = idx(Swrd.wordOns(idx)>0);
wrdidx = wrdidx(1:maxnum);

% find 5 lowest AAI that are also word onsets
[~, idx] = sort(Swrd.ambiguity, 'ascend');
wrdidx = [wrdidx; idx(Swrd.wordOns(idx)>0)];
wrdidx = wrdidx(1:maxnum*2);

xdata = -0.5:0.01:0.5;
% show spectrogram, and envelope and text
for wrd = wrdidx'
    figure;
    trial = Swrd(wrd, :);
    aud = trial.aud{1};
    env = trial.env{1};
    % combine precPhns and onsPhns with a | in the middle
    % text = [join(num2phns(trial.precPhns{1}, timit_details), '') '|' ...
    %     join(num2phns(trial.onsPhns{1}, timit_details), '')];
    % remove all 'cl'
    text = join(trial.phns{1}, '');
    % plot spectrogram on the left and envelope on the right
    subplot(1, 1, 1);
    imagesc(xdata, 1:80, aud);
    colormap(flipud(gray));
    brighten(0.5);
    xline(0, 'LineWidth', 2, 'Color', 'r');
    xlabel('Time (s)');
    xticks([-0.5 0 0.5]);
    yticks([1 80]);
    yticklabels({'0', '8'});
    ylabel('kHz');
    
    set(gca, 'YDir', 'normal', 'FontSize', 13);

    yyaxis right;
    plot(xdata, env, 'Color', 'k', 'LineWidth', 3);
    yticks([]);

    xlim([-0.2 0.2]);
    xticks([-0.2 0 0.2]);
    sgtitle([trial.precword{1} ' | ' trial.currword{1} ' ' ...
        'AAI=' num2str(trial.ambiguity(1))]);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% C - Average Envelope


befaft = [0.5 0.5];
x = -befaft(1):0.01:befaft(2);

colors = getColorsCrossComp(1);
colors = [0.2 0.2 0.2; 0.2 0.2 0.2];
addpath(genpath('shadederror/'))

figure;
ctr=1;
titles = {'Spanish', 'English'};
for i = {Dwrd, TDwrd}
    subplot(2, 1, ctr)
    ywrd = squeeze(cat(3, i{1}.env{logical(i{1}.wordOns)}))';
    ysyll = squeeze(cat(3, i{1}.env{~i{1}.wordOns}))';
    
    % Plot the average envelope for word and syllable
    shadedErrorBar(x, mean(ysyll, 1), nansem(ysyll, 1), ...
        'lineprops', {'color', [0.8 0.8 0.8], 'linewidth', 1.8});
    shadedErrorBar(x, mean(ywrd, 1), nansem(ywrd, 1), ...
        'lineprops', {'color', colors(ctr, :), 'linewidth', 2.5}); hold on;

    % Calculating the F-stat values over time
    [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
        reshape([ysyll; ywrd]', 1, size([ysyll; ywrd], 2), size([ysyll; ywrd], 1)), ...
        [ones(1, size(ysyll, 1)) 2*ones(1, size(ywrd, 1))], [1, 2]);
    fthresh = finv(1-0.0001, df1, df2);  

    % Visualizing the F-stat values on average envelope plot
    % scatter(x(fvals>fthresh), 0.1*ones(1, sum(fvals>fthresh)), 15, ...
    %     fvals(fvals>fthresh), 'filled');
    % cm = colormap("gray");
    % colormap(flipud(cm(1:200, :)));

    % Reference lines and formatting
    xline(0, 'LineWidth', 2, 'LineStyle', '-', 'Color','r');
    ylabel('Envelope(norm)'); 
    if ctr==2
        xlabel('Time (s)');
    end
    legend({'syllable', 'word'});
    set(gca, 'FontSize', 13);
    %ylim([-0.01 1]);
    yticks(0:0.5:1)
    title(titles{ctr}, 'FontWeight', 'normal');
    xlim([-0.2 0.2]);
    xticks([-0.2 0.2]);

    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;


%% D - Neural word boundary decoding split by AAI 

% fieldname = 'aud';
% window = 20;
% Dwrd.ambiguity = getAAI(window, Dwrd, fieldname, 20);
% TDwrd.ambiguity = getAAI(window, TDwrd, fieldname, 20);

% red for unfamiliar language, blue for native
colors = [0.8 0.2 0.2; 0.2 0.2 0.8];
corpi = {'DIMEX', 'TIMIT'};
titles = {'Spanish', 'English'};
Swrds = {Dwrd, TDwrd};

numbins = 3;
folds = 15;
figure; 
AUC = nan(2, numbins, folds);

for cr = 1:2
    corpus = corpi{cr};
    Swrd = Swrds{strcmpi(corpus, 'TIMIT')+1};
    time_label = '600ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'Figure3/decode/' filename], 'decode_details');

    % analysis is sensitive to these binedges!
    % bins constructed based on prctiles
    binedges = prctile(Swrd.ambiguity, 0:100/numbins:100);
    
    for ls = [1, 2]
        % get all the neural data for this bin
        pcaX = decode_details.pcaX{ls};
        yidx = decode_details.yidx{ls};
        y = decode_details.ys{ls};
        
        ambig = discretize(Swrd.ambiguity(yidx), binedges);
        Mdl =  decode_details.Mdl{ls};
        for bin = 1:numbins
            idx = find(ambig==bin);

            % shuffle the indices
            rng(1);
            idx = idx(randperm(length(idx)));
            for f = 1:folds
                % choose numind set of indices from the shuffle list
                numind = floor(length(idx)/folds);
                fidx = idx((f-1)*numind+1:f*numind);
                
                % get the pcaX for this bin
                pcaX_bin = pcaX(fidx, :);
                y_bin = y(fidx);
                
                % run the trained Mdl to get the AUC for the bin
                [~, scores] = predict(Mdl, pcaX_bin);
                [~, ~, ~, AUC(ls, bin, f)] = perfcurve(y_bin, scores(:, 2), 1);
            end
        end
    end
    
    AUCdiff = nan(folds*folds, numbins);
    for i = 1:bin
        if strcmpi(corpus, 'DIMEX')
            tmp = squeeze(AUC(1, i, :))-squeeze(AUC(2, i, :))';
        else % native - foreign
            tmp = squeeze(AUC(2, i, :))-squeeze(AUC(1, i, :))';
        end
        AUCdiff(:, i) = tmp(:);
    end

    subplot(1, 2, strcmpi(corpus, 'TIMIT')+1);
    b = boxchart(AUCdiff, 'Notch', 'on', 'MarkerStyle', ...
        'none', 'BoxWidth', 0.4, 'BoxFaceColor', [0.5 0.5 0.5], 'BoxEdgeColor', 'k'); hold on;
    b.WhiskerLineColor = [1 1 1];
    x = [ones(size(AUCdiff, 1), 1); 2*ones(size(AUCdiff, 1), 1); ...
        3*ones(size(AUCdiff, 1), 1)];
    y = [AUCdiff(:, 1); AUCdiff(:, 2); AUCdiff(:, 3)];

    % violin plot
    % violinplot(AUCdiff, {'low', 'med', 'high'}, 'MarkerSize', 5, ...
    %     'ShowWhiskers', false, 'EdgeColor', [1 1 1], 'BoxColor', [0 0 0], ...
    %     'ViolinColor', [0.5 0.5 0.5]);

    % swarmchart option
    % swarmchart(x, y, 10, [0.6 0.6 0.6], 'filled');
    % for j = 1:3
    %     line([j-0.2 j+0.2], [y(j) y(j)], 'LineWidth', 2, 'Color', 'k');
    % end

    % barchart not aligned with nature formatting
    % bar(y, 'FaceColor',  [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on;
    % errorbar(y, nansem(AUCdiff), '.', 'CapSize', 0, 'Color', 'k', 'LineWidth', 2);
    %
    % for j = 1:3
    %     x = j-0.1+0.2*rand(size(AUCdiff(:, j)));
    %     scatter(x, AUCdiff, 5, [0.5 0.5 0.5], 'filled');
    % end

    % linear model
    mdl = fitlm(y, x);
    % disp(mdl);

    % ANOVA for a categorical model since they are discrete bins
    aov = anova(AUCdiff); % 1-way anova
    disp(aov)
    disp(['P-value: ' num2str(aov.stats.pValue(1)) ', F: ' num2str(aov.stats.F(1))]);

    % plot the line trend
    % x = 1:numbins;
    % y = median(AUCdiff, 1);
    % plot(x, y, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'Marker', 'o'); hold on;
    
    ylim([-0.05 0.13]);
    yticks(0:0.05:0.5);
    yline(0);
    %xticks([1 numbins]);
    ylabel('Native - Foreign AUC');
    xticklabels({'low', '', 'high'});
    xlabel('Acoustic ambiguity (AAI)');
    set(gca, 'FontSize', 15);
    box off;
    title(titles{cr}, 'FontWeight', 'normal');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

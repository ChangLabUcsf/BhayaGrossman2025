%% Set up
out_crosscomp_startup;

% timit vowels are all categories with > 100 trials in stressed type
timit_vow = {'aa', 'ae', 'ao', 'ah', 'ay', 'ey', 'eh', 'ih', 'iy', 'ow'};
dimex_vow = {'a', 'e', 'i', 'o', 'u'};

if ~exist('Dvow', 'var')
    load("data/ExtendedFigures/ExtFigure6_DIMEXVow.mat");
end

if ~exist('TDvow', 'var')
    load("data/ExtendedFigures/ExtFigure6_TIMITVow.mat");
end

if ~exist('imgall', 'var')
    load([datapath '/Figure1/Figure1_ImgData.mat'], 'imgall');
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% A - Shaded error bar: Vowel ERP

SID = 'EC183'; 
el = 121; 

SID = 'EC222'; 
el = 101; % 118 

plotVowelErp(Dvow,SID, el, 'dimex', figure); 
plotVowelErp(TDvow,SID, el, 'timit', figure);

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% UNUSED: Vowel category counts

min_trials = 100;
% load in all vowel categories
Svows = {Dvow, TDvow};
fig = figure;
for vow = 1:2
    Svow = Svows{vow};
    corpus = Svow.corpus;
    idx = Svow.stress==1;

    labels = unique(Svow.vowel(idx));
    % get the counts for each vowel category using the cell array
    counts = cellfun(@(x) sum(strcmp(Svow.vowel(idx), x)), labels);

    % plot count for each vowel category
    subplot(1, 2, vow);
    bar(counts);
    title(corpus);
    xticks(1:length(labels));
    xticklabels(labels);
    ylabel('Count');
    set(gca, 'FontSize', 15);
    box off;
    yline(min_trials, 'Color', 'k', 'LineWidth', 2); % minimum number of trials for decent LDA

    % print out which vowels have more than 100 trials
    disp(corpus);
    disp(labels(counts > min_trials));
    disp('--------------------------------');
end
sgtitle('Vowel Counts');

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;


%% B - FIX COLORS Vowel formant plot

figure;
timit_vow = {'aa', 'ae', 'ao', 'ah', 'ay', 'ey', 'eh', 'ih', 'iy', 'ow'};
dimex_vow = {'a', 'e', 'i', 'o', 'u'};
incl_vows = {dimex_vow, timit_vow};
Svows = {Dvow, TDvow};
cols = {brewermap(length(dimex_vow), 'Dark2'), brewermap(length(timit_vow), 'Dark2')};

for s = 1:2
    subplot(1, 2, s)

    Svow = Svows{s};
    incl_vow = incl_vows{s};
    for j=1:length(incl_vow)       
        isVow = cellfun(@(x)isequal(x,incl_vow{j}),Svow.vowel);
        vowIdx = find(isVow & Svow.stress==1);
       
        if ~isempty(vowIdx)
            scatter_col=repmat(cols{s}(j, :), length(vowIdx), 1);
        end
        scatter(Svow.formantVals(2, vowIdx), Svow.formantVals(1, vowIdx), 25, scatter_col, 'filled', ...
        'HandleVisibility', 'off', 'MarkerFaceAlpha', 0.5); hold on;           
    end
    xlabel('F2 (Hz)'); 
    ylabel('F1 (Hz)');   
    set(gca, 'YDir', 'reverse', 'XDir', 'reverse');
    ylim([200 1000]);
    xlim([650 3100]);
    yticks([]);
    xticks([]);
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% C/D - Acoustic decoding of vowel categories (LDA)

% Note: Unsupervised dimensionality reduction is not used here because MDS scale is not
% able to handle the high dimensionality of the audio data 
% was set at 35
max_comp = 200;

feat = 'formants';
if strcmp(feat, 'aud')
    pcaflag = 1;
else
    pcaflag = 0;
end
kfolds = 10;

tps = 50:60;
% make sure the ordering is always the same for each variable
Svows = {Dvow, TDvow};

incl_vows = {dimex_vow, timit_vow};
corpora = {'dimex', 'timit'};

acc_acs = nan(2, 1);
acc_fold = nan(2, kfolds);
auc_fold = nan(2, kfolds);
ypredlabel = cell(2, 1);
for ctr = 1:2
    Svow = Svows{ctr};
    corpus = corpora{ctr};
    idx = Svow.stress==1 & ismember(Svow.vowel, incl_vows{ctr});

    if strcmp(feat, 'aud')
        allaud = Svow.aud(:, :, idx);
        % preprocessing and pca
        A = reshape(allaud(:, tps, :), [], size(allaud, 3));
    else
        allaud = Svow.formantVals(:, idx);
        A = allaud;
    end

    y =  Svow.vowelType(idx)';
    if pcaflag
        % with pca
        [coeff, score, ~, ~, exp] = pca(A', 'Algorithm', 'eig', ...
            'NumComponents', max_comp);
        n_comp = find(diff(cumsum(exp)>95));
    else
        % no pca
        coeff = A';
        n_comp = size(A, 1);
        score = A';
    end
    
    % run lda with kfold cross-validation
    y = arrayfun(@(x) find(x==unique(y)), y);
    if strcmp(feat, 'aud')
        [y_pred, mappedX, ~, ~] = LDAmap(score(:, 1:n_comp), y, 5);
        mappedX = squeeze(mean(mappedX(:, :, :), 'omitnan'));
    else
        [y_pred, ~, ~, ~] = LDAmap(score(:, 1:n_comp), y, 5);
        mappedX = A(1:3, :)';
    end
    [~, ~, auc, ~, ~, acc_tmp] = lda(score(:, 1:n_comp), y', ...
        1, [], [], kfolds);
    
    Svow.corpus = corpus;
    conf_acs = pdist2(mappedX(:, 1), mappedX(:, 2));
    ypredlabel{ctr} = [y, y_pred];
    
    % calculate silhouette score
    s = silhouette(mappedX,y);
    ri = RandIndex(y, y_pred);
    
    disp(['---------------------- ' upper(corpus) ' LDA -------------------------------']);
    % print number of components, print number of trials per category
    disp(['n components = ' num2str(n_comp)]);
    disp(['n trials = ' num2str(length(y))]);
    % minimum number of trials per category
    counts = cellfun(@(x) sum(strcmp(Svow.vowel(idx), x)), ...
        unique(Svow.vowel(idx)));
    disp(['min trials per category = ' num2str(min(counts))]);
    disp(['RandIndex = ' num2str(ri)]);
    disp(['mSI = ' num2str(mean(s))]);
    disp(['accuracy = ' num2str(sum(y_pred==y)/length(y))]);
    acc_acs(ctr) = sum(y_pred==y)/length(y);
    acc_fold(ctr, :) = acc_tmp;
    auc_fold(ctr, :) = auc;
    clear auc acc_tmp 
end

perfmetric = acc_fold;

figure;
subplot(1, 3, 1);
% scatter(1:2, acc_acs.*100, 105, 'k', 'LineWidth', 3); hold on;
% Create a scatter plot for the accuracy data
randjitter = @(x) x + randn(size(x))*0.05;
scatter(randjitter(ones(size(perfmetric, 2), 1)), perfmetric(1, :) * 100, ...
    45, 'filled', 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k'); hold on;
scatter(randjitter(ones(size(perfmetric, 2), 1) * 2), perfmetric(2, :) * 100, ...
    45, 'MarkerEdgeColor', 'k');
% Overlay a black box plot
boxplot(perfmetric' .* 100, 'Colors', 'k', 'Labels', {'Spanish', 'English'}, 'Widths', 0.5);
text(1, 80, num2str(median(perfmetric(1, :))), 'FontSize', 14);
text(2, 50, num2str(median(perfmetric(2, :))), 'FontSize', 14);

plot([0.75 1.25], [1/length(dimex_vow) 1/length(dimex_vow)].*100, 'Color', 'k', 'LineWidth', 2);
plot([1.75 2.25], [1/length(timit_vow) 1/length(timit_vow)].*100, ...
    'Color', 'k', 'LineWidth', 2);
text(1, (1/length(dimex_vow))*100, 'chance %', 'FontSize', 14);
text(2, (1/length(timit_vow))*100, 'chance %', 'FontSize', 14);
ylim([0 100]);
set(gca, 'FontSize', 15);
ylabel('Classifier Accuracy (%)');
yticks(0:50:100);
xticks(1:2);
xticklabels({'Spanish vowels', 'English vowels'});
box off;

for vow = 1:2
    subplot(1, 3, vow+1);
    vows = incl_vows{vow};
    ytrue = arrayfun(@(x) vows{x}, ypredlabel{vow}(:, 1), 'UniformOutput', false);
    ypred = arrayfun(@(x) vows{x}, ypredlabel{vow}(:, 2), 'UniformOutput', false);
    [C, order] = confusionmat(ytrue, ypred, 'order', vows);
    % make this be percentage
    C = C./ sum(C, 2)*100;   
    imagesc(C);
    colormap(flipud(magma));
    brighten(0.5);
    set(gca, 'FontSize', 15, 'XDir', 'normal', 'YDir', 'normal');
    cl = clim();
    
    yticks(1:length(vows));
    yticklabels(vows);
    xticks(1:length(vows));
    xticklabels(vows);
    title(corpus);
    xlabel('Predicted');
    ylabel('True');

    box off;
    cbh = colorbar;
    if vow==1
        clim([0 350]);
        cbh.Ticks = [0 300];
    else
        clim([0 150]);
        cbh.Ticks = [0 100];
    end
    clim([0 100])
    cbh.Ticks = [0 100];
    cbh.FontSize = 15;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;
%% E - Visualize single subject results from both corpora together

figure;
subj = 'monolingual';
timelabel = '100-400ms'; % after onset
corpora = {'dimex', 'timit'};
titles = {'Spanish vowels', 'English vowels'};
confmat = cell(2, 1);
field = 'AUC';
incl_vows = {dimex_vow, timit_vow};
vow_counts = cell(2, 1);
for ctr = 1:2
    corpus = corpora{ctr};
    subplot(1, 2, ctr);

    % load in data
    load([datapath 'ExtendedFigures/decode/' corpus ...
        '_vowel_decode_' subj '_' timelabel '.mat']);

    % plot the median AUC for each subject, colored by language
    confmat{ctr} = nan(height(vowel_decode_table), length(incl_vows{ctr}), ...
        length(incl_vows{ctr}));
    
    for s = 1:height(vowel_decode_table)
        numel = length(vowel_decode_table.Electrodes{s});
        randjitter = @(x) x + randn(size(x))*0.05;
        native = (vowel_decode_table.Language(s) == 1 && ctr==1) || ...
                (vowel_decode_table.Language(s) == 2 && ctr==2);
        if native
            scatter(1+randjitter(0), median(vowel_decode_table.(field){s}), numel*2, 'filled', ...
                'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5); hold on;
        else
            scatter(2+randjitter(0), median(vowel_decode_table.(field){s}), numel*2, 'filled', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5); hold on;
        end
        vowel_decode_table.medianAccuracy(s) = median(vowel_decode_table.(field){s});
        vowel_decode_table.native(s) = ~native;
        vowel_decode_table.hemi(s) = {imgall.(vowel_decode_table.SID{s}).hemi};

        tbl = tabulate(vowel_decode_table.Labels{s});
        [~, y] = ismember(incl_vows{ctr}, [tbl(:, 1)]);
        tbl = tbl(y, :);
        vow_counts{ctr}(s, :) = cell2mat(tbl(:, 2));

        % get confusion matrix
        y = vowel_decode_table.Labels{s};
        ypred = vowel_decode_table.Predicted{s};
        confmat_tmp = confusionmat(y, ypred, 'order', incl_vows{ctr}); 
        confmat{ctr}(s, :, :) = (confmat_tmp./sum(confmat_tmp, 2));
    end

    % add boxplot
    boxplot(vowel_decode_table.medianAccuracy, vowel_decode_table.native, ...
        'Colors', 'k', 'Symbol', '');

    % look at difference between languages (lme)
    lme_tbl = table();
    for s = 1:height(vowel_decode_table)
        nreps = length(vowel_decode_table.(field){s});
        tmp_tbl = table();
        tmp_tbl.auc = vowel_decode_table.(field){s};
        tmp_tbl.native = repmat(vowel_decode_table.native(s), nreps, 1);
        tmp_tbl.SID = repmat(vowel_decode_table.SID(s), nreps, 1);
        tmp_tbl.elecs = repmat(length(vowel_decode_table.Electrodes{s}), nreps, 1);
        tmp_tbl.hemi = repmat(length(vowel_decode_table.hemi{s}), nreps, 1);
        lme_tbl = [lme_tbl; tmp_tbl];
    end
    lme = fitlme(lme_tbl,'auc~1+native+(1|elecs)+(1|hemi)');
    text(1.5, 0.8, getSigStr(lme.Coefficients.pValue(2), 1), 'FontSize', 20);

    xticks(1:2);
    xticklabels({'Native', 'Foreign'});
    xlabel('Language known');
    
    if ctr==1
        ylabel('Classifier AUC');
        yticks([0.5 0.7 0.9]);
    else
        yticks([]);
    end

    ylim([0.4 0.9]);
    chance = 0.5;

    set(gca, 'FontSize', 15);
    box off;
    xlim([0.5 2.5]);
    yline(chance, 'LineStyle', '--', 'LineWidth', 2);
    title(titles{ctr}, 'Fontweight', 'normal');
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% F - Visualize single electrode vowel category decoding

corpora = {'dimex', 'timit'};
subj = 'monolingual';
timelabel = '100-400ms'; % after onset
for corp = 1:2
    corpus = corpora{corp};
    filename = [corpus '_vowel_decode_' subj '-elec_' timelabel '.mat']; % dimex filename
    load([datapath 'ExtendedFigures/decode/' filename], 'vowel_decode_table');

    p = nan(height(vowel_decode_table), 1); 
    for i = 1:height(vowel_decode_table)
        [~, p(i)] = ttest(vowel_decode_table.AUC{i}, 0.5, 'Tail', 'right');
    end
    vowel_decode_table.pAUC = p;

    subplot(1, 2, corp);
    medianauc = cellfun(@(x) median(x), vowel_decode_table.AUC);
    ls = [vowel_decode_table.ls{:}]==corp; % native
    violinplot({medianauc(ls), ...
        medianauc(~ls)}, {'native', 'foreign'}, ...
        'ViolinColor', {[0 0 1] ; [1 0 0]},  'MarkerSize', 6);
    
    xticks(1.5);
    xticklabels({corpus});
    ylim([0.43 0.7]);
    yticks(0.5:0.1:0.7);
    yline(0.5);
    if corp==1
        ylabel('Classifier AUC');
    end

    set(gca, 'FontSize', 13);
    box off;
    disp(upper(corpus));
    disp('Corrected---------')
    corrected_p = p'<(0.05/height(vowel_decode_table)); 
    disp(['Num sig in native: ' num2str(sum(ls==0&corrected_p)/sum(ls==0))])
    disp(['Num sig in foreign: ' num2str(sum(ls==1&corrected_p)/sum(ls==1))])
    disp('Uncorrected---------')
    uncorrected_p = p'<0.01; 
    disp(['Num sig in native: ' num2str(sum(ls==0&uncorrected_p)/sum(ls==0))])
    disp(['Num sig in foreign: ' num2str(sum(ls==1&uncorrected_p)/sum(ls==1))])

    % do a single tailed t-test for each language
    [p, h] = ranksum(medianauc(ls==0), medianauc(ls==1));
    title(getSigStr(p, 1));
    disp(['Ranksum for ' corpus ': p = ' num2str(p)]);
    disp(['Native mean AUC = ' num2str(mean(medianauc(ls==0)))]);
    disp(['Foreign mean AUC = ' num2str(mean(medianauc(ls==1)))]);
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% G - Single vowel pair analysis ('ih' vs. 'eh')

% plot the scatter of formant information for TIMIT, draw contours of vowels for DIMEx
incl_vows = {{'i', 'e'}, {'ey', 'ih', 'iy',}}; % dimex_vow, timit_vow 'iy', , 'ey'eh
Svow = TDvow;
cols = lines(5);

figure;
for j=1:length(incl_vows{2})       
    isVow = cellfun(@(x) isequal(x,incl_vows{2}{j}),Svow.vowel);
    vowIdx = find(isVow & Svow.stress==1);
    vowCol = find(strcmp(timit_vow, incl_vows{2}{j}));
    % scatter_col=repmat(cols(vowCol, :), length(vowIdx), 1);
    scatter_col=repmat(cols(j, :), length(vowIdx), 1);
    
    scatter(Svow.formantVals(2, vowIdx), Svow.formantVals(1, vowIdx), ...
        25, scatter_col, 'filled', 'HandleVisibility', 'on', ...
        'MarkerFaceAlpha', 0.6); hold on;
    xlabel('F2 (Hz)'); 
    ylabel('F1 (Hz)');
    clear scatter_col;
end
legend(incl_vows{2});

% plot the contour of the vowels for DIMEx on top of the scatter plot
Svow = Dvow;    
cols = brewermap(5, 'Dark2');
cols(5, :) = [0.2 0.6 0.9];
for j=1:length(incl_vows{1})
    isVow = cellfun(@(x)isequal(x,incl_vows{1}{j}),Svow.vowel);
    vowIdx = find(isVow);
    
    % get the histcounts2 of the points for the contour
    [Z, X, Y] = histcounts2(Svow.formantVals(2, vowIdx), Svow.formantVals(1, vowIdx), ...
        'Normalization', 'pdf');
    % make X be the center of the bins
    X = (X(1:end-1) + X(2:end)) / 2;
    % make Y be the center of the bins
    Y = (Y(1:end-1) + Y(2:end)) / 2;
    contour(X, Y, Z', [prctile(Z(:), 95) prctile(Z(:), 95)], ...
        'Color', cols(j, :), 'LineWidth', 2); % , , 'ShowText', 'on'
    hold on;
end
legend([incl_vows{2} incl_vows{1}]);
ylim([200 850]);
yticks([200 800]);
xticks([800 3000]);
xlim([800 3000]);

set(gca, 'YDir', 'reverse', 'Xdir', 'reverse');

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% H - Visualize single subject results from both corpora together

figure;
% imgall = load_allimgdata;
subj = 'monolingual';
timelabel = '100-400ms'; % after onset
corpora = {'timit'};
field = 'AUC';
incl_vows = {{'iy', 'ih'}, {'ih', 'ey'}, {'iy', 'ey'}};

vow_counts = cell(2, 1);
corpus = corpora{1};

ctr = 1;
for vowpair = incl_vows
    vowpairstr = join(vowpair{1}, '-');
    % load in data
    load([datapath 'ExtendedFigures/decode/' corpus ...
        '_vowel_decode_' subj '_' timelabel '_' vowpairstr{1} '.mat']);

    for s = 1:height(vowel_decode_table)
        numel = length(vowel_decode_table.Electrodes{s});
        randjitter = @(x) x + randn(size(x))*0.05;
        native = vowel_decode_table.Language(s) == 2;
        if native
            scatter(ctr+randjitter(0), mean(vowel_decode_table.(field){s}), numel*2, 'filled', ...
                'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5 ); hold on;
        else
            scatter(ctr+1+randjitter(0), mean(vowel_decode_table.(field){s}), numel*2, 'filled', ...
                'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5); hold on;
        end
        vowel_decode_table.medianAccuracy(s) = mean(vowel_decode_table.(field){s});
        vowel_decode_table.native(s) = ~native;
        vowel_decode_table.hemi(s) = {imgall.(vowel_decode_table.SID{s}).hemi};
    end

    % add boxplot
    boxplot(vowel_decode_table.medianAccuracy, vowel_decode_table.native, ...
        'Colors', 'k', 'Positions', [ctr ctr+1]);
    
    % look at difference between languages (lme)
    lme_tbl = table();
    for s = 1:height(vowel_decode_table)
        nreps = length(vowel_decode_table.(field){s});
        tmp_tbl = table();
        tmp_tbl.auc = vowel_decode_table.(field){s};
        tmp_tbl.native = repmat(vowel_decode_table.native(s), nreps, 1);
        tmp_tbl.SID = repmat(vowel_decode_table.SID(s), nreps, 1);
        tmp_tbl.elecs = repmat(length(vowel_decode_table.Electrodes{s}), nreps, 1);
        tmp_tbl.hemi = repmat(length(vowel_decode_table.hemi{s}), nreps, 1);
        lme_tbl = [lme_tbl; tmp_tbl];
    end
    lme = fitlme(lme_tbl,'auc~1+native+(1|elecs)+(1|hemi)');
    text(ctr+0.5, 0.75, getSigStr(lme.Coefficients.pValue(2), 1), 'FontSize', 15);
    ctr = ctr+3;
end

xticks(1.5:3:7.5);

xticklabels(cellfun(@(x) join(x, '-'), incl_vows));
xlabel('Language known');

ylabel('Classifier AUC');
ylim([0.33 0.8]);
yticks([0.5 0.7 0.9]);

set(gca, 'FontSize', 15);
box off;
yline(0.5, 'LineStyle', '--', 'LineWidth', 2); % chance

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% I - Visualize single electrode results from both corpora together

figure;

subj = 'monolingual';
timelabel = '100-400ms'; % after onset
corpora = {'timit'};
field = 'AUC';
incl_vows = {{'iy', 'ih'}, {'ih', 'ey'}, {'iy', 'ey'}};

vow_counts = cell(2, 1);
corpus = corpora{1};

ctr = 1;
for vowpair = incl_vows
    vowpairstr = join(vowpair{1}, '-');
    % load in data
    load([datapath 'ExtendedFigures/decode/' corpus ...
        '_vowel_decode_' subj '-elec_' timelabel '_' vowpairstr{1} '.mat']);

    for s = 1:height(vowel_decode_table)
        numel = length(vowel_decode_table.Electrodes{s});
        
        randjitter = @(x) x + randn(size(x))*0.05;
        vowel_decode_table.medianAccuracy(s) = mean(vowel_decode_table.(field){s});

        native = vowel_decode_table.Language(s) == 2;
        vowel_decode_table.native(s) = ~native;
        vowel_decode_table.hemi(s) = {imgall.(vowel_decode_table.SID{s}).hemi};
    end

    % add boxplot (half-half)
    subplot(1, length(incl_vows), ctr);
    violinplot({vowel_decode_table.medianAccuracy(vowel_decode_table.native), ...
        vowel_decode_table.medianAccuracy(~vowel_decode_table.native)}, ...
        {'native', 'foreign'}, 'ViolinColor', {[0 0 1] ; [1 0 0]},  'MarkerSize', 6);

    % look at difference between languages (ranksum)
    [p, h] = ranksum(vowel_decode_table.medianAccuracy(vowel_decode_table.native), ...
        vowel_decode_table.medianAccuracy(~vowel_decode_table.native));
    disp(['Ranksum for ' corpus ': p = ' num2str(p)]);
    disp(['Spanish mean AUC = ' num2str(mean(vowel_decode_table.medianAccuracy(~vowel_decode_table.native)))]);
    disp(['English mean AUC = ' num2str(mean(vowel_decode_table.medianAccuracy(vowel_decode_table.native)))]);
   
    ctr = ctr+1;

    xticks([1, 2]);
    xticklabels({'Foreign', 'Native'})
    xlabel(vowpairstr);

    ylim([0.3 0.78]);
    yticks([0.5 0.7 0.9]);
    chance = 0.5;
    if ctr==2
        ylabel('Classifier AUC');
    end

    set(gca, 'FontSize', 15);
    box off;
    yline(chance, 'LineStyle', '--', 'LineWidth', 2);
    sgtitle(corpus, 'Fontweight', 'normal');
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps;

%% Neural decoding of vowel categories (LDA), pooled subject

Svow = Dvow;
corpus = 'dimex';
incl_vows = dimex_vow;

Svow = TDvow;
corpus = 'timit';
incl_vows = timit_vow;

SIDs = {sSIDs, eSIDs};
SIDsexcl = {'EC222'};
for j = 1:2
    SIDs{j} = SIDs{j}(~contains(SIDs{j}, SIDsexcl));
end

%Svow = addtoDD(Svow, corpus, bef, aft, [eSIDs sSIDs]);
subj = 'monolingual';
Svow.corpus = corpus;

idx = ismember(Svow.vowel, incl_vows);  %  & Svow.stress==1
load(['select_elec/out_elecs_speechtypeftest_bychan_' corpus '_all.mat']);

% load in electrode responses
bef = 20;
aft = 50;
nfolds = 15;

% use neural window from onset to 300 ms after onset & baseline    
tps = bef+10:bef+40;
timelabel = '100-400ms'; % after onset

Xboth = cell(1, 2);
sidsboth = cell(1, 2);
for lang = 1:length(SIDs)
    % figure out total elecs for all SIDs in this language
    total_elecs = sum(cellfun(@(x) length(allidx.(x)), SIDs{lang}));
    % Initialize variables to store data for the table
    Xboth{lang} = nan(total_elecs, length(tps), length(Svow.vowel(idx)));
    sids = cell(total_elecs, 1);

    ctr = 1;
    for s = SIDs{lang}
        SID = s{1};
        % get speech responsive electrodes
        disp(['loading subject....' SID])
        elecs = allidx.(SID);
        for e = 1:length(elecs)
            Xboth{lang}(ctr, :, :) = squeeze(Svow.(SID).resp(elecs(e), tps, idx)); 
            sids(ctr) = {SID};
            ctr = ctr+1;
        end
    end

    % find electrodes with > 50 NaN trials
    if strcmp(corpus, 'dimex')
        nanelecs = sum(squeeze(mean(isnan(Xboth{lang}), 2)), 2)>7000;
    else
        nanelecs = sum(squeeze(mean(isnan(Xboth{lang}), 2)), 2)>50;
    end
    if sum(nanelecs)>0
        disp(['removing ' num2str(sum(nanelecs)) ' electrodes with > 50 NaN trials']);
        Xboth{lang} = Xboth{lang}(~nanelecs, :, :);
        sids = sids(~nanelecs);
    end
    sidsboth{lang} = sids;
end

% find trials in common across Xboth
trls = cellfun(@(x) squeeze(any(isnan(x), [1, 2])), Xboth, 'UniformOutput', false);
intertrls = trls{1} | trls{2};
idxtmp = find(idx);
idx(idxtmp(intertrls)) = 0;
Xboth = cellfun(@(x) x(:, :, ~intertrls), Xboth, 'UniformOutput', false);

% Initialize variables to store data for the table
acc = cell(length(SIDs), 1);
AUC = cell(length(SIDs), 1);
n_trials = nan(length(SIDs), 1);
n_classes = nan(length(SIDs), 1);
SID_list = cell(length(SIDs), 1);
elecs_list = cell(length(SIDs), 1);
confmat = cell(length(SIDs), 1);
confmat_ord = cell(length(SIDs), 1);
tps_list = cell(length(SIDs), 1);
lang_list = nan(length(SIDs), 1);
for lang = 1:length(SIDs)

    y = Svow.vowel(idx);
    [X, ~, ~, y, unisids] = makeDataMatrix(Xboth{lang}, y, sidsboth{lang}, 2000, 4, 1);
    % 950 if you require vowels to be stressed

    % accuracy is heavily impacted by distribution, make sure trials
    % between groups are largely the same
    sidstr = join(unisids(:), ' ');
    disp(['sids used...' sidstr{:}]);
    disp(['num trials used...' num2str(length(y))]);

    rng(2);
    [~, ~, auc, ~, ~, acc_tmp, ~, y_hat] = lda(X', y, 1, [], tps, nfolds);

    [C,order] = confusionmat(y, y_hat);

    % Store data for the table
    acc(lang) = {acc_tmp};
    AUC(lang) = {auc};
    n_trials(lang) = length(y);
    n_classes(lang) = length(unique(y));
    SID_list{lang} = unisids;
    tps_list{lang} = tps;
    lang_list(lang) = lang;
    confmat(lang) = {C};
    confmat_ord(lang) = {order};

    % display mean AUC
    disp(['ls = ' num2str(lang) ' mean AUC = ' num2str(mean(AUC{lang}))]);
end

% Create a table from the collected data
vowel_decode_table = table(SID_list, acc, AUC, n_trials, n_classes, tps_list, lang_list, ...
    confmat, confmat_ord, ...
    'VariableNames', {'SID', 'Accuracy', 'AUC', 'NumTrials', 'NumClasses', ...
    'TimePoints', 'Language', 'ConfMat', 'ConfMatOrd'});

filename = [corpus '_vowel_decode_' subj '_' timelabel '_pooled.mat']; % dimex filename
save([datapath 'ecog_decode/vowelCategory/' filename], 'vowel_decode_table');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% Visualize pooled results from both corpora together
  
figure;
subj = 'monolingual';
timelabel = '100-400ms'; % after onset
corpora = {'dimex', 'timit'};
field = 'AUC';
for ctr = 1:2
    corpus = corpora{ctr};
    subplot(1, 2, ctr);

    % load in data
    load([datapath 'ecog_decode/vowelCategory/' corpus ...
        '_vowel_decode_' subj '_' timelabel '_pooled.mat']);

    % plot the median AUC for each subject, colored by language
    spkrs = {'Spanish', 'English'};
    for s = 1:height(vowel_decode_table)
        % numel = length(vowel_decode_table.Electrodes{s});
        rng(2);
        randjitter = @(x, y) y + randn(x, 1)*0.05;
        numfolds = size(vowel_decode_table.(field){s}, 1);
        if vowel_decode_table.Language(s) == 1
            scatter(1+randjitter(numfolds, 0), vowel_decode_table.(field){s}, 100, 'filled', ...
                'MarkerFaceColor', 'b', 'MarkerFaceAlpha', 0.5); hold on;
            boxplot(vowel_decode_table.(field){s}, ...
            'Positions', 1, 'Colors', 'k', 'Labels', spkrs{s});
        else
            scatter(2+randjitter(numfolds, 0), vowel_decode_table.(field){s}, 100, 'filled', ...
                'MarkerFaceColor', 'r', 'MarkerFaceAlpha', 0.5); hold on;
            boxplot(vowel_decode_table.(field){s}, ...
            'Positions', 2, 'Colors', 'k', 'Labels', spkrs{s});
        end        
    end
    
    % look at difference between languages
    [~, p] = ttest2(vowel_decode_table.(field){vowel_decode_table.Language==1}, ...
        vowel_decode_table.(field){vowel_decode_table.Language==2});
    disp(['p = ' num2str(p)]);

    xticks(1:2);
    xticklabels({'Spanish', 'English'});
    xlabel('Language known');
    ylabel(field);
    set(gca, 'FontSize', 15);
    box off;
    
    xlim([0.5 2.5]);
    if strcmp(field, 'Accuracy')
        ylim([0 0.6]);
    else
        ylim([0.5 1]);
    end
    chance = 1/vowel_decode_table.NumClasses(s);
    yline(chance, 'LineStyle', '--', 'LineWidth', 2);
    title(corpus);
end

% show the confusion matrix
for ctr = 1:2
    corpus = corpora{ctr};
    subplot(1, 2, ctr);

    % load in data
    load([datapath 'ecog_decode/vowelCategory/' corpus ...
            '_vowel_decode_' subj '_' timelabel '_pooled.mat']);
    for s = 1:height(vowel_decode_table)
        imagesc(vowel_decode_table.ConfMat);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;


%% ------- Functions --------

% Two sample ttest over time
function [signtp] = tpttest(data, pthresh, debug)

    % data format each cell is a group with dim: time x reps
    assert(size(data{1}, 1)==size(data{2}, 1));
    tps = 1:size(data{1}, 1);    
    p = nan(1, length(tps));
    for t = tps
        [~, p(t)] = ttest2(data{1}(t, :), ...
        data{2}(t, :));
    end
    signtp = p<pthresh;
    
    % remove all single tps
    singles = strfind(signtp, [0 1 0]);
    for s = singles
        signtp(s+1) = 0;
    end
    
    % Plot ERPs for debugging purposes
    if debug
        figure;
        cols = brewermap(2, 'Set1');
        for p = 1:2
            x = -0.2:0.01:0.5;
            y = mean(data{p}, 2, 'omitnan');
            sem = std(squeeze(data{p}), 0, 2, 'omitnan')/...
                sqrt(size(squeeze(data{p}), 2));
            shadedErrorBar(x, y, sem, 'lineProps', {'markerfacecolor', cols(p, :)} ); 
            hold on;
            scatter(x(signtp), repmat(-0.2, sum(signtp), 1), 35, 'filled', 'r');
        end
    end
    signtp = sum(p<pthresh);
end

% rand index
function [ri] = RandIndex(y, y_pred)
    sC = y == y';
    sC(logical(tril(ones(length(y))))) = [];
    
    sK = y_pred == y_pred';
    sK(logical(tril(ones(length(y_pred))))) = [];
    
    % pairs of elements same sets in C and in the same set in K
    a = sum(sK & sC); 
    
    % pairs of elements different sets in C and in differnt set in K
    b = sum(~sK & ~sC); 
    ri = (a + b)/nchoosek(length(y),2);
    
    % adjusted Ri
    % ari = (ri - exp_ri);
end

% calculate accuracy over varying number of components retained from PCA
function [acc, Y_pred] = accOverComp(A, max_comp, y)

    acc = nan(max_comp, 1);
    Y_pred = nan(size(A, 2), max_comp);
    coeff = pca(A, 'Algorithm', 'eig', 'NumComponents', max_comp);    
    for n = 1:max_comp       
        MdlLinear = fitcdiscr(coeff(:, 1:n), y, 'DiscrimType', 'linear');
        % turns into partitioned model
        [Y_pred(:, n),~, ~] = predict(MdlLinear, coeff(:, 1:n));
        acc(n) = sum(Y_pred(:, n)==y)/length(y);
    end
end

% assumes no valid entry is zero
function [vec] = vectorRDM(rdm)
% process RDM so can be used for correlation analysis
    % ensure it is a dissimilarity matrix
    assert(all(diag(rdm)==0));
    % take only below triangular 
    mask = logical(tril(ones(size(rdm, 1), size(rdm, 2)), -1));
    vec = rdm(mask);
end

function [y_pred, mappedX] = testmapLDA(mdl, X, kfolds)
% uses a pre-trained LDA classifier to predict labels for test set
% and map test set to trained LD space

    dims = 3;
    mappedX = nan(kfolds, size(X, 1), dims);
    allW = nan(kfolds, size(X, 2), dims);
    for k = 1:length(mdl.Trained)   
        [W, LAMBDA] = eig(mdl.Trained{k}.BetweenSigma, ...
            mdl.Trained{k}.Sigma); 
        lambda = diag(LAMBDA);
        [~, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        allLDs = X*W;
        mappedX(k, :, :) = allLDs(:, 1:dims);
        % return corresponding weights for LDs
        allW(k, :, :) = W(:, 1:dims);
    end
    % find average mapping between trained models
    y_pred = ensemPredict(mdl, X, kfolds);

end

function [y_pred] = ensemPredict(partitionMdl, X, kfold)
    for k = 1:kfold
        mdl = partitionMdl.Trained{k};
        [tmp, ~, ~] = mdl.predict(X);
        y(k, :) = tmp;
    end
    y_pred = mode(y);
end

function [colors] = getColors(type, numel)
    switch type
        case 1 % language type
            % spanish, english, mandarin
            colors = brewermap(4, 'Dark2');
        case 2 % spanish language vowels
            colors = brewermap(5, 'Dark2');
            colors(5, :) = [0.2 0.6 0.9];
        case 3 % english language vowels
            colors=(brewermap(numel, 'Dark2'));
            colors(9, :) = [0.2 0.8 0.9];
    end
end

function [trls, selSIDs] = findOverlap(pres, SIDs, minTrl, minSubj)
    % pres as subjects x trls (1 indicates trial was presented to subject)

    varnames = {'SID set', 'trials', 'total'};
    overlapresp = array2table(zeros(0,3), 'VariableNames', varnames);
    % brute force, find as many overlaps as possible
    for k = minSubj:length(SIDs)
        C = nchoosek(1:length(SIDs), k);
        for r = 1:size(C, 1)
            % for this combination of subjects find overlap
            ovrlp = sum(pres(C(r, :), :));
            t2 = table({C(r, :)}, {ovrlp==k}, sum(ovrlp==k), ...
                'VariableNames', varnames);
            overlapresp = [overlapresp; t2];
        end
    end

    maxidx = find(overlapresp.total==max([overlapresp.total]), 1, 'last');
    selSIDs = SIDs(overlapresp.("SID set"){maxidx});
    trls = overlapresp.trials{maxidx};
end

function [dists] = calcClustDist(X, category)
    dists = nan(size(X, 1), 1);
    for cat = unique(category)'
        idx = find(category == cat);
        dists(idx) = pdist2(X(idx, :), median(X(idx, :)));
    end
end

function [C] = mat_prepro(B)
% preprocessing a matrix of form el x time x trials
    
    % find window with greatest variability
    tps = 3:size(B, 2)-2;
    for ctr = 1:length(tps)
        wind = tps(ctr)-2:tps(ctr)+2;
        tp_mean(ctr) = mean(B(:, wind, :), [1, 2, 3], 'omitnan');
    end

    [~, idx] = max(tp_mean);
    B = squeeze(mean(B(:, tps(idx)-2:tps(idx)+2, :), 2, 'omitnan'));
    C = rescale(B,'InputMin', min(B,[],1), 'InputMax', max(B,[],1));
end

function [sort_idx]=findOptSet(conf_acs, coord_acs, type, category, formants)

% finds the optimal set of tokens for acoustic dissimilarity analysis 
% (between English and Spanish)
% sorted_idx is the list of token indices

    % sort acoustic RDM by single col correlation to find set that will make up
    % lowest correlated acoustic RDMs
    cs = cell(2, 1);
    sort_idx = cell(2, 1);
    acs_corr = cell(2, 1);
    labels = {'Spanish', 'English'};
    for i = 1:2
        debug = 1;
        switch type
            case 1 % correlation method, finding least correlated columns of matrix                
                [~, sort_idx{i}] = sort(diag(corr(conf_acs{1, i}, conf_acs{2, i}, ...
                    'Type', 'Spearman')), 'ascend');

            case 2 % procrustes on both acoustic spaces
                y1 = coord_acs{1, i};
                y2 = coord_acs{2, i};
                [~, y2_trans] = procrustes(y1, y2);
                
                if debug
                    disp(['Original correlations: ' mat2str(diag(corr(y1, y2)))]);
                    disp(['New correlations: ' mat2str(diag(corr(y1, y2_trans)))]);
%                     figure;
%                     for x = 1:3
%                         subplot(1, 3, x);
%                         scatter(y1(:,x), y2_trans(:,x), 45, y2(:,x), 'filled');
%                     end
                    
                    ys = {y1, y2_trans};
                    %ys = {y1_form, y2_form};
                    cats = cell(2, 1);
                    figure;
                    for x = 1:2
                        ax(x) = subplot(1, 2, x);
                        cats{x} = arrayfun(@(x) find(x==unique(category{i})), category{i});
                        scatter(ys{x}(:,1), ys{x}(:,2), 45, cats{x}, 'filled');
                        colormap(getColors(i+1, length(unique(category{i}))));     
                        set(gca, 'FontSize', 15);
                    end
                    linkaxes(ax);

                    % comparing to formant space
                    [~, y1_form] = procrustes(formants{i}', y1);
                    [~, y2_form] = procrustes(formants{i}', y2);
                    corr(y1_form, formants{i}');
                    corr(y2_form, formants{i}');

                    figure;
                    violin([silhouette(ys{1}, cats{1}) ...
                        silhouette(ys{2}, cats{2})], 'medc', [], ...
                        'facecolor', [0.7 0.7 0.7], 'bw', [0.2 0.2]);
                    ylabel('silhouette score');
                    xticks([1 2]);
                    xticklabels({'cross-language', 'same language'});
                    set(gca, 'FontSize', 15);
                end
    
%                 [~, sort_idx{i}]=sort(diag(pdist2(y1, ...
%                     y2_trans)), 'ascend');   
                [~, sort_idx{i}]=sort(mean(pdist2(y1, ...
                    y2_trans)), 'ascend');
        end

        % remove first xx rows to see how much correlation drops
        rd1 = conf_acs{1, i}(sort_idx{i}, sort_idx{i});
        rd2 = conf_acs{2, i}(sort_idx{i}, sort_idx{i});           
        for j = 30:size(conf_acs{1, i}, 1)        
            cs{i}(j) = corr(vectorRDM(rd1(1:j, 1:j)), ...
                vectorRDM(rd2(1:j, 1:j)), ...
                'Type', 'Spearman');      
        end        
        sort_idx{i} = sort_idx{i}(1:find(cs{i}<0.3, 1, 'last'));
    
        % calculate correlation of new test set idx
        rd1_ct = rd1(1:length(sort_idx{i}), 1:length(sort_idx{i}));
        rd2_ct = rd2(1:length(sort_idx{i}), 1:length(sort_idx{i}));
        acs_corr{i} = corr(vectorRDM(rd1_ct), vectorRDM(rd2_ct), ...
            'Type', 'Spearman');
    
        % Acoustic correlation of this repetition
        disp([labels{i} ' tokens correlation: ' num2str(acs_corr{i})]); 
        disp(['num trials used: ' num2str(length(sort_idx{i}))]);  
    
        if debug % visualization of analysis step
            figure;     
            subplot(1, 3, 1); imagesc(rd1);
            xline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--'); 
            yline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--');
    
            subplot(1, 3, 2); imagesc(rd2);
            xline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--'); 
            yline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--');
    
            subplot(1, 3, 3); scatter(1:length(cs{i}), cs{i}, ...
                20, 'k', 'filled', 'MarkerFaceAlpha', 0.7); hold on
            ylabel('cumulative correlation between matrices');
    
            yyaxis right;
            scatter(1:length(cs{i}), max([var(rd1); var(rd2)]), ...
                 30, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
            xlim([20 length(cs{i})]);
            ylabel('variance in column');
            xlabel('tokens used');
            xline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--');
        end
    end
end
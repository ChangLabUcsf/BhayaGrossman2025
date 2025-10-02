% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;

% Note - EC202 has no STG coverage
SIDs = [sSIDs, eSIDs]; % , {'HS11', 'HS9', 'HS10'}
tps = 50:55;
load('data/Figure2/Figure2_WordUniVar.mat');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* *modelnames*;

%% A - Example annotated sentences 

% find a sentence you know the words for
wordlists = cellfun(@(x) join(x, ' '), {timit_details.sentdet.wordList});
x = find(strcmp(wordlists, 'the wagons were burning fiercely'));

sents = {42, x};

% Assign corpus details to a cell array
corpus_details = {dimex_details, timit_details};

% Load mexbet2ipa mapping
load('mexbet2ipa.mat');

% Create a new figure
figure;

% Set the x-axis limits for the plots
xl = [0, 1.8];
xlims = {[0.52 1.35], [0.15 1.3]};

% Iterate over the two corpora (DIMEx and TIMIT)
numcorp = 2;
for c = 1:numcorp
    details = corpus_details{c};
    xl = xlims{c};

    disp(details.sentdet(sents{c}).name)
    
    % Iterate over the selected sentences for the current corpus
    for s = 1:length(sents{c})
        % Setup for current sentence
        sent = sents{c}(s);
        sent_info = details.sentdet(sent);

        % Plot amplitude (waveform) of the sound
        subplot(5, numcorp, numcorp*1+c);        
        x = ((1:length(sent_info.sound))/sent_info.soundf)-0.5;
        if startsWith(sent_info.name, 'f')
            factor = 1;
        else
            factor = 3;
        end
        plot(downsample(x, 3), downsample(sent_info.sound, 3)*factor, ...
            'Color', [0.8 0.8 0.8]); hold on;

        % Resample x
        x = ((downsample(1:length(sent_info.sound), ...
            sent_info.soundf./sent_info.dataf))/sent_info.soundf)-0.5;
        offby = length(x)-length(sent_info.loudness);
        x = x(1+ceil(offby/2):end-floor(offby/2));

        % Plot amplitude envelope & syllable onse
        xline(x(sent_info.syltype>0), 'Color', 'k', 'LineStyle', ':', 'LineWidth', 1.5); hold on;
        plot(x, sent_info.loudness, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        
        % Formatting for the current subplot
        xticks(0:2);
        xlim(xl)
        ylim([0 1])
        box off;

        % Plot the spectrogram
        ax = subplot(5, numcorp, c);
        imagesc(x, 1:81, sent_info.aud);
        set(gca, 'YDir', 'normal');
        colormap(ax, flipud(gray));
        xlim(xl)
        xticks(0:2);
        str = join(sent_info.wordList, ' ');
        title([num2str(sent) ': ' str{:}]);   

        % Plot the phonetic features
        ax = subplot(5,numcorp, numcorp*4+c);

        numfeats = length(details.features.names);
        imagesc(x, 1:numfeats, sent_info.phnfeat);
        yticks(1:numfeats);
        yticklabels(details.features.names);
        set(gca, 'YDir', 'normal');

        colormap(ax, [1 1 1; 0 0 0]);
        xlim(xl);
        xticks(0:2);
        box off;
        
        % Plot the surprise 
        subplot(5, numcorp, numcorp*2+c);
        if startsWith(sent_info.name, 'f') % TIMIT
            field = 'engSurp';
        else % DIMEx
            field = 'spSurp';
        end

        if c<3 % only have surprise for Spanish and English
            stem(x(sent_info.(field)>0), sent_info.(field)(sent_info.(field)>0), ...
                'MarkerFaceColor', 'k', 'Color', 'k');
            xlim(xl);
            xticks([]);
            box off;
        end

        % Plot the phonemes
        if startsWith(sent_info.name, 'f') % TIMIT
            load('arpabet2ipa.mat');
            phnidx = cellfun(@(x) find(strcmp(arpabet, x)), ...
                sent_info.phnnames, 'UniformOutput', false);     
        elseif c<3 % DIMEx
            load('mexbet2ipa.mat');
            [phns, ~] = find(sent_info.phnmatonset);
            sent_info.phnnames = arrayfun(@(x) dimex_details.phnnames(x), ...
                phns, 'UniformOutput', false);
            phnidx = cellfun(@(x) find(strcmp(mexbet, x)), sent_info.phnnames, ...
                'UniformOutput', false);
        end

        phnidx(cellfun(@(x) isempty(x), phnidx)) = {-1};
        phnidx = cell2mat(phnidx);

        phns = cell(1, 1);
        phns(phnidx>0) = ipa(phnidx(phnidx>0));
        text(x(mean(sent_info.phnmatonset)>0), -0.8*ones(length(phns), ...
                1), phns);

        % Plot word boundaries
        subplot(5, numcorp, numcorp*3+c);
        %plot(x, sent_info.loudness, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5);
        xline(x(sent_info.wordOns>0), 'LineWidth', 2, 'HandleVisibility', 'off'); hold on;

        % plotting frequency
        wfidx = sent_info.wordL>0;
        scatter(x(wfidx), sent_info.wordL(wfidx), 40, 'k', 'filled', ...
            'DisplayName', 'wordLength', 'Marker', '<');

        wfidx = sent_info.wordFreqLog>0;
        scatter(x(wfidx), sent_info.wordFreqLog(wfidx), 30, 'k', ...
            'filled', 'DisplayName', 'wordFreq');
        xlim(xl);
        box off;
    end   
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% B - Native vs. Foreign: word / sequence surprisal unique variance (box-plots)

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 

thresh = 0.001;

features = {'wordO', 'wordF', 'wordL', 'surp', 'word+surp'}; % word+surp
titles = { 'Boundary',  'Frequency', 'Length', 'Phoneme Surprisal', 'Both'};
fields = {'eng_uv_all', 'sp_uv_all'}; 

ctr = 1;
native = cell(length(features), 1);
sid_native = cell(length(features), 1);
nonnative = cell(length(features), 1);
sid_nonnative = cell(length(features), 1);
sid = cell(length(features), 1);
el = cell(length(features), 1);
hemi = cell(length(features), 1);

for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);
    
    for ls = [1 2 4]
        % Extract subject id and electrode for linear mixed effect model
        if ismember(ls, [1, 2])
            sid{ctr} = [sid{ctr}; wordsurp_encoding.SID(wordsurp_encoding.ls==ls)];
            el{ctr} = [el{ctr}; wordsurp_encoding.el(wordsurp_encoding.ls==ls)];
            hemi{ctr} = [hemi{ctr}; wordsurp_encoding.hemi(wordsurp_encoding.ls==ls)];
        end

        for f = 1:2 % Aggregate over English and Spanish fields (English is first)
            field = fields{f};    
            y = wordsurp_encoding.(field)(:, index);   
            sid_tmp = wordsurp_encoding.SID;   
            
            % Removing bilinguals so the comparison is more balanced
            if (f == 1 && ismember(ls, 2)) || (f == 2 && ismember(ls, 1))
                native{ctr} = [native{ctr}; y(wordsurp_encoding.ls==ls)];  
                sid_native{ctr} = [sid_native{ctr}; sid_tmp(wordsurp_encoding.ls==ls)];
            elseif (f == 1 && ismember(ls, 1)) || (f == 2 && ismember(ls, 2))
                nonnative{ctr} = [nonnative{ctr}; y(wordsurp_encoding.ls==ls)];
                sid_nonnative{ctr} = [sid_nonnative{ctr}; sid_tmp(wordsurp_encoding.ls==ls)];
            end
        end
    end
    ctr = ctr + 1;
end

% Combine the consonant and vowel features into one
subplts = {1, 2, 3, 4}; % add 5 here to see word+surprisal
ynat_pie = nan(length(subplts), 2);
ynon_pie = nan(length(subplts), 2);

% split by sid
ynat_sid = cell(length(subplts), 2);
ynon_sid = cell(length(subplts), 2);

figure('Renderer', 'Painters');
for s = 1:length(subplts)
    
    % Combining native and unfamiliar subjects
    % subplot(1, length(subplts), s)
    rep = length(subplts{s});
    
    % Native boxplot
    y_nat = arrayfun(@(x) native{x}, [subplts{s}], 'UniformOutput', false);
    y_nat = cat(1, y_nat{:});
    ynat_pie(s, :) = [sum(y_nat<thresh); sum(y_nat>thresh)];

    % Per subject, how many electrodes meet threshold?
    y_sid = arrayfun(@(x) sid_native{x}, [subplts{s}], 'UniformOutput', false);
    y_sid = repmat(y_sid{1}, rep, 1);

    prct = cellfun(@(x) sum(y_nat>thresh & strcmp(y_sid, x)) ...
        ./sum(strcmp(y_sid, x)), unique(y_sid));
    cnt = cellfun(@(x) sum(strcmp(y_sid, x)), unique(y_sid));
    ynat_sid{s, 1} = [prct, cnt];
    clear y_sid prct cnt

    h=boxchart(ones(sum(y_nat>thresh), 1)+(s-1)*1.5, y_nat(y_nat>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.9], 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.3); hold on;
                    h.JitterOutliers = 'on';
                    h.MarkerStyle = '.';
                    h.MarkerColor = 'k';

    % Unfamiliar boxplot
    y_non = arrayfun(@(x) nonnative{x}, [subplts{s}], 'UniformOutput', false);
    y_non = cat(1, y_non{:});
    ynon_pie(s, :) = [sum(y_non<thresh); sum(y_non>thresh)];

    % Per sid, how many electrodes meet threshold?
    y_sid = arrayfun(@(x) sid_nonnative{x}, [subplts{s}], 'UniformOutput', false);
    y_sid = repmat(y_sid{1}, rep, 1);

    prct = cellfun(@(x) sum(y_non>thresh & strcmp(y_sid, x)) ...
        ./sum(strcmp(y_sid, x)), unique(y_sid));
    cnt = cellfun(@(x) sum(strcmp(y_sid, x)), unique(y_sid));
    ynon_sid{s} = [prct, cnt];
    clear y_sid prct cnt

    h=boxchart(ones(sum(y_non>thresh), 1)*1.5+(s-1)*1.5, y_non(y_non>thresh), ...
        'BoxFaceColor', [0.9 0.5 0.5], 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.3);
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([y_nat(y_nat>thresh); y_non(y_non>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    % xticks([1 2]);
    % xticklabels(conditionLabels);
    xlim([0.5 2.5]);
    ylim([0 0.05]);
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    idx = y_nat>thresh & y_non>thresh;
    tbl.rsq = [y_nat(idx); y_non(idx)];
    
    elecs = repmat(el{1}, rep, 1);
    sids = repmat(sid{1}, rep, 1);
    lss = cellfun(@(x) find(cellfun(@(y) ...
        ismember(x, y), {sSIDs, eSIDs, bSIDs})), sids);
    hemis = repmat(hemi{1}, rep, 1);

    tbl.sid = repmat(sids(idx), 2, 1);
    tbl.elec = repmat(elecs(idx), 2, 1);
    tbl.hemi = repmat(hemis(idx), 2, 1);
    tbl.ls = repmat(lss(idx), 2, 1);
    tbl.native = [ones(sum(idx), 1); ones(sum(idx), 1)*2];
    clear elecs sids

    lme2 = fitlme(tbl,'rsq~native+ls+(1|hemi)+(1|sid)+(1|elec:sid)');
    disp(titles{s})
    disp(lme2)
    p = lme2.Coefficients.pValue(3);

    disp(['LME native language p-value = ' num2str(p)]);
    
    % line([1, 2], [maxy(2)-0.1 maxy(2)-0.1], 'Color', 'k');
    if isempty(getSigStr(p, 2))
        %text(1, maxy(2)-0.02, 'n.s.', 'FontSize', 15);
        text((s-1)*1.5 + 1, maxy(2)-0.02, 'n.s.', 'FontSize', 25);
    else
        %text(1, maxy(2)-0.02, getSigStr(p, 1), 'FontSize', 15);
        text((s-1)*1.5 + 1, maxy(2)-0.02, getSigStr(p, 1), 'FontSize', 25);
    end
    %title(titles{s},  'FontWeight', 'normal');
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2 (log)');
    legend({'native', 'foreign'});
end
xlim([0.5 1.5*length(subplts)+0.5]);
xticks((1:length(subplts))*1.5-0.25);
xticklabels(titles);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% C - Native vs. Foreign: acoustic phonetic unique variance (box-plots)

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 

thresh = 0.001;

features = {'peakrate', 'formant', 'consonant'};
titles = {'PeakRate', 'Vowel Formants', 'Consonant'};

% features = {'pitch', 'env'};
% titles = {'Pitch', 'Envelope'};
% features = {'bisurp', 'trisurp'};
% titles = {'Biphone', 'Triphone'};
fields = {'eng_uv_all', 'sp_uv_all'}; 

ctr = 1;
native = cell(length(features), 1);
sid_native = cell(length(features), 1);
nonnative = cell(length(features), 1);
sid_nonnative = cell(length(features), 1);
sid = cell(length(features), 1);
el = cell(length(features), 1);
hemi = cell(length(features), 1);

for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);
    
    for ls = [1 2 4]
        % Extract subject id and electrode for linear mixed effect model
        if ismember(ls, [1, 2])
            sid{ctr} = [sid{ctr}; wordsurp_encoding.SID(wordsurp_encoding.ls==ls)];
            el{ctr} = [el{ctr}; wordsurp_encoding.el(wordsurp_encoding.ls==ls)];
            hemi{ctr} = [hemi{ctr}; wordsurp_encoding.hemi(wordsurp_encoding.ls==ls)];
        end

        for f = 1:2 % Aggregate over English and Spanish fields
            field = fields{f};    
            y = wordsurp_encoding.(field)(:, index);   
            sid_tmp = wordsurp_encoding.SID;   
            
            % Removing bilinguals so the comparison is more balanced
            if (f == 1 && ismember(ls, 2)) || (f == 2 && ismember(ls, 1))
                native{ctr} = [native{ctr}; y(wordsurp_encoding.ls==ls)];  
                sid_native{ctr} = [sid_native{ctr}; sid_tmp(wordsurp_encoding.ls==ls)];
            elseif (f == 1 && ismember(ls, 1)) || (f == 2 && ismember(ls, 2))
                nonnative{ctr} = [nonnative{ctr}; y(wordsurp_encoding.ls==ls)];
                sid_nonnative{ctr} = [sid_nonnative{ctr}; sid_tmp(wordsurp_encoding.ls==ls)];
            end
        end
    end
    ctr = ctr + 1;
end

% Combine the consonant and vowel features into one
subplts = cell(1, length(features));
for i = 1:length(subplts)
    subplts{i} = i;
end
ynat_pie = nan(length(subplts), 2);
ynon_pie = nan(length(subplts), 2);

% split by sid
ynat_sid = cell(length(subplts), 2);
ynon_sid = cell(length(subplts), 2);

figure('Renderer', 'Painters');
for s = 1:length(subplts)
    
    % Combining native and unfamiliar subjects
    % subplot(1, length(subplts), s)
    rep = length(subplts{s});
    
    % Native boxplot
    y_nat = arrayfun(@(x) native{x}, [subplts{s}], 'UniformOutput', false);
    y_nat = cat(1, y_nat{:});
    ynat_pie(s, :) = [sum(y_nat<thresh); sum(y_nat>thresh)];

    % Per subject, how many electrodes meet threshold?
    y_sid = arrayfun(@(x) sid_native{x}, [subplts{s}], 'UniformOutput', false);
    y_sid = repmat(y_sid{1}, rep, 1);

    prct = cellfun(@(x) sum(y_nat>thresh & strcmp(y_sid, x)) ...
        ./sum(strcmp(y_sid, x)), unique(y_sid));
    cnt = cellfun(@(x) sum(strcmp(y_sid, x)), unique(y_sid));
    ynat_sid{s, 1} = [prct, cnt];
    clear y_sid prct cnt

    h=boxchart(ones(sum(y_nat>thresh), 1)+(s-1)*1.5, y_nat(y_nat>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.9], 'MarkerColor', 'k', 'Notch','on', ...
        'BoxWidth', 0.3); hold on;
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Unfamiliar boxplot
    y_non = arrayfun(@(x) nonnative{x}, [subplts{s}], 'UniformOutput', false);
    y_non = cat(1, y_non{:});
    ynon_pie(s, :) = [sum(y_non<thresh); sum(y_non>thresh)];

    % Per sid, how many electrodes meet threshold?
    y_sid = arrayfun(@(x) sid_nonnative{x}, [subplts{s}], 'UniformOutput', false);
    y_sid = repmat(y_sid{1}, rep, 1);

    prct = cellfun(@(x) sum(y_non>thresh & strcmp(y_sid, x)) ...
        ./sum(strcmp(y_sid, x)), unique(y_sid));
    cnt = cellfun(@(x) sum(strcmp(y_sid, x)), unique(y_sid));
    ynon_sid{s} = [prct, cnt];
    clear y_sid prct cnt

    h=boxchart(ones(sum(y_non>thresh), 1)*1.5+(s-1)*1.5, y_non(y_non>thresh), ...
        'BoxFaceColor', [0.9 0.5 0.5], 'MarkerColor', 'k', 'Notch','on', ...
        'BoxWidth', 0.3);
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([y_nat(y_nat>thresh); y_non(y_non>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    % xticks([1 2]);
    % xticklabels(conditionLabels);
    % xlim([0.5 2.5]);
    ylim([0 0.1]);
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    idx = y_nat>thresh & y_non>thresh;
    tbl.rsq = [y_nat(idx); y_non(idx)];
    
    elecs = repmat(el{1}, rep, 1);
    sids = repmat(sid{1}, rep, 1);
    lss = cellfun(@(x) find(cellfun(@(y) ...
        ismember(x, y), {sSIDs, eSIDs, bSIDs})), sids);
    hemis = repmat(hemi{1}, rep, 1);

    tbl.sid = repmat(sids(idx), 2, 1);
    tbl.elec = repmat(elecs(idx), 2, 1);
    tbl.hemi = repmat(hemis(idx), 2, 1);
    tbl.ls = repmat(lss(idx), 2, 1);
    tbl.native = [ones(sum(idx), 1); ones(sum(idx), 1)*2];
    clear elecs sids

    lme2 = fitlme(tbl,'rsq~native+ls+(1|hemi)+(1|sid)+(1|elec:sid)');
    disp(titles{s})
    disp(lme2)
    p = lme2.Coefficients.pValue(3);

    disp(['LME native language p-value = ' num2str(p)]);
    
    % line([1, 2], [maxy(2)-0.1 maxy(2)-0.1], 'Color', 'k');
    if isempty(getSigStr(p, 2))
        %text(1, maxy(2)-0.02, 'n.s.', 'FontSize', 25);
        text((s-1)*1.5 + 1, maxy(2)-0.02, 'n.s.', 'FontSize', 20);
    else
        %text(1, maxy(2)-0.02, getSigStr(p, 1), 'FontSize', 15);
        text((s-1)*1.5 + 1, maxy(2)-0.02, getSigStr(p, 1), 'FontSize', 20);
    end
    % title(titles{s},  'FontWeight', 'normal');
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2 (log)');
end
xlim([0.5 1.5*length(subplts)+0.5]);
xticks((1:length(subplts))*1.5-0.25);
xticklabels(titles);
legend({'native', 'foreign'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% D - MNI Brain with Native / Unfamiliar word / surprisal electrodes highlighted (brain)

% initialize design electrode structure
fields = {'sp_uv_all', 'eng_uv_all'}; 
uv_thresh = 0;
% pval_thresh = 0.05;

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% uv feature order
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant',

figure();
for native = 0:1
    for f = 1:length(feats)
        feat = feats{f};
        index = find(ismember(wordsurp_details.featureOrd, feat));
    
        % Change the color depending on the feature 
        switch feat
            case 'peakrate'
                featcol = [0.3 0.8 0.9];
            case 'formant'
                featcol = [0.4 0.7 0.3];
            case 'consonant'
                featcol = [0.6 0.0 0.7];
            case 'word+surp'
                featcol = [0.7 0.1 0.2];
            case 'surp'
                featcol = [0.7 0.1 0.2];
            otherwise
                featcol = [0 0 0];
        end   
    
        % make desel structure
        desel=struct();
        desel.conds = 1:7;
        %ls = [1, 2]; % can only do 1-2
        
        % size and color
        % 1:20:200; %ones(1, 10)*0.00000001; %1
        desel.sz = [2; 35*ones(length(desel.conds), 1)]; 
        desel.sz = [5; 10; 25*desel.conds']; 
        
        % split up peak rate and phonetic features again for MNI plotting
        desel.labels = [];

        % use the native combination (ls = 1 & sp, ls = 2 & eng)
        desel.yval = nan(height(wordsurp_encoding.ls), 1);
        if native
            % Spanish natives
            sidx = find(wordsurp_encoding.ls==1);
            desel.yval(sidx) = arrayfun(@(x) wordsurp_encoding.(fields{1})(x, index), sidx);

            % English natives
            eidx = find(wordsurp_encoding.ls==2);
            desel.yval(eidx) = arrayfun(@(x) wordsurp_encoding.(fields{2})(x, index), eidx);
        else % unfamiliar case
            % Spanish natives
            sidx = find(wordsurp_encoding.ls==1);
            desel.yval(sidx) = arrayfun(@(x) wordsurp_encoding.(fields{2})(x, index), sidx);

            % English natives
            eidx = find(wordsurp_encoding.ls==2);
            desel.yval(eidx) = arrayfun(@(x) wordsurp_encoding.(fields{1})(x, index), eidx);
        end
        
        % discretize values
        % 0.005:0.005:0.01
        % binedges = [-1 uv_thresh:0.005:0.01 0.015:0.015:0.1];
        binedges = [-1 uv_thresh:0.01:0.1];
        for s=unique(wordsurp_encoding.SID)'
            SID = s{1};
            idx = strcmp(wordsurp_encoding.SID, SID);
            desel.(SID).elid = wordsurp_encoding.el(idx);
            desel.(SID).condition = discretize(desel.yval(idx), ...
                binedges);
        end

        if native
            cls = flipud(blues(8));
        else
            cls = flipud(reds(8));
        end
        desel.cols = [0.3 0.3 0.3; cls(3:end, :)];

        figure;
        colormap(desel.cols);
        colorbar;
        
        mni_lh = plotMNIElec(unique(wordsurp_encoding.SID), desel, 'lh', 0, 1, imgall);
        
        sgtitle(native);
        l = light;
        view(270, 0);   
        set(l,'Style', 'infinite', 'Position', [-1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.85;
        % add a pie
        % axes('Position',[.6 .15 .3 .3])
        % p = pie([sum(mni_lh.cond>1), sum(mni_lh.cond==1)], [1 1]); 
        % p(1).FaceColor = [desel.cols(5, :)];
        % p(1).EdgeColor = 'none';
        % p(3).FaceColor = [0.6 0.6 0.6];
        % p(3).EdgeColor = 'none';
        % p(2).Color = 'w';
        % p(2).FontWeight = 'bold';
        % p(2).FontSize = 13;
        % p(4).FontWeight = 'bold';
        % p(4).Color = 'w';
        % p(4).FontSize = 13;

        mni_rh = plotMNIElec(unique(wordsurp_encoding.SID), desel, 'rh', 0, 1, imgall);
        sgtitle(native);
        l = light;
        view(90, 0);
        set(l,'Style', 'infinite', 'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.85;
     end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% E - Brain counts with Native / Unfamiliar word / surprisal electrodes highlighted (brain)

pthresh = 0.05;
% make empty tables with all these at 0
anat_labels = {'superiortemporal', 'postcentral', ...
    'middletemporal', 'supramarginal', 'precentral'};
tbl = table([0 0 0 0 0]', anat_labels', 'VariableNames', {'count', 'area'});
anat_counts = {tbl, tbl};

for native = [0, 1]
    for ls = 1:2
        if (ls == 1 && ~native) || (ls == 2 && native) 
            idx = wordsurp_encoding.ls==ls ...
                & wordsurp_encoding.eng_wordsurp_pval<pthresh;
        else
            idx = wordsurp_encoding.ls==ls ...
                    & wordsurp_encoding.sp_wordsurp_pval<pthresh;
        end
        [counts, ~, ~, labels] = crosstab(wordsurp_encoding.anatomy(idx));

        % if anat_counts{native+1} is already a table, add to existing, otherwise make a new table
        if istable(anat_counts{native+1})
            % find equivalent area in table and add count
            for i = 1:length(labels)
                area_idx = find(strcmp(anat_counts{native+1}.area, labels(i)));
                if ~isempty(area_idx)
                    anat_counts{native+1}.count(area_idx) = ...
                        anat_counts{native+1}.count(area_idx) + counts(i);
                else
                    anat_counts{native+1} = [anat_counts{native+1}; ...
                        table(counts(i), labels(i), 'VariableNames', ...
                        {'count', 'area'})];
                end
            end
        else
            anat_counts{native+1} = table(counts, labels, ...
                'VariableNames', {'count', 'area'});
        end
    end
end

% Create horizontal bar charts for speech-responsive areas
figure;
subplot(1, 2, 1);

% remove areas with count < 5 and sort by count
% anat_counts{1} = anat_counts{1}(anat_counts{1}.count >= 2, :);
% anat_counts{1} = sortrows(anat_counts{1}, 'count', 'descend');
sortidx = cellfun(@(x) find(ismember([anat_counts{1}.area], x)), ...
    anat_labels);
barh(anat_counts{1}.area(sortidx), anat_counts{1}.count(sortidx), 'FaceColor', ...
    [1, 0, 0], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Red color for non-native
title('Foreign');
xlabel('Count');
ylabel('Anatomical Area');
set(gca, 'FontSize', 14);
xlim([0 120]);
xticks([0 100]);

subplot(1, 2, 2);
% remove areas with count < 5 and sort by count
% anat_counts{2} = anat_counts{2}(anat_counts{2}.count >= 2, :);
% anat_counts{2} = sortrows(anat_counts{2}, 'count', 'descend');
sortidx = cellfun(@(x) find(ismember([anat_counts{2}.area], x)), ...
    anat_labels);
barh(anat_counts{2}.area(sortidx), anat_counts{2}.count(sortidx), 'FaceColor', ...
    [0, 0, 1], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Blue color for native
title('Native');
xlabel('Count');
ylabel('Anatomical Area');
set(gca, 'FontSize', 14);
sgtitle('Speech-responsive areas');
xlim([0 155]);
xticks([0 100]);

figure;
% same thing but together for native and foreign (one plot)
y_positions = 1:length(sortidx);
bar_width = 0.35;  % Width of each bar

% Plot bars side by side
sortidx = cellfun(@(x) find(ismember([anat_counts{1}.area], x)), ...
    anat_labels);
barh(y_positions - bar_width/2, anat_counts{1}.count(sortidx), bar_width, ...
    'FaceColor', [1, 0, 0], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Red for non-native
hold on;
sortidx = cellfun(@(x) find(ismember([anat_counts{2}.area], x)), ...
    anat_labels);
barh(y_positions + bar_width/2, anat_counts{2}.count(sortidx), bar_width, ...
    'FaceColor', [0, 0, 1], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Blue for native

% Set y-axis labels
yticks(y_positions);
yticklabels(anat_counts{1}.area(sortidx));

% add legend
legend('Foreign', 'Native');
title('Word-level areas');

% add xlabel
xlabel('Count');
xticks([0 100]);
xlim([0 155]);
yticklabels({'superiortemp', 'postcent', ...
    'middletemp', 'supramarg', 'precent'});
set(gca, 'FontSize', 14);
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% F - Percentage of word-level responses across native and foreign

% electrode selection details
timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");
asccd_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat");

figure;
pthresh = 0.05;
native_elecs = {dimex_elecs, timit_elecs, asccd_elecs}; % Spanish, English, Mandarin
foreign_elecs = {timit_elecs, dimex_elecs, timit_elecs}; % English, Spanish, English
native_fields = {'sp_wordsurp_pval', 'eng_wordsurp_pval'};
foreign_fields = {'eng_wordsurp_pval', 'sp_wordsurp_pval'};
prct_word = nan(2, 2);
for ls = [1, 2] % languages : spanish, english, mandarin
    
    % find number of speech responsive electrodes in table
    idx = sent_encoding.ls == ls;
    sids = unique(sent_encoding.SID(idx));
    % 1st column : native, 2nd column : foreign
    numspeech = zeros(length(sids), 2);
    numword = zeros(length(sids), 2);

    % for each subject, identify the number of total elecs
    ctr = 1;
    for s = sids'
        sid = s{1};
        numspeech(ctr, :) = [length(native_elecs{ls}.allidx.(sid)), ...
            length(foreign_elecs{ls}.allidx.(sid))];
        nativeidx = strcmp(wordsurp_encoding.SID, sid) & ...
            wordsurp_encoding.(native_fields{ls})<pthresh;
        foreignidx = strcmp(wordsurp_encoding.SID, sid) & ...
            wordsurp_encoding.(foreign_fields{ls})<pthresh;
        numword(ctr, :) = [sum(nativeidx), sum(foreignidx)];
        clear nativeidx foreignidx;
        ctr = ctr + 1;
    end
    prct_word(ls, 1) = sum(numword(1, :))/ sum(numspeech(1, :));
    prct_word(ls, 2) = sum(numword(2, :))/ sum(numspeech(2, :)); 
end

% make a figure for pie chart combining across languages
figure;
subplot(1, 2, 1);
pc = pie([mean(prct_word(:, 1)), 1-mean(prct_word(:, 1))], [0 1]);
pc(1).FaceColor = [0.6 0.6 1];
pc(3).FaceColor = [0.8 0.8 0.8];
pc(1).EdgeColor = 'none';
pc(3).EdgeColor = 'none';
pc(2).FontSize = 14;
pc(4).FontSize = 14;

subplot(1, 2, 2);
pc = pie([mean(prct_word(:, 2)), 1-mean(prct_word(:, 2))], [0 1]);
pc(1).FaceColor = [1 0.6 0.6];
pc(3).FaceColor = [0.8 0.8 0.8];
pc(1).EdgeColor = 'none';
pc(3).EdgeColor = 'none';
pc(2).FontSize = 14;
pc(4).FontSize = 14;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% G - Native vs. Foreign: word and sequence surprisal unique variance (scatter)

% Specify the labels for analysis
labels = {'word+surp'};

% Create a figure for the scatter plots
figure('Position', [200, 100, 900, 900]);

% Set a UV threshold for filtering
pval_thresh = 0.05;

% Initialize subplot counter and axis array
ctr = 1;
ax = nan(3, 1);

cols = [256 256 256; 28 117 188; 237 28 36; 200 135 236]./256;
% make a patch in the upper right quadrant thats light purple
%patch([0 0.08 0.08 0], [0 0 0.08 0.08], cols(4, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35); hold on;
% make a patch in the lower right quadrant thats light red
patch([0 0 0.08 0.08], [-0.01 0.08 0.08 -0.01], cols(3, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);
%patch([0 0 0.08 0.08], [-0.01 0 0 -0.01], cols(3, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);
% make a patch in the upper left quadrant thats light blue
patch([0.08 0.08 -0.01 -0.01], [0 0.08 0.08 0], cols(2, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35); hold on;
%patch([0 0 -0.01 -0.01], [0 0.08 0.08 0], cols(2, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);

% Iterate over the labels
for label = labels
    x_all = [];
    y_all = [];
    
    % Find the index of the current label in the featureOrd array
    index = find(ismember(wordsurp_details.featureOrd, label));
    
    % Iterate over the ls values
    for ls = [1, 2]
        % Create a subplot for the current label and ls value
        ax(ctr) = subplot(length(labels), 1, ctr);
        lsidx = wordsurp_encoding.ls==ls;
        
        % Get the corresponding unique variance values for English and Spanish
        if ls == 1
            x = wordsurp_encoding.eng_uv_all(lsidx, index);
            y = wordsurp_encoding.sp_uv_all(lsidx, index);
            pvals_x = wordsurp_encoding.eng_wordsurp_pval(lsidx);
            pvals_y = wordsurp_encoding.sp_wordsurp_pval(lsidx);
        else
            x = wordsurp_encoding.sp_uv_all(lsidx, index);
            y = wordsurp_encoding.eng_uv_all(lsidx, index);
            pvals_x = wordsurp_encoding.sp_wordsurp_pval(lsidx);
            pvals_y = wordsurp_encoding.eng_wordsurp_pval(lsidx);
        end
        
        % Get the subject IDs for the current ls value
        sid = cellfun(@(x) str2double(x(3:end)), ...
            wordsurp_encoding.SID(wordsurp_encoding.ls==ls));
        
        % Remove data points that do not meet the UV threshold or contain NaN values
        % all([x,y]<uv_thresh, 2)
        neg = all([pvals_x,pvals_y]>pval_thresh, 2) | ...
            all([x,y]<0, 2) | isnan(x) | isnan(y);
        %neg = all([x,y]<uv_thresh, 2) | isnan(x) | isnan(y);
        x(neg) = [];
        y(neg) = [];
        sid(neg) = [];

        % Perform permutation testing to compute correlation coefficient and p-value
        maxlim = prctile([x; y], 100);
        minlim = prctile([x; y], 3);  

        % Plot the scatter plot
        colors = x-y;
        x_all = [x_all; x];
        y_all = [y_all; y];
        %scatter3(x, y, sid, 20, colors, 'filled', ...
        %            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8); hold on;
        
        scatter(x, y, 15, colors, 'filled', ...
                    'MarkerEdgeColor', 'k', ...
                    'MarkerFaceAlpha', 0.8, ...
                    'LineWidth', 0.25); hold on;
        
        % colorbar should be blue to red going through white (1, 1, 1)
        % first create blue to white
        cols = [linspace(43/256, 1, 50); linspace(122/256, 1, 50); linspace(186/256, 1, 50)]';
        % then add white to red
        cols = [cols; [linspace(1, 213/256, 50); linspace(1, 37/256, 50); ...
            linspace(1, 39/256, 50)]'];

        colormap(cols);
        clim([-0.05 0.05]);
        view(2);
        
        % Set the x and y axis limits and add labels
        xlim([minlim maxlim]);
        ylim([minlim maxlim]);
        xlabel(['Foreign ' label{1} ' \Delta R^2']);
        ylabel(['Native ' label{1} ' \Delta R^2']);
    end

    [r, p] = corr(x_all(x_all>0&y_all>0), y_all(x_all>0&y_all>0), 'Rows', ...
        'complete', 'type', 'Spearman');
    title({['r(' num2str(sum(x_all>0&y_all>0)) ')= ' num2str(r) ','], ['p=' num2str(p, 4)]});
     % Increment the subplot counter
    ctr = ctr + 1; 
    hold on;
end

% Link axes for phonetic features and word onset/surp
linkaxes(ax(1:length(labels)));
h = refline(1, 0);
set(gca, 'FontSize', 15);

% Add reference lines and lines at 0
for i = 1:length(labels)
    subplot(length(labels), 1, i);
    xline(0, 'Color', 'k', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
    xlim([-0.01 0.08]);
    ylim([-0.01 0.08])
    xticks([0 0.08]);
    yticks([0 0.08]);
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    % make colors all black
    colormap(repmat([0 0 0], 256, 1));
end

% add inlaid quadrant count on the top right corner
% add a inlaid quadrant plot
ax = axes('Position',[.65 .65 .3 .3]);
quads = rot90(histcounts2(x_all, y_all, [-100 0 100], [-100 0 100]));

imagesc([1, 3; 0, 2]);
% change the colormap of this axis
cols = [256 256 256; 28 117 188; 237 28 36; 200 135 236]./256;
cols = brighten(cols, 0.8);
colormap(ax, cols);

% colormap(ax, flipud(gray));
%colormap([1 1 1; 1 0 0; 1 0 1;0 0 1])
% add in the text overlaid
quads = flipud(rot90(quads));
for x = 1:2
    for y = 1:2
        text(x, y, num2str(quads(x, y)), 'Color', 'w', ...
            'FontSize', 16, 'HorizontalAlignment', 'center');
    end
end
box off;
xticks([1 2]);
yticks([1 2]);
xticklabels({'-UV', '+UV'});
yticklabels({'+UV', '-UV'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% ----------------------- Supplementary Figures --------------------------
%% S7 - Correlating word and sequence surprise in native speech condition (scatter)

feats = {'surp', 'wordO', 'wordF', 'wordL'};
fields = {'sp_uv_all', 'eng_uv_all'};
featnames = {'surprisal', 'onset', 'frequency', 'length'};

% Iterate over the features combinations
figure;
ctr = 1;
for f= 1:length(feats)
    for j=f+1:length(feats)
        % Find the index of the current label in the featureOrd array
        index = [];
        index(1) = find(ismember(wordsurp_details.featureOrd, feats{f}));
        index(2) = find(ismember(wordsurp_details.featureOrd, feats{j}));
        
        % Get the corresponding unique variance values for English and Spanish
        x_all = zeros(height(wordsurp_encoding.ls), 1);
        y_all = zeros(height(wordsurp_encoding.ls), 1);
        for ls = 1:2
            uvs = wordsurp_encoding.(fields{ls});
            nativeidx = wordsurp_encoding.ls==ls;
            x_all(nativeidx) = uvs(nativeidx, index(1));
            y_all(nativeidx) = uvs(nativeidx, index(2));
        end

        % Remove data points that contain NaN values
        neg = isnan(x_all) | isnan(y_all) | (x_all<0 & y_all<0);
        x_all(neg) = [];
        y_all(neg) = [];

        subplot(2, 3, ctr);
        andpos = x_all>0 & y_all>0;
        scatter(x_all, y_all, 15, 'k', 'filled', 'MarkerFaceAlpha', 0.8); hold on;
        ylim([-0.01 0.05]);
        xlim([-0.01 0.05]);
        yticks(0:0.025:0.05);
        xticks(0:0.025:0.05);

        % add best fit line
        h = lsline;
        h.Color = 'r';
        h.LineWidth = 2;

        orpos = x_all>0 & y_all<0 | x_all<0 & y_all>0;
        scatter(x_all(orpos), y_all(orpos), 10, [0.5 0.5 0.5], ...
            'filled','MarkerFaceAlpha', 0.8);
        ylabel([featnames{j} ' \Delta R^2']);
        xlabel([featnames{f} ' \Delta R^2']);

        % add a reference line
        xline(0, 'Color', 'k', 'LineWidth', 1.5);
        yline(0, 'Color', 'k', 'LineWidth', 1.5);

        % Perform permutation testing to compute correlation coefficient and p-value
        [r, p] = corr(x_all, y_all, 'Rows', 'complete', 'Type', 'Spearman');
        title(['r(' num2str(length(x_all)) ')= ' num2str(r, 2) ', ' getSigStr(p, 2)], 'FontWeight', 'normal');
        ctr = ctr + 1;
    end
end

% xlim([-0.005 0.03]);
% xticks(0:0.01:0.03);
% ylim([-0.005 0.07]);
% yticks(0:0.02:0.06);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% S8 - Native vs. Foreign: word + sequence surprisal unique variance (by subject box-plots)

pthresh = 0.01;
numsubj = length(unique(wordsurp_encoding.SID));
native_unfamiliar_prct = nan(3, numsubj);

ctr=1;
lss = nan(numsubj, 1);
for s = unique(wordsurp_encoding.SID)'
    SID = s{1};
    idx = strcmp(wordsurp_encoding.SID, SID);
    lss(ctr) =  wordsurp_encoding.ls(find(idx, 1));

    % compare proportion of significant word-surp electrodes across speech
    % conditions
    if wordsurp_encoding.ls(find(idx, 1))==1 % Spanish speaker
        
        % first column is native language (Spanish), second column is
        % unfamiliar
        native_unfamiliar_prct(1, ctr) = ...
            sum(wordsurp_encoding.sp_wordsurp_pval(idx)<pthresh)/sum(idx);
        native_unfamiliar_prct(2, ctr) = ...
            sum(wordsurp_encoding.eng_wordsurp_pval(idx)<pthresh)/sum(idx);

    elseif wordsurp_encoding.ls(find(idx, 1))==2 % English speaker
        native_unfamiliar_prct(1, ctr) = ...
            sum(wordsurp_encoding.eng_wordsurp_pval(idx)<pthresh)/sum(idx);
        native_unfamiliar_prct(2, ctr) = ...
            sum(wordsurp_encoding.sp_wordsurp_pval(idx)<pthresh)/sum(idx);
    end
    % third column is number of electrodes
    native_unfamiliar_prct(3, ctr) = sum(idx);
    ctr=ctr+1;
end

figure;
rng(2);
sz = native_unfamiliar_prct(3, :)*2;
jitter = rand(length(lss), 1)*0.4-0.2;
cols = [1, 1, 1; 0.5, 0.5 0.5];
for ls  = unique(lss)' 
    lsidx = lss == ls;
    
    scatter(jitter(lsidx)+1, native_unfamiliar_prct(1, lsidx), sz(lsidx), ...
        cols(ls, :), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5]); hold on;
    scatter(jitter(lsidx)+2, native_unfamiliar_prct(2, lsidx), sz(lsidx), ...
        cols(ls, :), 'filled', 'MarkerEdgeColor', [0.5 0.5 0.5], 'HandleVisibility', ...
        'off');
end
legend({'Spanish', 'English'});

% formatting
xlim([0.5 2.5]);
ylim([0 1]);
xticks([1 2]);
xticklabels({'native', 'foreign'});
ylabel('significant word+surp \Delta R^2 (%)');
set(gca, 'FontSize', 13);
yticks([0 0.5 1]);
yticklabels({'0', '50', '100'});

% draw lines connecting all the scatters
plot([jitter+1 jitter+2]', [native_unfamiliar_prct(1, :); ...
    native_unfamiliar_prct(2, :)], 'Color', ...
        [0.7 0.7 0.7], 'HandleVisibility', 'off');

% lme model stats for unpaired difference
tbl = table(native_unfamiliar_prct(1, :)', ones(length(lss), 1), ...
    native_unfamiliar_prct(3, :)', lss, 'VariableNames', ...
    {'prct', 'native', 'ls', 'numel'});
% add the unfamiliar
tbl = [tbl; table(native_unfamiliar_prct(2, :)', zeros(length(lss), 1), ...
    native_unfamiliar_prct(3, :)', lss, 'VariableNames', ...
    {'prct', 'native', 'ls', 'numel'})];
lme = fitlme(tbl, 'prct ~ native + ls + (1|numel)');
%disp(lme);

% paired
[p, h] = signrank(native_unfamiliar_prct(1, :), ...
    native_unfamiliar_prct(2, :), 'Tail', 'right');

disp(p);
disp(['one-tailed paired sign rank ' num2str(p) ', observations = ' ...
    num2str(size(native_unfamiliar_prct, 2))]);

% lme model stats for paired difference
% tbl = table(native_unfamiliar_prct(1, :)'-native_unfamiliar_prct(2, :)', ...
%     lss, native_unfamiliar_prct(3, :)', 'VariableNames', ...
%     {'prct', 'ls', 'numel'});
% 
% lme = fitlme(tbl, 'prct ~ 1 + (1|ls) + (1|numel)');
% disp(lme);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

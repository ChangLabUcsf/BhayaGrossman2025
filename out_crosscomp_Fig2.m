% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;

% Note - EC202 has no STG coverage
SIDs = [sSIDs, eSIDs]; % , {'HS11', 'HS9', 'HS10'}
tps = 50:55;

% could add _wordFreqLog to the full model
% load in surprisal values
% using phonetic feature instead of splitting by vowel
modelnames_timit={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % 2 remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % 3 remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_engSurpNoOnsBin', ... % 4 remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... % 5 remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % 6 remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... % 7 remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % 8 full
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns', ... % 9 remove surprise and length
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordFreqLog', ... 
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin_wordFreqLog'}; % add in frequency
%'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... %remove word

modelnames_dimex={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_spSurpNoOnsBin', ... %remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_spSurpNoOnsBin', ... % remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin' ... % full
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordFreqLog', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin_wordFreqLog'}; % full v2

[wordsurp_encoding] = loadUniqueVarTbl(modelnames_timit, modelnames_dimex, SIDs);
wordsurp_details.featureOrd ...
    = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp', 'wordO', 'wordL', 'wordF'};
% the word feature contain the word onset and the word length, wordO is
% only onset, wordL is only length, wordF is wordFrequency
imgall = load_allimgdata;
wordsurp_encoding.hemi = cellfun(@(x) imgall.(x).hemi, ...
    wordsurp_encoding.SID, 'UniformOutput', false);
wordsurp_details.models_dimex = modelnames_dimex;
wordsurp_details.models_timit = modelnames_timit;

% load in p-values from permutation testing
sp_wordsurp_pval = nan(1, height(wordsurp_encoding));
en_wordsurp_pval = nan(1, height(wordsurp_encoding));
prefix = 'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_';
permodel = {[ prefix 'engSurpNoOnsBin_wordFreqLog'], ...
    [prefix 'spSurpNoOnsBin_wordFreqLog']};
for s = unique(wordsurp_encoding.SID)'
    SID = s{1};
    idx = strcmp(wordsurp_encoding.SID, SID);
    elidx = wordsurp_encoding.el(idx);

    for c = 1:2 %{'timit', 'dimex'}
        pvalpath=fullfile(datapath, 'permTest_wordSurp', SID); % c{1}, 
        cmod=dir(fullfile(pvalpath, '*_zX*_*mat')); 
        permidx = find(contains({cmod.name}, permodel{c}));
        permfname=cmod(permidx).name;

        % load pvalues from permutation testing
        pvals = load(fullfile(pvalpath, permfname), 'pval');
        pvals = pvals.pval;
        
        if c==1 % strcmp(c{1}, 'timit')
            en_wordsurp_pval(idx) = pvals(elidx);
        else
            sp_wordsurp_pval(idx) = pvals(elidx);
        end
    end
end
wordsurp_encoding.sp_wordsurp_pval = sp_wordsurp_pval';
wordsurp_encoding.eng_wordsurp_pval = en_wordsurp_pval';

rsq_thresh = 0.05;
idx = (wordsurp_encoding.eng_base_rsq>rsq_thresh& wordsurp_encoding.sp_base_rsq>rsq_thresh);
wordsurp_encoding = wordsurp_encoding(idx, :);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Alternate models
% Alternate models
% modelnames_timit={'phnfeatonset_maxDtL_wordOns_wordL_engSurpBin', ... % no onset feat
%     'onset_maxDtL_wordOns_wordL_engSurpBin', ... % no consonant  feat
%     'onset_phnfeatonset_wordOns_wordL_engSurpBin', ... % no peakrate feat
%     'onset_maxDtL_wordOns_wordL_engSurpBin', ... % no formant feat
%     'onset_phnfeatonset_maxDtL', ... % no word/surp feat (base)
%     'onset_phnfeatonset_maxDtL_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpBin', ... % full model feat
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_engSurpBin' % larger full model
%     };  
% 
% modelnames_dimex={'maxDtL_phnfeatonset_wordOns_wordL_spSurpBin', ... % no onset feat
%     'onset_maxDtL_wordOns_wordL_spSurpBin', ... % no consonant  feat
%     'onset_phnfeatonset_wordOns_wordL_spSurpBin', ... % no peakrate feat
%     'onset_maxDtL_wordOns_wordL_spSurpBin', ... % no formant feat
%     'onset_phnfeatonset_maxDtL', ... % no word/surp feat (base)
%     'onset_phnfeatonset_maxDtL_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpBin', ... % full model feat
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_spSurpBin' % larger full model
%     }; 

% modelnames_timit={'phnfeatConsOnset_maxDtL_formant_wordOns_wordL_engSurpBin', ... % no onset feat
%     'onset_maxDtL_formant_wordOns_wordL_engSurpBin', ... % no consonant  feat
%     'onset_phnfeatConsOnset_formant_wordOns_wordL_engSurpBin', ... % no peakrate feat
%     'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_engSurpBin', ... % no formant feat
%     'onset_phnfeatConsOnset_maxDtL_formant', ... % no word/surp feat (base)
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_engSurpBin', ... % full model feat
%     };  
% 
% modelnames_dimex={'phnfeatConsOnset_maxDtL_formant_wordOns_wordL_spSurpBin', ... % no onset feat
%     'onset_maxDtL_formant_wordOns_wordL_spSurpBin', ... % no consonant  feat
%     'onset_phnfeatConsOnset_formant_wordOns_wordL_spSurpBin', ... % no peakrate feat
%     'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_spSurpBin', ... % no formant feat
%     'onset_phnfeatConsOnset_maxDtL_formant', ... % no word/surp feat (base)
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_spSurpBin', ... % full model feat
%     }; 

%% A - Example annotated sentences 

% find a sentence you know the words for
wordlists = cellfun(@(x) join(x, ' '), {timit_details.sentdet.wordList});
x = find(strcmp(wordlists, 'the wagons were burning fiercely'));

sents = {42, x};

% Call makeSurprisal function for TIMIT and DIMEx details
[timit_details.sentdet] = makeSurprisal(timit_details.sentdet, 8, 'timit');
[dimex_details.sentdet] = makeSurprisal(dimex_details.sentdet, 8, 'dimex');

% Add in frequency and length
[timit_details.sentdet] = makeWordFreq(timit_details.sentdet, 8, 'timit');
[dimex_details.sentdet] = makeWordFreq(dimex_details.sentdet, 8, 'dimex');

% Assign corpus details to a cell array
corpus_details = {dimex_details, timit_details};

% Load mexbet2ipa mapping
load('mexbet2ipa.mat');

% Create a new figure
figure;

% Set the x-axis limits for the plots
xl = [0, 1.8];

% Iterate over the two corpora (DIMEx and TIMIT)
numcorp = 2;
for c = 1:numcorp
    details = corpus_details{c};

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
        xline(x(sent_info.wordOns>0), 'LineWidth', 2);
        xlim(xl);
        box off;
    end   
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% C - Native vs. Foreign: acoustic phonetic unique variance (box-plots)

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
fieldname = {'English', 'Spanish'};
conditionLabels = {'native', 'foreign'};
thresh = 0.001;

features = {'peakrate', 'formant', 'consonant'};
titles = {'PeakRate', 'Vowel Formants', 'Consonant'};
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
subplts = {1, 2, 3};
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

%% B - Native vs. Foreign: word / sequence surprisal unique variance (box-plots)

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
fieldname = {'English', 'Spanish'};
conditionLabels = {'native', 'foreign'};
thresh = 0.001;

features = {'wordO', 'wordF', 'wordL', 'surp'};
titles = { 'Boundary',  'Frequency', 'Length', 'Phoneme Surprisal'};
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
subplts = {1, 2, 3, 4};
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

%% D - MNI Brain with Native / Unfamiliar word / surprisal electrodes highlighted (brain)

% initialize design electrode structure
fields = {'sp_uv_all', 'eng_uv_all'}; 
uv_thresh = 0.001;

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
        ls = [1, 2]; % can only do 1-2
        
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
        binedges = [-1 uv_thresh:0.005:0.01 0.015:0.015:0.1];
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
        
        mni_lh = plotMNIElec(unique(wordsurp_encoding.SID), desel, 'lh', 0);
        
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

        mni_rh = plotMNIElec(unique(wordsurp_encoding.SID), desel, 'rh', 0);
        sgtitle(native);
        l = light;
        view(90, 0);
        set(l,'Style', 'infinite', 'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.85;

        % add a pie
        % axes('Position',[.6 .15 .3 .3])
        % p = pie([sum(mni_rh.cond>1), sum(mni_rh.cond==1)], [1 1]); 
        % p(1).FaceColor = [desel.cols(5, :)];
        % p(1).EdgeColor = 'none';
        % p(3).FaceColor = [0.6 0.6 0.6];
        % p(3).EdgeColor = 'none';
        % p(2).FontWeight = 'bold';
        % p(2).Color = 'w';
        % p(2).FontSize = 13;
        % p(4).FontWeight = 'bold';
        % p(4).Color = 'w';
        % p(4).FontSize = 13;

        % new pie figure that combines across hemispheres
        figure;
        wordsurp = sum(mni_lh.cond>1) + sum(mni_rh.cond>1);
        nonwordsurp = sum(mni_lh.cond==1) + sum(mni_rh.cond==1);
        p = pie([wordsurp, nonwordsurp], [1 1]);
        p(1).FaceColor = [desel.cols(5, :)];
        p(1).EdgeColor = 'none';
        p(3).FaceColor = [0.6 0.6 0.6];
        p(3).EdgeColor = 'none';
        p(2).FontWeight = 'bold';
        p(2).Color = 'w';
        p(2).FontSize = 13;
        p(4).FontWeight = 'bold';
        p(4).Color = 'w';
        p(4).FontSize = 13;
     end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% E - Native vs. Foreign: word and sequence surprisal unique variance (scatter)

% Specify the labels for analysis
labels = {'word+surp'};

% Create a figure for the scatter plots
figure('Position', [200, 100, 900, 900]);

% Set a UV threshold for filtering
uv_thresh = 0.001;
pval_thresh = 0.1;

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
        %neg = any([pvals_x,pvals_y]>pval_thresh, 2) | isnan(x) | isnan(y);
        neg = all([x,y]<uv_thresh, 2) | isnan(x) | isnan(y);
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
    title({['r= ' num2str(r) ','], ['p=' num2str(p, 4)]});
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
%% S6 - Correlating word and sequence surprise in native speech condition (scatter)

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
        scatter(x_all(orpos), y_all(orpos), 10, [0.5 0.5 0.5], 'filled','MarkerFaceAlpha', 0.8);
        ylabel([featnames{j} ' \Delta R^2']);
        xlabel([featnames{f} ' \Delta R^2']);

        % add a reference line
        xline(0, 'Color', 'k', 'LineWidth', 1.5);
        yline(0, 'Color', 'k', 'LineWidth', 1.5);

        % Perform permutation testing to compute correlation coefficient and p-value
        [r, p] = corr(x_all, y_all, 'Rows', 'complete', 'Type', 'Spearman');
        title(['r= ' num2str(r, 2) ', ' getSigStr(p, 2)], 'FontWeight', 'normal');
        ctr = ctr + 1;
    end
end


% xlim([-0.005 0.03]);
% xticks(0:0.01:0.03);
% ylim([-0.005 0.07]);
% yticks(0:0.02:0.06);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% S9 - Correlating unique variance of word / surp to unique variance of other features

% Extract unique variance of word+surp from the wordsurp table

% Top row of figure is English (TIMIT), bottom row is Spanish (DIMEx)
fields = {'eng_uv_all', 'sp_uv_all'};

labels = {{'word+surp', 'onset'}, {'word+surp', 'peakrate'}, ...
    {'word+surp', 'formant'}, {'word+surp', 'consonant'}};

for field = fields
    figure('Position', [200, 100, 900, 500]);

    % Select only native speakers of the language presented
    if strcmp(field, 'eng_uv_all')
        lsidx = ismember(wordsurp_encoding.ls, [2, 4]);
        corpus = 'English';
    else
        lsidx = ismember(wordsurp_encoding.ls, [1, 4]);
        corpus = 'Spanish';
    end

    ctr=1;
    for l = labels
        label = l{1};
        index = [];
        index(1) = find(ismember(wordsurp_details.featureOrd, label{1}));   
        index(2) = find(ismember(wordsurp_details.featureOrd, label(2)));   
        
        ax(ctr)=subplot(2, length(labels), ctr+length(labels));
       
        x = wordsurp_encoding.(field{1})(lsidx, index(2));
        y = wordsurp_encoding.(field{1})(lsidx, index(1));
        sid = cellfun(@(x) str2double(x(3:end)), ...
            wordsurp_encoding.SID(wordsurp_encoding.ls==ls));

        % Ensure positive unique variance for both features being compared
        pos = x>0 & y>0;
        uv_thresh = 0.001;

        [r(ctr), p(ctr)] = corr(x(pos), y(pos), 'Type', 'Spearman'); 
        scatter(x(pos), y(pos), 20, 'k', 'filled'); hold on;
%         scatter(x(~pos), y(~pos), 10, [0.5 0.5 0.5], 'filled');
        
        xticks(-0.05:0.05:0.15);
        yticks(-0.05:0.05:0.15);

        % Ensure limits capture 97% of scatter
        maxlim = prctile([x; y], 100);
        minlim = prctile([x; y], 3);  
        xlim([0 maxlim]);
        ylim([0 0.05]);
    
        % Add unique variance labels
        xlabel([label{2} ' \Delta R^2']);

        if ctr==1
            ylabel([label{1} ' \Delta R^2']);
        end
        %title({['r= ' num2str(r(ctr), 3) ',']},{['p=' num2str(p(ctr), 3)]});

        set(gca, 'FontSize', 13);
        ctr=ctr+1;
    end

    subplot(2, length(labels), 1:length(labels)); 
    bar(r, 'EdgeColor', 'none', 'FaceColor', [0.8 0.8 0.8]);
    % print the pvalue on top of every bar
    for i=1:length(r)
        text(i, r(i)+0.01, getSigStr(p(i), 2), 'HorizontalAlignment', ...
                'center', 'FontSize', 15);
    end
    
    xlim([0.25 4.75]);
    ylim([-0.5 0.5]);
    yticks([0, 0.4]);
    box off;
    ylabel('Spearman corr')
    set(gca, 'FontSize', 13)

    for i=(1:length(labels))+length(labels)
        subplot(2, length(labels), i);
        h = lsline;
        h.LineWidth = 2;
        h.Color = 'k';
        xline(0, 'Color', 'k', 'LineWidth', 1.5);
        yline(0, 'Color', 'k', 'LineWidth', 1.5);
    end
    sgtitle([corpus ' speech']);
end
%linkaxes(ax(1:3));   % phonetic feat

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% ----------------------------- UNUSED Panels ----------------------------
%% Comparing acoustic phonetic representations across word-boundary electrodes

uv_thresh = 0.0005;

% aggregate uvs from ONLY native speech conditions
native_uv = [wordsurp_encoding.sp_uv_all(wordsurp_encoding.ls==1, :); ...
    wordsurp_encoding.eng_uv_all(wordsurp_encoding.ls==2, :)];
native_rsq = [wordsurp_encoding.sp_full_rsq(wordsurp_encoding.ls==1, :); ...
    wordsurp_encoding.eng_full_rsq(wordsurp_encoding.ls==2, :)];

% get the index of word / surprisal electrodes in the native speech
% condition
wordidx = find(strcmp(wordsurp_details.featureOrd, 'word+surp'));
word_elecs = native_uv(:, wordidx) > uv_thresh;

% compare the unique variance of each feature family across the two groups
figure; 
feats = 1:4;
% light green, orange, teal, yello
cols = [114 171, 66;  44 170 170; 234, 177, 32;240 90 41;]/256;

for f = feats

    % get the index of the feature
    feat = wordsurp_details.featureOrd{f};
    featidx = find(strcmp(wordsurp_details.featureOrd, feat));
    pos_uv = native_uv(:, featidx) > uv_thresh; % positive uv for violin plot

    % subplot
    subplot(1, length(feats), f);
%     scatter(ones(sum(word_elecs), 1)-0.025+randn(sum(word_elecs), 1)*0.05, ...
%         native_uv(word_elecs, featidx), 5, 'k', 'filled'); hold on;
%     scatter(2*ones(sum(~word_elecs), 1)-0.025+randn(sum(~word_elecs), 1)*0.05, ...
%         native_uv(~word_elecs, featidx), 5, 'k', 'filled'); hold on;
    
    % violin plot for the uv of word elecs and non word elecs
% violin plot for the uv of word elecs and non word elecs
    violinplot([native_uv(word_elecs & pos_uv, featidx); native_uv(~word_elecs & pos_uv, featidx)], ...
        [ones(sum(word_elecs & pos_uv), 1); 2*ones(sum(~word_elecs & pos_uv), 1)], ...
        'MarkerSize', 5, 'ViolinColor', [cols(f, :); brighten(cols(f, :), 0.8)] , 'MedianColor', [0 0 0], ...
        'Bandwidth', 0.001, 'Width', 0.3);
    hold on;
    hold on;

    [h, p] = ttest2(native_uv(word_elecs & pos_uv, featidx), native_uv(~word_elecs & pos_uv, featidx));
    text(1, 0.08, getSigStr(p, 2));
    
    set(gca, 'YScale', 'log', 'FontSize', 13);
    ylim([0.000 0.1]);
    xlim([0.5 2.5]);
    xticks([1 2]);
    xticklabels({'word elecs', 'non-word elecs'});
    ylabel('\Delta R^2');
    %don't show legend
    legend('off');
    box off;
    title(feat);
end

% rsq value
% figure;
% violinplot([native_rsq(word_elecs); native_rsq(~word_elecs)], ...
%         [ones(sum(word_elecs), 1); 2*ones(sum(~word_elecs), 1)], ...
%         'MarkerSize', 5, 'ViolinColor', [cols(f, :); brighten(cols(f, :), 0.8)] , 'MedianColor', [0 0 0], ...
%         'Bandwidth', 0.1, 'Width', 0.3);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Native vs. Foreign: word + sequence surprisal unique variance (by subject box-plots)
thresh = 0.001;
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
            sum(wordsurp_encoding.sp_wordsurp_pval(idx)<thresh)/sum(idx);
        native_unfamiliar_prct(2, ctr) = ...
            sum(wordsurp_encoding.eng_wordsurp_pval(idx)<thresh)/sum(idx);

    elseif wordsurp_encoding.ls(find(idx, 1))==2 % English speaker
        native_unfamiliar_prct(1, ctr) = ...
            sum(wordsurp_encoding.eng_wordsurp_pval(idx)<thresh)/sum(idx);
        native_unfamiliar_prct(2, ctr) = ...
            sum(wordsurp_encoding.sp_wordsurp_pval(idx)<thresh)/sum(idx);
    end
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

% lme model
% tbl = table(native_unfamiliar_prct(1, :)', ones(length(lss), 1), ...
%     native_unfamiliar_prct(3, :)', lss, 'VariableNames', ...
%     {'prct', 'native', 'ls', 'numel'});
% % add the unfamiliar
% tbl = [tbl; table(native_unfamiliar_prct(2, :)', zeros(length(lss), 1), ...
%     native_unfamiliar_prct(3, :)', lss, 'VariableNames', ...
%     {'prct', 'native', 'ls', 'numel'})];

tbl = table(native_unfamiliar_prct(1, :)'-native_unfamiliar_prct(2, :)', ...
    lss, native_unfamiliar_prct(3, :)', 'VariableNames', ...
    {'prct', 'ls', 'numel'});

lme = fitlme(tbl, 'prct ~ 1 + (1|ls) + (1|numel)');
disp(lme);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Single subject electrode unique variance trends

SID = 'EC183';
ls = wordsurp_encoding.ls(find(strcmp(wordsurp_encoding.SID, SID), 1));
fieldname = {'English', 'Spanish'};
conditionLabels = {'native', 'foreign'};
thresh = 0.001;

% features = {'peakrate', 'formant', 'consonant'};
% titles = {'PeakRate', 'Vowel Formants', 'Consonant'};

features = {'wordO', 'wordF', 'wordL', 'surp'};
titles = { 'Boundary',  'Frequency', 'Length', 'Phoneme Surprisal'};
fields = {'eng_uv_all', 'sp_uv_all'}; 

ctr = 1;
native = cell(length(features), 1);
sid_native = cell(length(features), 1);
nonnative = cell(length(features), 1);
sid_nonnative = cell(length(features), 1);
el = cell(length(features), 1);
hemi = cell(length(features), 1);

figure; 
for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);
    
    for f = 1:2 % Aggregate over English and Spanish fields
        field = fields{f};    
        y = wordsurp_encoding.(field)(:, index);   
        sid_tmp = wordsurp_encoding.SID;   
        
        % Removing bilinguals so the comparison is more balanced
        if (f == 1 && ismember(ls, 2)) || (f == 2 && ismember(ls, 1))
            native = y(strcmp(SID, wordsurp_encoding.SID));  
        elseif (f == 1 && ismember(ls, 1)) || (f == 2 && ismember(ls, 2))
            nonnative = y(strcmp(SID, wordsurp_encoding.SID));  
        end
    end

    % scatter plot
    idx = native>thresh & nonnative>thresh;
    subplot(1, length(features), ctr);
    scatter(ones(sum(idx), 1), ...
        native(idx), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.8); hold on;
    scatter(2*ones(sum(idx), 1), ...
        nonnative(idx), 20, 'k', 'filled', 'MarkerFaceAlpha', 0.8);
    % plot a line for each electrode
    for i = find(idx)
        plot([1 2], [native(i), nonnative(i)], 'Color', [0.7 0.7 0.7]);
    end
    xlim([0.5 2.5]);
    [~, p] = ttest2(native, nonnative);
    title(p);
    xticks([1 2]);
    xticklabels(conditionLabels)

    ctr = ctr + 1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;
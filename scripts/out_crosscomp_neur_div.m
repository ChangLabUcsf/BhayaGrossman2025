%% Set up

% Ilina Bhaya-Grossman
% 03.23.2023
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('util'))

zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;

% Note - EC202 has no STG coverage
dSIDs = {'HS11', 'EC237', 'EC266', 'EC197', 'EC261', 'HS8', 'HS9', 'HS10'};
dLang_all = {'mandarin', 'russian', 'catalan', 'arabic', 'korean',  ...
    'mandarin', 'mandarin', 'mandarin'};
% english proficiency
prof_all = [4, 4, 3, 2, 0, 0, 0, 0];

% get color map
spec = spectral(8);
spec = flipud(spec([1:3 5:6], :));
cols = arrayfun(@(x) spec(x+1, :), prof_all, 'UniformOutput',false);
cols_all = cat(1, cols{:});

% changed colormap
spec = flipud([35, 100, 170; 61, 165, 217; 115, 191, 184; 254, 198, 1; 234, 115, 23;]./256);
cols = arrayfun(@(x) spec(x+1, :), prof_all, 'UniformOutput',false);
cols_all = cat(1, cols{:});

timit_details = load('out_sentence_details_timit_all_loudness.mat');
% asccd_details = load('stim_info/out_sentence_details_acssd_loudness.mat');
% tps = 50:55;

bef=50;
aft=50;

% loading in subject data
TDwrd = loadDwrdMulti('timit', bef, aft, dSIDs, timit_details);

% load language specific word data

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% --------------------------- English (TIMIT) ----------------------------
%% aggregate unique variance values from strf models

% using phonetic feature instead of splitting by vowel
% to do: add noOns surprisal feature

modelnames_timit={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... %remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... %remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_engSurpNoOnsBin', ... %remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... % remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin'}; % full
%'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... %remove word

modelnames = {modelnames_timit};

% 'consfeat', 'formant', 
lang = {'eng'};
corpus = {'timit'};

timit_elecs = load('diverse_lang/out_elecs_speechtypeftest_bychan_timit_DIV.mat');

% determine unique variance per feature and primary encoding
varnames = {'SID', 'el', 'ls', 'prof', ...
    [lang{1} '_base_rsq'], ...
    [lang{1} '_full_rsq'], ...
    [lang{1} '_uv_all'], ...
    [lang{1} '_full_beta'], ...
    };
wordsurp_encoding =  array2table(zeros(0, length(varnames)), 'VariableNames', varnames);

for s = dSIDs
    SID = s{1}; 
    ls = dLang_all(ismember(dSIDs, SID)); 
    prof = prof_all(ismember(dSIDs, SID));
    corpusStrf{1} = loadMultModelStrf(SID, modelnames{1}, corpus{1}, ...
        datapath, 1, 'v5');

    if ~any(cellfun(@(x) isempty(x), [corpusStrf{1}]))

        % find minimum test R
        minel = min(cellfun(@(x) length(x{3}.meanTestR), corpusStrf));
        els = timit_elecs.allidx.(SID);

        base = cell(2, 1);
        full = cell(2, 1);
        uvall = cell(2, 1);
        fullBeta = cell(2, 1);

        % preallocating cell space
        % full models
        noOns = cell(2, 1);
        noCons = cell(2, 1);
        noPeakr = cell(2, 1);
        noForm = cell(2, 1);
        noWordSurp = cell(2, 1);
        noWord = cell(2, 1);
        noSurp = cell(2, 1);

        % unique variance
        uvOns = cell(2, 1);
        uvCons = cell(2, 1);
        uvPeakr = cell(2, 1);
        uvForm = cell(2, 1);
        uvWordsurp = cell(2, 1);
        uvWord = cell(2, 1);
        uvSurp = cell(2, 1);

        for l = 1:1
            % full models
            base{l} = (corpusStrf{l}{5}.meanTestR.^2)';
            full{l} = (corpusStrf{l}{8}.meanTestR.^2)';
            fullBeta{l} = corpusStrf{l}{8}.strf;

            % single feature excluded models
            noOns{l} = (corpusStrf{l}{1}.meanTestR.^2)';
            noCons{l} = (corpusStrf{l}{2}.meanTestR.^2)';
            noPeakr{l} = (corpusStrf{l}{3}.meanTestR.^2)';
            noForm{l} = (corpusStrf{l}{4}.meanTestR.^2)';
            
            % single feature excluded models (surprisal)
            noWordSurp{l} = (corpusStrf{l}{5}.meanTestR.^2)';        
            noSurp{l} = (corpusStrf{l}{6}.meanTestR.^2)';
            noWord{l} = (corpusStrf{l}{7}.meanTestR.^2)';
            
            % unique variances
            uvOns{l} = full{l}-noOns{l};
            uvPeakr{l} = full{l}-noPeakr{l};
            uvForm{l} = full{l}-noForm{l};
            uvCons{l} = full{l}-noCons{l};
            uvWordsurp{l} = full{l}-noWordSurp{l};

            % this one is calculated without the surprisal feature
            % uvWord{l} = noSurp{l}-noWordSurp{l};

            % now calculated with the surprisal feature
            uvWord{l} = noSurp{l}-noWordSurp{l};
            uvSurp{l} = noWord{l}-noWordSurp{l};
        end

        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        profs = repmat(prof, length(els), 1);

        for l = 1:1
            betaCell{l} = squeeze(mat2cell(fullBeta{l}{1}(:, :, els), ...
                size(fullBeta{l}{1}, 1), size(fullBeta{l}{1}, 2), ones(length(els),1)));
        end


        tmp = table(sids, els, lss, profs, ...
            base{1}(els), full{1}(els), ...
            [uvOns{1}(els), uvPeakr{1}(els), uvForm{1}(els), uvCons{1}(els), ...
                uvWordsurp{1}(els), uvWord{1}(els), uvSurp{1}(els)], betaCell{1}, ...
            'VariableNames', varnames);

        wordsurp_encoding = [wordsurp_encoding; tmp];
    else
        warning(['Missing subject ' SID]);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *wrd;

%% run model uv comparison
%% vis: uv comparison across subject

thresh = 0.001;

featureOrd = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
features = {'word+surp', 'formant', 'consonant', 'peakrate'}; %'peakrate', , 
% titles = {'peakrate', 'phonetic features', 'word+surp'};

figure;
ctr=1;
fields = {'eng_uv_all'}; 
for feat = features
    index = ismember(featureOrd, feat);

    % show individual electrodes per subject
    field = fields{1};    
    y = wordsurp_encoding.(field)(:, index);   
    meethresh = y>thresh;

    uvmed=nan(length(dSIDs), 1);
    subplot(2, length(features), ctr)
    for i = 1:length(dSIDs)
        SID = dSIDs{i};
        
%         boxchart(repmat(i, sum(sidx), 1), y(sidx), ...
%             'MarkerStyle','none', ...
%             'BoxFaceColor', cols_all(i, :)); hold on;
        sidx = strcmp(wordsurp_encoding.SID, SID);
        h = bar(i, sum(sidx & meethresh)); 
        hold on; box off;
        h.FaceColor = cols_all(i, :);
        h.EdgeColor = "none";
        set(gca, 'FontSize', 15, 'XTickLabelRotation', 60)

        uvmed(i) = sum(sidx & meethresh);

%         % plot those electrodes that do not meet threshold
%         scatter(repmat(i, sum(sidx & ~meethresh), 1) ...
%             - 0.1 + 0.2 *rand(sum(sidx & ~meethresh), 1), ...
%             y(sidx & ~meethresh), 5, [0.7 0.7 0.7], 'filled', ...
%             'MarkerFaceAlpha', 0.6);
%         hold on;
% 
%          % plot those electrodes that do meet threshold
%         scatter(repmat(i, sum(sidx & meethresh), 1) ...
%             - 0.1 + 0.2*rand(sum(sidx & meethresh), 1), ...
%             y(sidx & meethresh), 10, cols_all(i, :), 'filled', ...
%             'MarkerFaceAlpha', 0.6);
% 
%         % find median for uvs greater than zero
%         uvmed(i) = median(y(sidx & meethresh));
%         line([i-0.25 i+0.25], [uvmed(i) uvmed(i)], 'LineWidth', ...
%             3, 'Color', cols_all(i, :));   
%         ylim([-0.01 0.025]);
    end

    % Plot a best fit line through subject medians
%     coefs = polyfit(prof_all, uvmed, 1);
    
    % Use polyval to calcualte function values of that fit
%     plot(1:0.1:length(dSIDs), fliplr(polyval(coefs, 1:0.1:length(dSIDs))),'-', ...
%         'LineWidth',2.5, 'Color',[0.6 0.6 0.6]);
    clear coefs

    title(feat);
    xticks(1:8);
    xticklabels(dLang_all);
    ctr=ctr+1;

    % show proportion of electrodes per subject
    ax = subplot(2, length(features), ctr+length(features)-1);

    perc = nan(1, length(dSIDs));
    count = nan(1, length(dSIDs));
    total = perc;
    for i = 1:length(dSIDs)
        SID = dSIDs{i};
        sidx = strcmp(wordsurp_encoding.SID, SID);

        % percentage of electrodes with positive unique values
        perc(i) = sum(meethresh & sidx) / sum(sidx);
        count(i) = sum(meethresh & sidx);
        total(i) = sum(sidx);
    end
    scatter(1:length(dSIDs), perc, total*5, cols_all, ...
            'filled', 'MarkerFaceAlpha', 0.5); hold on;
    lsline;
    title(feat);
    xticks(1:8);
    xlim([0 9])
    xticklabels(dLang_all);

    tbl = wordsurp_encoding(:, {'SID', 'prof', 'el'});
    tbl.uv = y;
    lme = fitlme(tbl,'uv~1+prof+(1|SID:el)');
    disp(lme)
    disp([feat{1} ': ' num2str(lme.Coefficients.pValue(2))]);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: uv comparison within subject
thresh = 0.00;

featureOrd = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
features = {'peakrate', 'formant', 'consonant', 'word+surp'}; %'peakrate', , 
% titles = {'peakrate', 'phonetic features', 'word+surp'};

fields = {'eng_uv_all'}; 
cols = lines(length(features));
for i = 1:length(dSIDs)
    SID = dSIDs{i};
    sidx = strcmp(wordsurp_encoding.SID, SID);

    figure;
    j = 1;
    for feat = features
        index = ismember(featureOrd, feat);

        % show individual electrodes per subject
        field = fields{1};    
        y = wordsurp_encoding.(field)(:, index);   
        meethresh = y>thresh;

        boxchart(repmat(j, sum(sidx & meethresh), 1), y(sidx & meethresh), ...
            'MarkerStyle','none', 'BoxFaceColor',cols(j, :)); hold on;

        scatter(repmat(j, sum(sidx & meethresh), 1) - 0.1 + ...
            0.2*rand(sum(sidx&meethresh), 1), ...
            y(sidx & meethresh), 35, cols(j, :), 'filled');

        j = j+1;
    end
    set(gca, 'YScale', 'log');
    title([upper(dLang_all{i}) ', english proficiency: ' ...
        num2str(prof_all(i)) '/5']);
    xticks(1:length(features));
    xticklabels(features);  
    ylabel('feature \Delta R^2');
    xlabel('speech feature');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: native brains

% initialize design electrode structure
fieldnames = {'English'};
fields = {'eng_uv_all'}; 
thresh = 0.001;

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
feat = 'formant';
% uv feature order
featureOrd = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word'};
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant',

fig = figure();
for lang = 1:1
    for f = 1:length(feats)
        feat = feats{f};
        index = find(ismember(featureOrd, feat));
    
        switch feat
            case 'peakrate'
                featcol = [0.3 0.8 0.9];
            case 'formant'
                featcol = [0.4 0.7 0.3];
            case 'consonant'
                featcol = [0.6 0.0 0.7];
            case 'word+surp'
                featcol = [0.7 0.1 0.2];
            otherwise
                featcol = [0 0 0];
        end   
    
        % make desel structure
        desel=struct();
        desel.conds = 1:7;
        ls = [1, 2, 4]; % can only do 1-2
        
        % size and color
        desel.sz = [35; 35*ones(length(desel.conds), 1)]; %1:20:200; %ones(1, 10)*0.00000001; %1
        
        % split up peak rate and phonetic features again for MNI plotting wordsurp_encoding.ls(x)
        desel.labels = [];
        desel.yval = arrayfun(@(x) wordsurp_encoding.(fields{lang})(x, index), ...
            1:height(wordsurp_encoding));
        
        % discretize values
        % yvals = sort(desel.yval(desel.yval>0 & ismember(wordsurp_encoding.ls,ls)'));
        % binedges = yvals(1:ceil(length(yvals)/8):length(yvals));
        % [~, binedges] = discretize(desel.yval(desel.yval>0 & ismember(wordsurp_encoding.ls,ls)'), ...
        %     length(desel.conds)-1);
        % manual non-linear edges
        binedges = [0:0.005:0.01 0.015:0.015:0.045];
        binedges = [-1, binedges];
        
        for s=unique(wordsurp_encoding.SID)'
            SID = s{1};
            idx = strcmp(wordsurp_encoding.SID, SID);
            desel.(SID).elid = wordsurp_encoding.el(idx);
            desel.(SID).condition = discretize(desel.yval(idx), ...
                binedges);
        end
        
        desel.cols = [1 1 1; [linspace(0, featcol(1), length(binedges)); ...
            linspace(1, featcol(2), length(binedges)); ...
            linspace(1, featcol(3), length(binedges))]'];
        if strcmp(feat, 'word+surp')
                
            gns = flipud(fpurple(length(binedges)-2));
            desel.cols = [0 0 0 ; gns];
    
            % colorbar
            figure;
            colormap(desel.cols);
            colorbar;            
        end
        
%         lsid = find(ismember(wordsurp_encoding.ls, ls));
        SIDs = dSIDs;
    
        [native_plot] = plotNativeElec(SIDs, desel, 1);
        sgtitle(lang)
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;


%% run word boundary decoding

load('diverse_lang/out_elecs_speechtypeftest_bychan_timit_DIV.mat', 'allidx');
Swrd = TDwrd;
corpus = 'timit';

% single window before onset
startp = 31;
timing = 57;
timelabel = '600ms';

tps = startp(1):startp(end)+max(timing);
nreps = 20;

tic
varnames = {'SID', 'elecs', 'ls', 'acc', 'weights', 'auc'}; 

wordBoundary_logisticwrapper(Swrd, dSIDs, wordsurp_encoding, corpus, ...
    startp, timing, timelabel, 'diverse', 500, nreps);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;


%% vis: word boundary decoding

corpus = 'TIMIT';
timelabel = '600ms';

filename = [corpus '_word_decode_diverse_' timelabel '_bysubj.mat']; % dimex filename
load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

for i = 1:length(dSIDs)
    nreps = length(decode_details.tbl.auc(i, :));
%     boxchart(i*ones(nreps,1 ), decode_details.tbl.auc(i, :), ...
%         'BoxFaceColor', cols_all(i, :), 'MarkerColor','k', ...
%         'LineWidth',2, 'JitterOutliers','on', 'MarkerStyle','.'); hold on;
    scatter(i*ones(nreps,1)-0.1+rand(nreps, 1)*0.2, decode_details.tbl.auc(i, :), ...
        15,cols_all(i, :), "filled", 'HandleVisibility','off', ...
        'MarkerFaceAlpha', 0.6); hold on;
    aucmean = mean(decode_details.tbl.auc(i, :));
    line([i-0.25 i+0.25], [aucmean aucmean], 'LineWidth', 3, 'Color', cols_all(i, :));
end

xticks(1:length(dSIDs));
yline(0.5);
xticklabels(decode_details.tbl.ls);
l = legend({'', '4', '3', '2', '0'});

title(l, 'English Proficiency');

set(gca, 'FontSize', 13);
ylabel('AUC');
ylim([0.3 0.9]);
yticks([0.3 0.9]);
xlim([0 9]);
xlabel('Native Language');

% run stats
auc = [decode_details.tbl.auc([1:4, 6:8], :)];
prof = repmat(prof_all([1:4, 6:8]), 20, 1)';

lme_tbl = table();
lme_tbl.auc = auc(:);
lme_tbl.prof = prof(:);
lme = fitlme(lme_tbl,'auc~1+prof');
disp(lme)

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% --------------------------- Other languages ----------------------------
%% korean

corpus = 'korean';
SID = 'EC261';
[out] = loadOutFile(datapath, corpus, SID, 2);

% for i = 1:size(out, 2)
%     out(i).soundOns = [0, out(i).duration-sum(out(i).befaft)];
% end
% [out,fn] = out_mk_envelope(out, 'loudness', [0.5, 0.5]);

% plot envelope at word and syllable boundary
a.sentdet = out;
KDwrd = loadDwrd_diverse('korean', 50, 50, {SID}, a);
idx = KDwrd.sentOns<1 & (KDwrd.syll>1 | isnan(KDwrd.syll));
KDwrd = KDwrd(idx, :);

ywrd = squeeze(cat(3, KDwrd.env{logical(KDwrd.wordOns)}))';
ysyll = squeeze(cat(3, KDwrd.env{~KDwrd.wordOns}))';
plotEnv(ywrd, ysyll, [0.5 0.5]);
xlim([-0.3 0.3])

% plot Korean HFA peak to English HFA peak
outs = cell(1, 2);
outs{1} = out;
[outs{2}] = loadOutFile(datapath, 'timit', SID, 0);

plotSentenceResp(outs, {corpus, 'timit'}, SID)

%% korean word onset decoding

% single window before onset
startp = 31;
timing = 57;
timelabel = '600ms';
Swrd =  KDwrd;
index = find(strcmp(dLang_all, 'korean'));

corpus = 'korean';

% subset to only the monolinguals
mono_encoding = wordsurp_encoding(strcmp(wordsurp_encoding.ls, corpus), :);

uSIDs = [unique(mono_encoding.SID)];
[decode_details] = wordBoundary_logisticwrapper(Swrd, uSIDs, mono_encoding, corpus, ...
    startp, timing, timelabel);

% Plot korean word boundary decoding results
figure;
nreps = length(decode_details.tbl.auc(1, :));
scatter(1*ones(nreps,1)-0.1+rand(nreps, 1)*0.2, decode_details.tbl.auc(1, :), ...
        15,cols_all(index, :), "filled", 'HandleVisibility','off', ...
        'MarkerFaceAlpha', 0.6); hold on;
aucmean = mean(decode_details.tbl.auc(1, :));
line([1-0.25 1+0.25], [aucmean aucmean], 'LineWidth', 3, 'Color', cols_all(index, :));

% Plot against english decoding results
corpus = 'TIMIT';
timelabel = '600ms';

filename = [corpus '_word_decode_diverse_' timelabel '_bysubj.mat']; % dimex filename
load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

nreps = length(decode_details.tbl.auc(index, :));
scatter(2*ones(nreps,1)-0.1+rand(nreps, 1)*0.2, decode_details.tbl.auc(index, :), ...
        15, cols_all(index, :), "filled", 'HandleVisibility','off', ...
        'MarkerFaceAlpha', 0.6); hold on;
aucmean = mean(decode_details.tbl.auc(index, :));
line([2-0.25 2+0.25], [aucmean aucmean], 'LineWidth', 3, 'Color', cols_all(index, :));

xticks([1, 2]);
xticklabels({'Korean', 'English'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% russian

corpus = 'russian';
SID = 'EC237';

[outs{1}] = loadOutFile(datapath, corpus, SID, 2);
[outs{2}] = loadOutFile(datapath, 'timit', SID, 1);

% % fix this
% for i = 1:size(outs{1}, 2)
%     outs{1}(i).soundOns = [0, outs{1}(i).duration];
% end
% [outs{1},~] = out_mk_envelope(outs{1}, 'loudness', [0.2, 0.2]);

% plot Russian HFA peak to English HFA peak
plotSentenceResp(outs, {corpus, 'timit'}, SID);


%% russian word onset decoding
%% arabic

corpus = 'arabic';
SID = 'EC197';

[outs{1}] = loadOutFile(datapath, corpus, SID, 55);
[outs{2}] = loadOutFile(datapath, 'timit', SID, 2);

% plot Arabic HFA peak to English HFA peak
plotSentenceResp(outs, {corpus, 'timit'}, SID)

% for i = 1:size(out, 2)
%     out{1}(i).soundOns = [0, out(i).duration];
% end
% [out{1},~] = out_mk_envelope(out{1}, 'loudness', [0.2, 0.2]);


%% catalan

corpus = 'catalan';
SID = 'EC266';

[outs{1}] = loadOutFile(datapath, corpus, SID, 26);
[outs{2}] = loadOutFile(datapath, 'timit', SID, 16);

% for i = 1:size(outs{1}, 2)
%     outs{1}(i).soundOns = [0, outs{1}(i).duration];
% end
% [outs{1},~] = out_mk_envelope(outs{1}, 'loudness', [0.2, 0.2]);

plotSentenceResp(outs, {corpus, 'timit'}, SID);


%% catalan word onset decoding

%% mandarin

corpus = 'asccd';
SID = 'HS8';

[outs{1}] = loadOutFile(datapath, corpus, SID, 1, 'noCAR_');
[outs{2}] = loadOutFile(datapath, 'timit', SID, 10);

% for i = 1:size(outs{1}, 2)
%     outs{1}(i).soundOns = [0, outs{1}(i).duration];
% end
% [outs{1},~] = out_mk_envelope(outs{1}, 'loudness', [0.2, 0.2]);

plotSentenceResp(outs, {corpus, 'timit'}, SID);

%% mandarin 2
corpus = 'asccd';
SID = 'HS9';

[outs{1}] = loadOutFile(datapath, corpus, SID, 1, 'noCAR_');
[outs{2}] = loadOutFile(datapath, 'timit', SID, 10, 'noCAR_');

% for i = 1:size(outs{1}, 2)
%     outs{1}(i).soundOns = [0, outs{1}(i).duration];
% end
% [outs{1},~] = out_mk_envelope(outs{1}, 'loudness', [0.2, 0.2]);

plotSentenceResp(outs, {corpus, 'timit'}, SID);

%% mandarin 3
corpus = 'asccd';
SID = 'HS10';

[outs{1}] = loadOutFile(datapath, corpus, SID, 1, 'noCAR_');
[outs{2}] = loadOutFile(datapath, 'timit', SID, 10, 'noCAR_');

% for i = 1:size(outs{1}, 2)
%     outs{1}(i).soundOns = [0, outs{1}(i).duration];
% end
% [outs{1},~] = out_mk_envelope(outs{1}, 'loudness', [0.2, 0.2]);

plotSentenceResp(outs, {corpus, 'timit'}, SID);


%% mandarin word boundary decoding
%% vis: word boundary decoding
%% ------------------------------ functions -------------------------------

function [out] = loadOutFile(datapath, corpus, SID, block, CARflag)

    if nargin<5, CARflag=''; end

    filedir= fullfile([datapath SID '/' corpus '/block_z/']);
    if nargin>3
        filename = [SID ...
        '_B' num2str(block) '_HilbAA_70to150_8band_' CARflag 'all_0_mel_' ...
        upper(corpus) '_zflag_global_out_resp_log'];
    else
        sd = dir(filedir);
        filename = sd(3).name;
    end

    load([filedir filename], 'out');
end


function plotSentenceResp(outs, corps, SID)

    dur = nan(1, length(corps));
    befs = cell(1, length(corps));
    dataf = nan(1, length(corps));
    for i = 1:length(corps)
        % find minimum number of time points for both corpora
        dur(i) = min(cellfun(@(x) size(x, 2), {outs{1, i}.resp}));
        befs(i) = {outs{1, i}(1).befaft};
        dataf(i) = outs{1, i}.dataf;             
    end
    befs = arrayfun(@(x) befs{x}(1)*dataf(x), 1:length(corps));
    mintp = min(arrayfun(@(x) dur(x)-befs(x), 1:length(corps)));

    % plot imagesc
    figure;
    maxresp =cell(1, length(corps));
    for i = 1:length(corps)
        
        resp = cellfun(@(x) x(:, befs(i):befs(i)+mintp), ...
            {outs{i}.resp}, 'UniformOutput', false);
        resp = cat(3, resp{:});

        ax = subplot(1, 3, i);
        imagesc(mean(resp, 3));
        title(upper(corps{i}));
        colormap(ax, flipud(piyg));
        clim([-2 2]);
        box off;
        colorbar();

        % maximum HFA response
        tmp = cellfun(@(x) max(smoothdata(x, 2, "gaussian"), [], [2, 3]), ...
            {outs{i}.resp}, 'UniformOutput', false);
        maxresp{i} = cat(2, tmp{:});        
    end
    
    % edit this so that we can see both brains like in Fig. 1
    if length(corps)>1
        % plot HFA difference
        ax = subplot(1, 3, 3);
        peakresp = [rescale(mean(maxresp{1}, 2)), rescale(mean(maxresp{2}, 2))];
        scatter(peakresp(:, 1), peakresp(:, 2), 25, diff(peakresp'), 'filled');
    
        xlabel([upper(corps{1}) ' HFA peak (rescale)']);
        ylabel([upper(corps{2}) ' HFA peak (rescale)']);
        [r, pval]= corr(peakresp(:, 1), peakresp(:, 2));
        title({['r=' num2str(r)] ,['p=' num2str(pval)]});
        
        clim([-0.5 0.5]);
        colormap(ax, flipud(spectral));
        refline(1, 0);
    
        % plot native brain with difference map
        [imgall] = load_allimgdata;
        if isfield(imgall, SID)
            img_native = imgall.(SID).img_native;
        
            figure;
            ax1 = axes;
            hold(ax1, 'all');   
            ctmr_gauss_plot(img_native.cortex,[0 0 0], 0, imgall.(SID).hemi);
            alpha 0.6;
        
            % plot electrodes on native brain
            % Hide the top axes
             ax2 = ax1;
        
            % use the minimum, either from out struct or from elecmatrix
            elec = 1:min(size(peakresp, 1), size(img_native.elecmatrix, 1));
            scatter3(ax2, img_native.elecmatrix(elec, 1), ...
                img_native.elecmatrix(elec, 2), ...
                img_native.elecmatrix(elec, 3), 30, diff(peakresp(elec, :)'), 'o', ...
                'filled', 'MarkerFaceAlpha', 0.85, 'MarkerEdgeColor', [0 0 0]);  hold on;
    
            ax2.Visible = 'off';
            ax2.XTick = [];
            ax2.YTick = [];   
            caxis(ax2, [-0.5 0.5]);
            colormap(ax2, flipud(spectral));
    
             %%Link them together
            linkaxes([ax1,ax2], 'xyz');   
            linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', ...
                'CameraTarget', 'XLim', 'YLim', 'ZLim'});
            % this will error out if 3rd child is not brain patch
            ax2.Children(3).FaceColor = [0.8 0.8 0.8];
    
            l = light;
            if strcmp(imgall.(SID).hemi, 'lh')
                view(270, 0);   
                set(l,'Style', 'infinite', 'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
            elseif strcmp(imgall.(SID).hemi, 'rh')
                view(90, 0);
                set(l,'Style', 'infinite', 'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
            end   

            %% alternative plotting

            % plot native brain
            cm = [0 0 1; 1 0 0];   
            load('diverse_lang/out_elecs_speechtypeftest_bychan_timit_DIV', ...
                'allidx')
            % first in native language, second is english
            for lang = 1:2
                figure;
                ax1 = axes;
                hold(ax1, 'all');   
                ctmr_gauss_plot(img_native.cortex,[0 0 0], 0, imgall.(SID).hemi);
                alpha 0.6;
            
                % plot electrodes on native brain
                % Hide the top axes
                ax2 = ax1;
                
                % use the minimum, either from out struct or from elecmatrix
                elec = allidx.(SID);
                scatter3(ax2, img_native.elecmatrix(elec, 1), ...
                    img_native.elecmatrix(elec, 2), ...
                    img_native.elecmatrix(elec, 3), ...
                    round(peakresp(elec, lang)*70)+10, cm(lang, :), ...
                    'o', 'filled', 'MarkerFaceAlpha', 0.65, ...
                    'MarkerEdgeColor', [1 1 1]);  hold on;

                ax2.Visible = 'off';
                ax2.XTick = [];
                ax2.YTick = [];   
                caxis(ax2, [-0.5 0.5]);
        
                 %%Link them together
                linkaxes([ax1,ax2], 'xyz');   
                linkprop([ax1, ax2],{'CameraUpVector', 'CameraPosition', ...
                    'CameraTarget', 'XLim', 'YLim', 'ZLim'});
                % this will error out if 3rd child is not brain patch
                ax2.Children(3).FaceColor = [0.8 0.8 0.8];

                % plot electrodes that did not meet threshold
                elec_all = 1:256;
                elec = elec_all(~ismember(elec_all, allidx.(SID)));
                scatter3(ax2, img_native.elecmatrix(elec, 1), ...
                    img_native.elecmatrix(elec, 2), ...
                    img_native.elecmatrix(elec, 3), ...
                    10, [0 0 0], ...
                    'o', 'filled', 'MarkerFaceAlpha', 0.65, ...
                    'MarkerEdgeColor', [1 1 1]);  hold on;
        
        
                l = light;
                if strcmp(imgall.(SID).hemi, 'lh')
                    view(270, 0);   
                    set(l,'Style', 'infinite', 'Position',[-1 0 0],'Color', ...
                        [0.8 0.8 0.8]);
                elseif strcmp(imgall.(SID).hemi, 'rh')
                    view(90, 0);
                    set(l,'Style', 'infinite', 'Position',[1 0 0],'Color', ...
                        [0.8 0.8 0.8]);
                end   
            end
            
        end
    end
end
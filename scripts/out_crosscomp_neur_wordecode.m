%% Set up

% Ilina Bhaya-Grossman
% 01.08.2022
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('util'))

zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
SIDs = [sSIDs eSIDs bSIDs];

% corpus details
timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
% asccd_details = load('stim_info/out_sentence_details_acssd_loudness.mat');
% tps = 50:55;

% selected electrodes
timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

% before and after word time points
bef=50;
aft=50;
%
% loading in word-level subject data % eSIDs, sSIDs
%TDwrd = loadDwrdMulti('timit', bef, aft, [eSIDs sSIDs], timit_details);
Dwrd = loadDwrdMulti('dimex',  bef, aft, [eSIDs sSIDs], dimex_details);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% -------------- NEURAL CLASSIFICATION OF WORD BOUNDARY ------------------
%% Load in unique variance of word features over spectrogram

% alternate TRF models
% with pitch change feature
% modelnames={'onset_phnfeatonset_maxDtL', ...
%     'onset_phnfeatonset_maxDtL_engSurp', ...
%     'onset_phnfeatonset_maxDtL_engSurpBin', ...
%     'onset_phnfeatonset_maxDtL_engSurp_engSurpBin'};   
% modelnames={'aud', ...
%     'aud_engSurp', ...
%     'aud_engSurpBin', ...
%     'aud_engSurp_engSurpBin'};    
% modelnames={'onset_aud_phnfeatonset_maxDtL', ...
%     'onset_aud_phnfeatonset_maxDtL_engSurp', ...
%     'onset_aud_phnfeatonset_maxDtL_engSurpBin'};    
% , ...'onset_aud_phnfeatonset_maxDtL_engSurp_engSurpBin'

% load in TRF models   
varnames = {'SID', 'el', 'ls', 'base_rsq', 'wordOns_rsq', 'wordLn_rsq', ...
    'uv_wordOns', 'uv_word', 'uv_peakRate'}; % , 'uv_peakRate'
word_encoding = array2table(zeros(0, 9), 'VariableNames', varnames);
modelnames={'onset_maxDtL_aud', ...    
    'onset_aud_maxDtL_wordOns', ...
    'onset_aud_maxDtL_wordOns_wordL', ...
    'onset_aud'};  %     'onset_aud'
corpus = {'timit', 'dimex'};

% determine unique variance per feature and primary encoding
for s = SIDs
    SID = s{1}; 
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs}));

    corpusStrf = cell(1, 2);
    for l = 1:2
        corpusStrf{l} = loadMultModelStrf(SID, modelnames(1:end), corpus{l}, ...
            datapath, 1);
    end

    % timit has more subject coverage (includes mandarin patients)
    if ~any(cellfun(@(x) isempty(x), corpusStrf{1}))
        numel = length(corpusStrf{1}{1}.meanTestR);
        basersq = ones(numel, 2)*-1;
        wrdrsq = ones(numel, 2)*-1;
        wrdlenrsq = ones(numel, 2)*-1;
        audrsq = ones(numel, 2)*-1;

        % models without pitch
        basersq(:, 1)   = (corpusStrf{1}{1}.meanTestR.^2)';
        wrdrsq(:, 1)    = (corpusStrf{1}{2}.meanTestR.^2)';
        wrdlenrsq(:, 1) = (corpusStrf{1}{3}.meanTestR.^2)';
        audrsq(:, 1)    = (corpusStrf{1}{4}.meanTestR.^2)';

        % Extract mean test R-squared values for each model
        if ~any(cellfun(@(x) isempty(x), corpusStrf{2}))
            numel_dimex = length(corpusStrf{2}{1}.meanTestR);
            basersq(1:numel_dimex, 2)   = (corpusStrf{2}{1}.meanTestR.^2)';
            wrdrsq(1:numel_dimex, 2)    = (corpusStrf{2}{2}.meanTestR.^2)';
            wrdlenrsq(1:numel_dimex, 2) = (corpusStrf{2}{3}.meanTestR.^2)';
            audrsq(1:numel_dimex, 2) = (corpusStrf{2}{4}.meanTestR.^2)';
        end
    
        % No longer finding electrodes that meet the R-squared threshold for the full model
        % els = find(sum(basersq > 0.1, 2)>0);

        % Include all speech responsive electrodes instead
        if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            els = union(timit_elecs.allidx.(SID), dimex_elecs.allidx.(SID));
        elseif isfield(dimex_elecs.allidx, SID)
            els = dimex_elecs.allidx.(SID);
        elseif isfield(timit_elecs.allidx, SID)
            els = timit_elecs.allidx.(SID);
        end

        % Calculate unique variances for each feature and primary encoding
        uvwrdons= wrdrsq - basersq;
        uvwrd = wrdlenrsq - basersq;
        uvpeakrate = basersq - audrsq;
        
        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        
        tmp = table(sids, els, lss, basersq(els, :), wrdrsq(els, :), ...
            wrdlenrsq(els, :), uvwrdons(els, :), uvwrd(els, :),  ...
            uvpeakrate(els, :), 'VariableNames', ...
            varnames);
        word_encoding = [word_encoding; tmp];
    else
        warning(['Missing subject ' SID])
    end
end

clearvars -except *all subj *vow* *cons* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh Dcons *wrd* *elecs;

%% --------------------------- all subjects -------------------------------
%% Load in unique variance of phoneme surprisal + word boundary over feature model

modelnames_dimex={'onset_maxDtL_wordOns_wordL_spSurpNoOnsBin', ... % remove phonetic feature
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... % remove word
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... % remove peakrate
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ...% full model
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; % remove surprise and word

modelnames_timit = {'onset_maxDtL_wordOns_wordL_engSurpNoOnsBin', ... % remove phonetic feature
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... % remove word
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % remove peakrate
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ...% full model
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; % remove surprise and word

% modelnames_dimex={'onset_maxDtL_wordOns_wordL_spSurpBin', ... % remove phonetic feature
%     'onset_phnfeatonset_maxDtL_spSurpBin', ... % remove word
%     'onset_phnfeatonset_maxDtL_wordOns_wordL', ... % remove surprise
%     'onset_phnfeatonset_wordOns_wordL_spSurpBin', ... % remove peakrate
%     'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpBin', ...% full model
%     'onset_phnfeatonset_maxDtL'}; % remove surprise and word
% 
% modelnames_timit = {'onset_maxDtL_wordOns_wordL_engSurpBin', ... % remove phonetic feature
%     'onset_phnfeatonset_maxDtL_engSurpBin', ... % remove word
%     'onset_phnfeatonset_maxDtL_wordOns_wordL', ... % remove surprise
%     'onset_phnfeatonset_wordOns_wordL_engSurpBin', ... % remove peakrate
%     'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpBin', ...
%     'onset_phnfeatonset_maxDtL'}; % remove surprise and word};

modelnames = {modelnames_timit, modelnames_dimex};
corpus = {'timit', 'dimex'};

% Determine unique variance per feature and primary encoding
for s = unique(word_encoding.SID)'
    SID = s{1}; 
    
    for l = 1:2
        % Load STRF models for the specific SID, model names, and corpus
        corpusStrf = loadMultModelStrf(SID, modelnames{l}(1:end), corpus{l}, ...
            datapath, 1);
        
        if ~any(cellfun(@(x) isempty(x), corpusStrf))
            % Extract mean test R-squared values from the loaded STRF models
            minel = length(corpusStrf{1}.meanTestR);
            rmphn(1:minel, 1)  = (corpusStrf{1}.meanTestR.^2)';
            rmwrd(1:minel, 1)  = (corpusStrf{6}.meanTestR.^2)';
            rmsurp(1:minel, 1) = (corpusStrf{3}.meanTestR.^2)';
            rmpr(1:minel, 1)   = (corpusStrf{4}.meanTestR.^2)';
            full(1:minel, 1)   = (corpusStrf{5}.meanTestR.^2)';
    
            % Compute unique variances for different features
            uv_word = full - rmwrd;
            uv_peakrate = full - rmpr;
            uv_phonetic = full - rmphn;
            uv_surp = rmsurp - rmwrd;
    
            % Combine unique variances into a single matrix
            uvs = [uv_word, uv_peakrate, uv_phonetic, uv_surp];
    
            % Find the indices in the word_encoding structure for the current SID
            tblidx = find(strcmp(word_encoding.SID, SID));
            
            % Filter out indices that exceed the size of uvs matrix
            els = word_encoding.el(tblidx);
            els = els(els < size(uvs, 1));
            
            % Assign the unique variance values to the corresponding fields in word_encoding
            word_encoding.(['uv_' corpus{l}])(tblidx(els < size(uvs, 1))) = ...
                mat2cell(uvs(els, :), ones(1, length(els)), 4);
            word_encoding.(['acsfeat_rsq_' corpus{l}])(tblidx(els < size(uvs, 1)))...
                = rmwrd(els);
            
            % Clear temporary variables to free up memory
            clear rm* full uv*
        else
            % Display a warning if STRF information is missing for the current SID
            warning([SID ' strf info missing.']);
        end
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% Running the logistic classification on neural data for different times around onset

% requires you to load word responses for all SIDs
Swrd = Dwrd;
window = 20; % 200ms window around boundary
Swrd.ambiguity = getAAI(window, Swrd, 'aud');

% english/spanish - [1, 2], mandarin - 3, bilingual - 4
ls = [1, 2]; %[1, 2];% [1, 2]; 
subtype = ''; % using different model, aud + word/wordOnset _word
trialset = '_highAAI'; % all trials, high ambiguity or low ambiguity _highAAI

if startsWith(Swrd.sentid{1}, 's')
    corpus = 'DIMEX';
    lowAAIthresh = 0.15;
    highAAIthresh = 0.5;
    column = 2;
else
    corpus = 'TIMIT';
    lowAAIthresh = 0.2;
    highAAIthresh = 0.65;
    column = 1;
end

idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));

% subset to both trialsets to ensure same subjects etc are used across sets
if any(strcmp(trialset, {'_highAAI', '_lowAAI'}))
    idx = idx & (Swrd.ambiguity>highAAIthresh | Swrd.ambiguity<lowAAIthresh);
end

% expanding window
% timing = [2, 5, 10, 20];
% startp = repmat(51, 4, 1);
% timelabel = '20-200ms';

% sliding window
% startp = 50:5:90;
% timing = repmat(4, 1, length(startp));
% timelabel = '50ms-sliding';

% single window
% startp = 51;
% timing = 47;
% timelabel = '500ms';

% single window before onset ****
startp = 31;
timing = 57;
timelabel = '600ms';

% sliding window
% startp = 30:5:90;
% timing = ones(length(startp), 1)*10;
% timelabel = '10ms-sliding';

% sliding window
% startp = 30:5:90;
% timing = ones(length(startp), 1)*5;
% timelabel = '5ms-sliding';

% subset to only the monolinguals
% mono_encoding = word_encoding(word_encoding.ls<3, :);

encoding_tbl = word_encoding;
word_uv = nan(height(encoding_tbl ), 1);
for i = 1:height(encoding_tbl )
    if ~isempty(encoding_tbl .(['uv_' lower(corpus)]){i})
        % true word unique variance versus word+surprise (1 or 4)
        word_uv(i) = encoding_tbl.(['uv_' lower(corpus)]){i}(4);
    end
end
wordAud_uv = encoding_tbl.uv_word(:, column);

switch subtype
% was using surprisal+word model
    case '_acs' % only "acoustic" electrodes
        tblidx = ismember(encoding_tbl.ls,ls) & (word_uv<0 & ~isnan(word_uv));
        mono_encoding = encoding_tbl (tblidx, :);
    case '_word' % only word and sequence electrodes
        tblidx = ismember(encoding_tbl .ls,ls) & (word_uv>0 & ~isnan(word_uv));
        mono_encoding = encoding_tbl (tblidx, :);
    case '_acsAud' % stringent on word_uv
        tblidx = ismember(word_encoding.ls,ls) & ...
            ((word_uv<0 & wordAud_uv<0) & ~isnan(word_uv));
        mono_encoding = word_encoding(logical(tblidx), :);
    case '_wordAud' % top 50 acs electrodes
        tblidx = ismember(word_encoding.ls,ls) & ...
            ((word_uv<0 | wordAud_uv>0) & ~isnan(word_uv));
        mono_encoding = word_encoding(logical(tblidx), :);
    case '_acsTop'
        thresh = 0.015;
        % threshold on acoustic feature rsq
        tblidx = ismember(encoding_tbl.ls,ls) & ...
            (word_uv<0 & ~isnan(word_uv)) & ...
            encoding_tbl.(['acsfeat_rsq_' lower(corpus)])>thresh;
        mono_encoding = encoding_tbl (tblidx, :);
    case '' % all
        mono_encoding = word_encoding(ismember(word_encoding.ls,ls), :);
end
clear word_uv

tps = startp(1):startp(end)+max(timing);
X_neural = nan(height(mono_encoding), length(tps), sum(idx));

% aggregation of neural data
disp('Loading all neural data....');
uSIDs = [unique(mono_encoding.SID)];
tic
for i = 1:length(uSIDs)
    disp(['loading subject ' num2str(i) ' of ' num2str(length(uSIDs))])
    SID = uSIDs{i};
    elidx = find(strcmp(mono_encoding.SID, SID));
    els = mono_encoding.el(elidx);

    % find all relevant neural responses
    tmpidx = cellfun(@(x) ~isempty(x), [Swrd.(SID)]);
    X_tmp = nan(size(Swrd.(SID){find(tmpidx, 1)}, 1), length(tps), ...
        height(Swrd));
    
    % include all non-nan responses
    tmp = cat(3, Swrd.(SID){:});
    X_tmp(:, :, tmpidx) = tmp(:, tps, :);
    % remove all sentence onset responses

    % make sure all electrodes are under maximum electrode number
    els = els(els<size(X_tmp, 1));
    elidx = elidx(els<size(X_tmp, 1));

    X_neural(elidx, :, :) = X_tmp(els, :, idx);
    clear X_tmp
end
toc
disp('Loading neural data complete!');

% logistic regression with 5-fold cross-validation
% idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));
if startsWith(Swrd.sentid{1}, 's')
    cols = brewermap(length(startp)+1, 'YlGnBu'); % dimex color map
else
    cols = brewermap(length(startp)+1, 'YlOrRd'); % timit color map
end

cols = cols(2:end, :);
if ls == 3
    subjs = {'Mandarin'};
    nfolds = 5;
elseif ls == 4
    subjs = {'Bilingual'};
    nfolds = 20;
else
    subjs = {'Spanish', 'English'};
    nfolds = 20;
end

if strcmp(trialset, '_highAAI') || strcmp(trialset, '_lowAAI')
    nfolds = 15; % highAAI / lowAAI have fewer trials
end

AUC = nan(length(subjs), length(timing), nfolds);
acc = cell(length(subjs), 1);
% electrode weighting through the logistic regression and PCA 
comp_weighted = cell(length(subjs), length(timing));
elecs = cell(length(subjs), length(timing));
scores = cell(length(subjs), length(timing));
yidx = cell(length(subjs), 1);

f = figure;
for ls = ls 

    for t = 1:length(timing)
        % reset idx for each subject group
        idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));
        % subset to both trialsets to ensure same subjects etc are used across sets
        if any(strcmp(trialset, {'_highAAI', '_lowAAI'}))
            idx = idx & (Swrd.ambiguity>highAAIthresh | Swrd.ambiguity<lowAAIthresh);
        end

        
        % assumes that X_neural starts at startp(1)
        start = startp(t)-startp(1)+1;
        tps = start:start+timing(t);

        elidx = mono_encoding.ls==ls;
        elshape = reshape(repmat(find(elidx), 1, length(tps)), ...
        sum(elidx)*length(tps), []);

        X = X_neural(elidx, tps, :);        
        y = Swrd.wordOns(idx)>0;        
         
        if ls == 3
            mintrl = 500;
        else
            mintrl = 1900;
        end            

        [X, nanrow, nancol, y] = makeDataMatrix(X, y, ...
            mono_encoding.SID(elidx), mintrl);
        tmp = find(idx);
        yidx{ls} = tmp(~nancol); % this won't be correct in AAI cases
        elshape(nanrow) = [];
        elecs(ls, t) = {unique(elshape)};

        idx(tmp(nancol))=0;
        switch trialset
            case '_highAAI'
                X = X(:, Swrd.ambiguity(idx)>highAAIthresh);
                idx = idx & Swrd.ambiguity>highAAIthresh;
                y = Swrd.wordOns(idx)>0;
            case '_lowAAI'
                X = X(:, Swrd.ambiguity(idx)<lowAAIthresh);
                idx = idx & Swrd.ambiguity<lowAAIthresh;
                y = Swrd.wordOns(idx)>0;
        end

         % Display how many electrodes and subjects are included in the final matrix
        disp(['Electrodes in ls ' num2str(ls) ...
            ', t' num2str(timing(t)) ': ' num2str(size(X, 1)/length(tps))])
        disp(['Unique subjects included: ' ...
            num2str(length(unique(mono_encoding.SID(elecs{ls, t}))))...
            ', trials: ' num2str(size(X, 2))])
 
        % Computes bootstrapped accuracy across 50 repeats of 80-20 train
        % test splits, with pca computed
        % Weights are significant model components x electrode
        [fp, tp, AUC(ls, t, :), ~, tmp_score, acc{ls}(:, t), weights] = ...
            logistic(X', y, 1, [], tps, nfolds);
        scores{ls, t} = nan(sum(idx), 1);
        if ~strcmp(trialset, '_highAAI') && ~strcmp(trialset, '_lowAAI')
            scores{ls, t}(~nancol) = tmp_score;
        end
        comp_weighted(ls, t) = {weights};

        % plot the current ROC curve for the last fold (as example)
        set(0, 'CurrentFigure', f)
        plot(fp, tp, 'LineWidth', 2, 'Color', cols(t, :)); hold on;   
    end

    h = refline(1, 0);
    h.LineStyle = '--';
    h.Color = 'k';
    xlabel('False positive rate') 
    ylabel('True positive rate')
    
    box off; set(gca, 'FontSize', 15);  
end

decode_details = struct();
decode_details.acc = acc; % accuracy for both subject group types
decode_details.timing = timing;
decode_details.startp = startp;
decode_details.weights = comp_weighted;
decode_details.AUC = AUC;
decode_details.elecs = elecs;
decode_details.encoding = mono_encoding; % for threshold 
decode_details.score = scores;
decode_details.y = Swrd.wordOns(idx)>0;
decode_details.yidx = yidx;

if isscalar(ls) && ls == 3
    subj = 'mandarinmono';
elseif isscalar(ls) && ls == 4
    subj = 'bilingual';
else
    subj = 'monolingual';
end

filename = [corpus '_word_decode_' subj '_' timelabel subtype trialset '.mat']; % dimex filename
save([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% Visualize the unique variance from word and others features

cols = getColorsCrossComp(1);

ctr=1;
featlabels = {'peak rate', 'phonetic features'};
field = 'uv_dimex';

plt_title = 'DIMEX';

figure;
for feat = 2:3    
    for ls = [1, 2, 4]
        subplot(2, 3, ctr)
        lsidx = word_encoding.ls==ls;
        uvs_timit = cell2mat(word_encoding.(field)(lsidx));
        idx = uvs_timit(:, 1)>0 & uvs_timit(:, feat)>0;

        scatter(uvs_timit(idx, 1), uvs_timit(idx, feat), 45, ...
            cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.6); hold on;

        % Calculating the unique variance correlations
        [r, p] = corr(uvs_timit(idx, 1), ...
            uvs_timit(idx, feat), 'type', 'pearson');
        text(0.01, 0.08, [num2str(r, 2) ':' num2str(p, 2)])

        % Figure formatting
        ylabel(featlabels{feat-1});
        xlabel('word+surprisal');
        ylim([0.00001 0.1]);
        xlim([0.00001 0.1]);

        % Adding reference lines
        h=refline(1, 0);
        h.LineWidth = 1.5;
        h.Color = 'k';

        ctr = ctr+1;
    end
end

sgtitle(plt_title);

figure;

% look at weights across low and high ambiguity decoders



clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;
%% Visualize decoding accuracy results (fixed, 500 ms)

Swrd = Dwrd;
ls = [1 2];
subj = 'monolingual';
type = '_word';

if isscalar(ls)
    switch ls 
        case 3
            subj = 'mandarinmono';
        case 4
            subj = 'bilingual';    
    end
end

if startsWith(Swrd.sentid{1}, 's')
    corpus = 'DIMEX';
else
    corpus = 'TIMIT';
end
time_label = '500ms';

filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
% idx = Swrd.sentOns<1;

load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

styles = {':', '-', '', '-'};
% subjs = {'Spanish', 'English'};
startp = decode_details.startp;
if startsWith(Swrd.sentid{1}, 's')
    cols = brewermap(length(startp)+1, 'YlGnBu'); % dimex color map
else
    cols = brewermap(length(startp)+1, 'YlOrRd'); % timit color map
end
cols = cols(2:end, :);

figure('Position', [100, 300, 150, 300]);
for l = ls
    accls = squeeze(decode_details.AUC(l, :, :));
    for t = 1:length(startp) % t
        if length(startp)<2
            boxplot(accls(:, t)*100, 'Position', l, 'Color', cols(t, :), ...
                'BoxStyle', 'filled'); hold on;
            xlabel('Subject Group');
        else
            boxplot(accls(:, t)*100, 'Position', t, 'Color', cols(t, :), ...
                'BoxStyle', 'filled'); hold on;
            xlabel('Post-onset (ms)');
        end
    end
end

% Formatting
if strcmp(subj, 'monolingual')
    xticklabels({'Spanish', 'English'});
    xticks([1 2]);
    xlim([0.5 2.5]);
end
legend('off');

% for sliding window
% labels = startp+round(decode_details.timing(1)/2)-startp(1);
% xticklabels(split(num2str(labels*10)));
yline(0.5, 'Color', 'k');
yticks();
ylim([50 90]);
ylabel('AUC')    
set(gca, 'FontSize', 15);
box off;

% score distributions
% if ls == [1, 2]
%     [~, p]= ttest2(decode_details.AUC(2, :, :), decode_details.AUC(1, :, :));
%     text(1.5, 85, getSigStr(p), 'FontSize', 15); 
% 
%     % score distributions
%     figure;
%     subplot(2, 2, [1 3]);
%     isz = decode_details.score{1}==0 | decode_details.score{2}==0;
%     idx = find(idx);
%     scatter(decode_details.score{1}(~isz), ...
%         decode_details.score{2}(~isz), 55, Swrd.wordOns(idx(~isz)), ...
%         'filled', 'MarkerFaceAlpha', 0.4);
%     colormap([0.6 0 0.6; 0.1 0.7 0.2]);
%     h = lsline;
%     h.LineWidth = 2.5;
%     h.Color = [0.5 0.5 0.5];
%     yticks([0 0.5 1]);
%     xticks([0 0.5 1]);
%     xlabel('Score (Spanish subj)');
%     ylabel('Score (English subj)');
%     set(gca, 'FontSize', 15);
%     
%     subplot(2, 2, 2);
%     [f, x] = ksdensity(decode_details.score{1}(~isz & ~decode_details.y));
%     plot(x, f, 'LineWidth', 2, 'Color', [0.6 0 0.6]); hold on;
%     [f, x] = ksdensity(decode_details.score{1}(~isz & decode_details.y));
%     plot(x, f, 'LineWidth', 2, 'Color', [0.1 0.7 0.2]);
%     
%     ylabel('density');
%     xlim([-0.2 1.2])
%     xline(0); xline(1);
%     set(gca, 'fontsize', 15);
%     title('Spanish subjects');
%     xlabel('score');
%     
%     yyaxis right
%     histogram(decode_details.score{1}(~isz & ~decode_details.y), 20, ...
%         'FaceColor', [0.6 0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
%         'Normalization', 'Probability');
%     histogram(decode_details.score{1}(~isz & decode_details.y), 20, ...
%         'FaceColor', [0.1 0.7 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
%         'Normalization', 'Probability');
%     
%     subplot(2, 2, 4);
%     [f, x] = ksdensity(decode_details.score{2}(~isz & ~decode_details.y));
%     plot(x, f, 'LineWidth', 2, 'Color', [0.6 0 0.6]); hold on;
%     [f, x] = ksdensity(decode_details.score{2}(~isz & decode_details.y));
%     plot(x, f, 'LineWidth', 2, 'Color', [0.1 0.7 0.2]);
%     ylabel('density');
%     xlim([-0.2 1.2]);
%     xline(0); xline(1);
%     set(gca, 'fontsize', 15);
%     title('English subjects');
%     xlabel('score');
%     
%     yyaxis right
%     histogram(decode_details.score{2}(~isz & ~decode_details.y), 20, ...
%         'FaceColor', [0.6 0 0.6], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
%         'Normalization', 'Probability');
%     histogram(decode_details.score{2}(~isz & decode_details.y), 20, ...
%         'FaceColor', [0.1 0.7 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.3, ...
%         'Normalization', 'Probability');
% end

clsplt = [1, 5; 3, 7];
ctr=1;
for l = ls

    numel = 40;
    f = figure;
    
    weights = squeeze(mean(decode_details.weights{l}));

    % Cluster the data
    eva = evalclusters(weights','kmeans', 'CalinskiHarabasz', 'klist', 1:5);
    
    for c = 1:eva.OptimalK
        subplot(2, eva.OptimalK, c);

        clusterweights = weights(:, eva.OptimalY==c);
        prct = (squeeze(sum(abs(clusterweights), [1, 2]))...
            ./squeeze(sum(abs(weights), [1, 2])))*100;
        prctelec = (prct./size(clusterweights, 2));

        [~, sortidx] = sort(squeeze(sum(abs(clusterweights))), ...
        'descend');
        elnum = min(numel, sum(eva.OptimalY==c));
        imagesc(clusterweights(:, sortidx(1:elnum))');     
        ylabel('Electrodes'); 
        title(['Cluster: ' num2str(c) ', ' ...
            num2str(size(clusterweights, 2)) ':' num2str(prct, 2) '%'])
    
        % 
        yyaxis right
        plot(smoothdata(mean(clusterweights(:, sortidx(1:elnum)), 2)), ...
            'Color', [0.5 0.5 0.5], 'linewidth', 3);
        ylim([-0.004 0.003]);
        yticks([]);
        
        colormap(vik); set(gca, 'FontSize', 15);
        caxis([-0.003 0.003])
        brighten(0.5);

        % find absolute index of els
        elidx = decode_details.elecs{l}(eva.OptimalY==c);
        elidx = elidx(sortidx);
        SIDs = decode_details.encoding.SID(elidx);
        els = decode_details.encoding.el(elidx);

        % accumulate resp
        allresp = nan(numel, bef+aft+1, height(Dwrd));
        yidx = decode_details.yidx{l};
        for e = 1:elnum
            SID = SIDs{e};
            
            nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
            resp = nan(size(Swrd.(SID){1}, 1), bef+aft+1, height(Swrd));
            resp(:, :, ~nanidx) = cat(3, Swrd.(SID){~nanidx});

            % Weighted response
            trlidx = intersect(yidx, find(~nanidx));
            allresp(e, :, yidx) = resp(els(e), :, ...
                trlidx).*sum(abs(weights(:, e)));
        end
    
        % Weight response and plot ERP
        dummy = struct();
        SID = 'avg';
        dummy.(SID).resp = mean(allresp(:, :, yidx), 'omitnan');
        
        dummy.wordOns = Swrd.wordOns(yidx);
        dummy.syllOns = ones(length(yidx), 1);
    
        subplot(2, eva.OptimalK, c+eva.OptimalK);
        plotWordErp(dummy, SID, 1, [], f, getColors(3), 1, 0.5); hold on;
        yticks(0);
        ylabel('HGA (z-score)');
        set(gca, 'FontSize', 13)
    end

    figure;
    tmp = corr(weights, weights); 
    [~, sortidx] = sort(eva.OptimalY); 
    imagesc(tmp(sortidx, sortidx));
    colormap(inferno);
    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd* X_neural;

%% Visualize decoding accuracy across high and low AAI trials

alltbl = cell(2, 1);
for swrd = {TDwrd, Dwrd}
    Swrd = swrd{1};

    ls = [1 2];
    subj = 'monolingual';
    type = '';
    
    if isscalar(ls)
        switch ls 
            case 3
                subj = 'mandarinmono';
            case 4
                subj = 'bilingual';    
        end
    end
    
    if startsWith(Swrd.sentid{1}, 's')
        corpus = 'DIMEX';
    else
        corpus = 'TIMIT';
    end
    time_label = '600ms';
    trialsets = {'_lowAAI', '_highAAI'};
    weights = cell(1, 2);
    encoding_subtbl = cell(1, 2);
    
    figure('Position', [100, 300, 150, 300]);
    ctr=0;
    for trialset = trialsets
        filename =  [corpus '_word_decode_' subj '_' time_label type trialset{1} '.mat']; 
        % idx = Swrd.sentOns<1;
        
        load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
        
        styles = {':', '-', '', '-'};
        % subjs = {'Spanish', 'English'};
        startp = decode_details.startp;
        if startsWith(Swrd.sentid{1}, 's')
            cols = brewermap(length(startp)+1, 'YlGnBu'); % dimex color map
        else
            cols = brewermap(length(startp)+1, 'YlOrRd'); % timit color map
        end
        cols = cols(2:end, :);
        
        for l = ls
            accls = squeeze(decode_details.AUC(l, :, :));
            for t = 1:length(startp) % t
                if length(startp)<2
                    boxplot(accls(:, t)*100, 'Position', l+ctr, 'Color', cols(t, :), ...
                        'BoxStyle', 'filled'); hold on;
                    xlabel('Subject Group');
                else
                    boxplot(accls(:, t)*100, 'Position', t, 'Color', cols(t, :), ...
                        'BoxStyle', 'filled'); hold on;
                    xlabel('Post-onset (ms)');
                end
            end
            encoding_subtbl{l} = decode_details.encoding(decode_details.elecs{l}, :);
        end
        [~, p] = ttest2(squeeze(decode_details.AUC(1, :, :)), squeeze(decode_details.AUC(2, :, :)));
        text(ctr+1, 100, getSigStr(p, 2))
        weights{ctr+1} = decode_details.weights;
        ctr=ctr+3;
    end
    
    % % Formatting
    % if strcmp(subj, 'monolingual')
    %     xticklabels({'Spanish', 'English'});
    %     xticks([1 2]);
    %     xlim([0.5 2.5]);
    % end
    % legend('off');
    % 
    % % for sliding window
    % % labels = startp+round(decode_details.timing(1)/2)-startp(1);
    % % xticklabels(split(num2str(labels*10)));
    % yline(0.5, 'Color', 'k');
    % yticks();
    % ylim([50 90]);
    ylabel('AUC')    
    set(gca, 'FontSize', 15);
    box off;
    
    figure; 
    if strcmp(corpus, 'DIMEX')
        native_ls = 1;
        field = 'uv_dimex';
    else
        native_ls = 2;
        field = 'uv_timit';
    end
    
    low_x = zscore(squeeze(max(mean(abs(weights{1}{native_ls, :}),1)))); 
    high_y = zscore(squeeze(max(mean(abs(weights{4}{native_ls, :}), 1))));
    
    scatter(low_x, high_y, 25, low_x-high_y, 'filled', 'MarkerEdgeColor', 'k');
    colormap(puor);
    caxis([-3 3]);
    ylabel('high AAI weight (z)');
    xlabel('low AAI weight (z)');
    [r, p] =corr(low_x, high_y, 'type','Spearman');
    title(['rval=' num2str(r) ', pval=' num2str(p)]);
    xl = xlim();
    yl = ylim();
    xlim([min([xl, yl]), max([xl, yl])]);
    refline(1, 0);
    
    % bin by whether show stronger weight in high AAI case or low AAI case
    discdiff = discretize(low_x-high_y, [-5 -1, 1 5]);
    discdiff(discdiff==2) = nan; % remove the 0s
    discdiff(discdiff==3) = 2;
    
    % bin 1 is greater high AAI weight, bin 3is greater low AAI weight
    subtbl = encoding_subtbl{native_ls}.(field);
    tbl = table();
    tbl.uv_wordsurp = cellfun(@(x) x(1), encoding_subtbl{native_ls}.(field));
    tbl.discdiff = discdiff;
    tbl.realdiff = low_x-high_y;
    tbl.SID = encoding_subtbl{native_ls}.SID;
    tbl.el = encoding_subtbl{native_ls}.el;
    % remove the nans
    tbl = tbl(~isnan(tbl.discdiff), :);
    
    figure; 
    h=boxchart(tbl.discdiff, tbl.uv_wordsurp); hold on;
    h.MarkerStyle = 'none';
    
    [~, p] = ttest2(tbl.uv_wordsurp(tbl.discdiff==1), tbl.uv_wordsurp(tbl.discdiff==2));
    text(1.5, 0.5, getSigStr(p, 2), 'FontSize', 15);
    
    % also scatter the randomly jittered points over the boxchart
    puor_tmp = puor;            
    for i = 1:2
        idx = tbl.discdiff==i;
        x = ones(sum(idx), 1)*i -0.05 + 0.1*randn(size(tbl.discdiff(idx)));
        y = tbl.uv_wordsurp(idx);
    
        % color by true difference
        scatter(x, y, 25, tbl.realdiff(idx), 'filled');
        colormap(puor);
        caxis([-3 3]);
    end
    xticks([1, 2]);
    xticklabels({'high AAI elecs', 'low AAI elecs'});
    title(getSigStr(p, 2));
    set(gca, 'FontSize', 13);

    alltbl{native_ls} = tbl;
end
alltbl = [alltbl{1};alltbl{2}];

%
% plotting MNI scatter
desel=struct();
desel.cols = [puor_tmp(30, :); puor_tmp(220, :)];
desel.labels = {'-1', '1'};
desel.conds = [1, 2]; % change to accomodate subject groups
desel.sz = repmat(30, 1, 2);
for s = unique(alltbl.SID)'
    SID = s{1};
    idx = strcmp(alltbl.SID, SID);

    desel.(SID).condition = alltbl.discdiff(idx);
    desel.(SID).elid = alltbl.el(idx); 
end
[mni_lh] = plotMNIElec(unique(alltbl.SID), desel, 'lh', 0);
light("Style","infinite","Position",[100 100 0]);

[mni_rh] = plotMNIElec(unique(alltbl.SID), desel, 'rh', 0);
light("Style","infinite","Position",[-100 100 0]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% Visualize decoding accuracy across all subject groups (fixed, 600 ms)

figure('Position', [100, 300, 150, 300], 'renderer', 'painters');
Swrds = {Dwrd, TDwrd};
decode_type = 'word'; % 'word' 
ctr=1;
for s = Swrds
    Swrd = s{1};
    time_label = '600ms';
    % time_label = '10ms-sliding';
    lss = [1, 2, 4];
    type = '';
    colors = getColorsCrossComp(1);
    
    switch time_label
        case '10ms-sliding'
            accls = nan(4, 13, 20);
        otherwise
            accls = nan(4, 20);
    end 
    
    subplot(2, 1, ctr)
    
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
        load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
        styles = {':', '-', '', '-'};
    %     subjs = {'Spanish', 'English'};
    %     startp = decode_details.startp;
    
        switch time_label
            case '10ms-sliding'
                accls(ls, :, :) = squeeze(decode_details.AUC(ls, :, :));
            otherwise
                accls(ls, :) = squeeze(decode_details.AUC(ls, :, :));
        end 
    end
    
    % if using sliding window use maximum timepoint across subject groups
    if strcmp(time_label, '10ms-sliding')
        [~, maxtp] = max(mean(accls, [1 3], 'omitnan'));
        accls = squeeze(mean(accls(:, maxtp-4: maxtp+4, :), 2));
    end
    
    if ctr == 1 % dimex
        position = [1, 4, 0, 2];
    else
        position = [4, 1, 0, 2];
    end

    for ls = lss
        h=boxchart(ones(size(accls, 2), 1)*position(ls), ...
            accls(ls, :)*100, 'BoxFaceColor', colors(ls, :)); 
        h.JitterOutliers = 'on';
        h.MarkerStyle = '.';
        h.MarkerColor = colors(ls, :);
        hold on; % 'BoxStyle', 'filled'
        
    end  

    xticks();
    xlim([0.5 4.5]);
    legend('off');
    
    %pos = [1, 2; 1, 3; 2, 3;];
    combo = [1, 2; 1, 4; 2, 4;];
    pos = nan(size(combo, 1), size(combo, 2));
    for r = 1:size(combo, 2)
        for c = 1:size(combo, 1)
            pos(c, r) = position(combo(c, r));
        end
    end

    for i = 1:3        
        [h, p] = ttest2(accls(combo(i, 1), :), accls(combo(i, 2), :));
        if ~isempty(getSigStr(p))
            plot([pos(i, 1)-0.02 pos(i, 2)+0.02], ...
                [70+5*i 70+5*i], 'LineWidth', 1.5, 'Color', 'k');
            text(mean(pos(i, :))-0.25, 70+5*i+1, getSigStr(p, 2), 'FontSize', 15);
        end
    end
    
    % for sliding window
    % labels = startp+round(decode_details.timing(1)/2)-startp(1);
    % xticklabels(split(num2str(labels*10)));
    
    yticks();
    yline(0.5, 'Color', 'k');
   
    ylim([48 90]);
    ylabel('AUC')    
    set(gca, 'FontSize', 15);
    box off;
    ctr = ctr+1;
end

%formatting
xticks([1 2 4]);
xlabel('Subject Group');
xticklabels({'Monolingual', 'Bilingual', 'Unfamiliar' });

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% Visualize decoding accuracy native vs. unfamiliar (fixed, 600 ms)

figure('Position', [100, 300, 150, 300], 'renderer', 'painters');

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
        load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
        switch time_label
            case '10ms-sliding'
                accls{ctr}(ls, :, :) = squeeze(decode_details.AUC(ls, :, :)); % Stores AUC values in the matrix
            otherwise
                accls{ctr}(ls, :) = squeeze(decode_details.AUC(ls, :, :));
        end 
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
labels = {{'native', 'unfamiliar'}, {'monolingual', 'bilingual'}};

cols = {'b', 'r'};
for t = 1:2 % Type
    subplot(1, 2, t);
    plt = plts{t};
    for c = 1:2
        h = boxchart(ones(size(plt, 2), 1)*c, plt(c, :), ...
            'BoxFaceColor', [0.5 0.5 0.5]); % Creates boxplots
        h.JitterOutliers = 'on';
        h.MarkerStyle = '.';
        h.MarkerColor = 'k';
        hold on;
        scatter(randn(size(plt, 2), 1)*0.1+c-0.05, plt(c, :), 25, 'filled', ...
            cols{c})

        hold on;
    end

    % Perform a Wilcoxan ranksum on the native vs. non-native decoding
    [p, ~] = ranksum(plt(1, :), plt(2, :)); 
    line([1.25 1.75], [.85, .85], 'Color', 'k', 'LineWidth', 1.5); % Draws a line
    text(1.35, .87, getSigStr(p, 2), 'FontSize', 13); % Adds text to the plot
    
    % Formatting
    xlabel('Group Type');
    ylabel('AUC');

    ylim([0.45 0.9]);
    yticks([0.5 0.9])
    yline(0.5, 'Color', 'k');

    xlim([0.5 2.5])
    xticks([1 2]);  
    xticklabels(labels{t})
    
    set(gca, 'FontSize', 15);
    box off;
end

figure; 
bil = [accls{1}(4, :), accls{2}(4, :)];
plt = [native; bil; nonnative];
pos = [1, 2, 4]; % X positions for boxplots
labels = {'monolingual', 'bilingual', 'non-native'};
for c = 1:3
    h = boxchart(ones(size(plt, 2), 1)*pos(c), plt(c, :), ...
        'BoxFaceColor', [0.5 0.5 0.5]);

    disp([labels{c} ': ' num2str(mean(plt(c, :))) ', std:' num2str(std(plt(c, :)))])
    
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';
    hold on;
end

combs = [1, 2; 2, 3; 1, 3];
for c = 1:3
    [p, ~] = ranksum(plt(combs(c, 1), :), plt(combs(c, 2), :)); % Performs a rank sum test
    line([pos(combs(c, 1)) pos(combs(c, 2))], [75+c*5, 75+c*5], 'Color', 'k', 'LineWidth', 1.5);
    text(mean([pos(combs(c, 1)) pos(combs(c, 2))]), 76+c*5, getSigStr(p, 2), 'FontSize', 13); % Adds text to the plot
end

% Formatting
xlabel('Group Type');
ylabel('AUC');

ylim([50 95]);
yticks([50 90])
yline(0.5, 'Color', 'k');

xticks([1, 2, 4]);  
xticklabels(labels)

set(gca, 'FontSize', 15);
box off;

% another version of the plot
figure; 
bil = [accls{1}(4, :), accls{2}(4, :)];
plt = {[native, bil], nonnative};
pos = [1, 2]; % X positions for boxplots
labels = {'native', 'unfamiliar'};
for c = 1:2
    h = boxchart(ones(size(plt{c}, 2), 1)*pos(c), plt{c}*100, ...
        'BoxFaceColor', [0.5 0.5 0.5]);
    
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';
    hold on;
end

% add bilingual and monolingual to single box plot
scatter(randn(size(native))*0.1+pos(1)-0.1, native*100, 25, 'k', ...
    'filled', 'HandleVisibility', 'off');
scatter(randn(size(bil))*0.1+pos(1)+0.1, bil*100, 35, [0.5 0.5 0.5], 'HandleVisibility', 'off');

combs = [1, 2;];
for c = 1
    [p, ~] = ranksum(plt{1}, plt{2}); % Performs a rank sum test
    line([pos(combs(c, 1)) pos(combs(c, 2))], [75+c*5, 75+c*5], 'Color', 'k', 'LineWidth', 1.5);
    text(mean([pos(combs(c, 1)) pos(combs(c, 2))]), 76+c*5, getSigStr(p, 2), 'FontSize', 13); % Adds text to the plot
end

% Formatting
xlabel('Group Type');
ylabel('AUC');

ylim([45 90]);
yticks([50 90])
yline(0.5, 'Color', 'k');

xticks([1, 2, 4]);  
xticklabels(labels)

set(gca, 'FontSize', 15);
box off;


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% Visualize decoding accuracy all subject groups (sliding, 10 ms)

addpath(genpath('shadederror/'))
subtype = '';

switch subtype
    case '_acs'
        linestyle = '--';
     case '_acsTop'
        linestyle = '--';
    case ''
        linestyle = '-';
    case '_word'
        linestyle = '-.';
end

figure;
ctr=1;
for c = {'DIMEX', 'TIMIT'}
    
    corpus = c{1};
    disp(['-----------' corpus '--------------']);
    filename = [corpus '_word_decode_monolingual_5ms-sliding' subtype '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    x = (decode_details.startp-50)./100;

    subplot(2, 1, ctr)
    colors = getColorsCrossComp(1);
    data = nan(2, length(x), size(decode_details.AUC, 3));
    for ls = 1:2
        data(ls, :, :) = smoothdata(squeeze(decode_details.AUC(ls, :, :)), ...
            'SmoothingFactor', 0.1).*100;
        %data = squeeze(decode_details.AUC(ls, :, :));
        y = squeeze(mean(data(ls, :, :), 3));
        err = squeeze(nansem(data(ls, :, :), 3));
        shadedErrorBar(x, y, err, ...
            'lineprops', {'color', colors(ls, :), 'linewidth', 2, ...
            'linestyle', linestyle}); hold on;

        elecs = decode_details.elecs{ls};
        disp(['ls: ' num2str(ls) ', elecs:' num2str(length(elecs)) ...
            ', num subj: ' num2str(length(unique(decode_details.encoding.SID(elecs))))]);

        % print max peak +- values
        maxtp = argmax(mean(squeeze(data(ls, :, :)), 2));
        disp(datastats(data(ls, maxtp, :)))
    end

    X = reshape(cat(2, squeeze(data(1, :, :)), squeeze(data(2, :, :))), 1, ...
        length(x), size(data, 3)*2);
    y = [ones(1, size(data, 3)) 2*ones(1, size(data, 3))];
    [fvals, ~, ~, df1, df2] = Fstat_TIMIT(X, y, [1, 2]);
    fthresh = finv(1-0.01, df1, df2);  

    scatter(x(fvals>fthresh), 52*ones(1, sum(fvals>fthresh)), 45, ...
        int32(fvals(fvals>fthresh))', 'filled', 'HandleVisibility', 'off');
    cm = colormap("gray");
    colormap(flipud(cm(1:200, :)))

    ylim([45 75]);
    yticks((0.5:0.1:0.8)*100);
    set(gca, 'FontSize', 15);
    
    h=title(['Listening to ' corpus]);
    h.FontWeight='normal';
    xline(0, '-k', 'linewidth', 1, 'handlevisibility', 'off');
    yline(50, '-k', 'linewidth', 1, 'handlevisibility', 'off');
    xlabel('Time (s)')
    ylabel('AUC');

    ls=4;
    filename = [corpus '_word_decode_bilingual_10ms-sliding' subtype '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
 
    data(4, :, :) = smoothdata(squeeze(decode_details.AUC(ls, :, :)), ...
        'SmoothingFactor', 0.01).*100;
    y = squeeze(mean(data(4, :, :), 3));
    err = squeeze(nansem(data(4, :, :), 3));
    shadedErrorBar(x, y, err, ...
        'lineprops', {'color', colors(ls, :), 'linewidth', 2, ...
        'linestyle', linestyle}); hold on;
    legend({'Spanish', 'English', 'Bilingual'});
    elecs = decode_details.elecs{ls};
        disp(['ls: ' num2str(ls) ', elecs:' num2str(length(elecs)) ...
            ', num subj: ' num2str(length(unique(decode_details.encoding.SID(elecs))))]);

    % bilingual versus native monolingual
    X = reshape(cat(2, squeeze(data(ctr, :, :)), squeeze(data(4, :, :))), 1, ...
        length(x), size(data, 3)*2);
    y = [ones(1, size(data, 3)) 2*ones(1, size(data, 3))];
    [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(X, y, [1, 2]);
    fthresh = finv(1-0.001, df1, df2);  

    scatter(x(fvals>fthresh), 50*ones(1, sum(fvals>fthresh)), 45, ...
        int32(fvals(fvals>fthresh)), 'filled', 'HandleVisibility', 'off');
    cm = colormap("gray");
    colormap(flipud(cm(1:200, :)))

%     if strcmp(corpus, 'TIMIT')
%         ls=3;
%         filename = [corpus '_word_decode_mandarinmono_10ms-sliding.mat']; 
%         load(filename);
%         data = smoothdata(squeeze(decode_details.AUC(ls, :, :)));
%         y = squeeze(mean(data, 2));
%         err = squeeze(nansem(data, 2));
%         shadedErrorBar(x, y, err, ...
%             'lineprops', {'color', colors(ls, :), 'linewidth', 1.5}); hold on;
%         legend({'Spanish', 'English', 'Bilingual', 'Mandarin'});
%     end
    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd* X_neural;


%% Visualize acoustic versus word-level electrode decoding (fixed, 600 ms)
% the full decoding model (fixed, 600 ms) 

time_label = '600ms';
% time_label = '10ms-sliding';
lss = [1, 2];

Swrds = {Dwrd, TDwrd};
types = {'_acsTop', '_word'};
colors = getColorsCrossComp(1);
position = [1, 2, 0, 3];

accall = cell(length(types), length(Swrds));

% Collect AUC values for each subject group, type, and time_label
for s = 1:length(Swrds)
    Swrd = Swrds{s};
    for t = 1:length(types)
        type = types{t};
        
        switch time_label
            case '10ms-sliding'
                accls = nan(4, 13, 20);
            case '500ms'
                accls = nan(4, 20);
        end     

        for ls = lss
            % Determine the subject group and corpus
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
            
            % Load AUC values based on the specified time_label and type
            filename = [corpus '_word_decode_' subj '_' time_label type '.mat']; 
            load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

            switch time_label
                case '10ms-sliding'
                    accls(ls, :, :) = squeeze(decode_details.AUC(ls, :, :));
                case '500ms'
                    accls(ls, :) = squeeze(decode_details.AUC(ls, :, :));
                case '600ms'
                    accls(ls, :) = squeeze(decode_details.AUC(ls, :, :));
            end 

            % If using sliding window, use the maximum timepoint across subject groups
            if strcmp(time_label, '10ms-sliding')
                [~, maxtp] = max(mean(accls, [1 3], 'omitnan'));
                accls = squeeze(mean(accls(:, maxtp-4: maxtp+4, :), 2));
            end
            
            accall(s, t) = {accls};
        end
    end
end

% For linear mixed effect model
nfolds = size(accls, 2);

ctr=1;
AUC = nan(nfolds*length(Swrds)*length(types), 1);
comp = nan(nfolds*length(Swrds)*length(types), 1);
type = nan(nfolds*length(Swrds)*length(types), 1);
lang = nan(nfolds*length(Swrds)*length(types), 1);

for s = 1:length(Swrds)
    figure('Renderer', 'Painters'); 
    
    for t = 1:length(types)
        accls = accall{s, t};
        
        for ls = lss
            h = boxchart(ones(size(accls, 2), 1)*position(ls)+(t-1)*3, ...
                accls(ls, :)*100, 'BoxFaceColor', colors(ls, :)); 
            h.JitterOutliers = 'on';
            h.MarkerStyle = '.';
            h.MarkerColor = colors(ls, :);
            hold on;
            xlabel('Subject Group');
        end

        % Formatting
        xticklabels({'Spanish', 'English', 'Bilingual' });
        xticks([1 2 3]);
        xlim([0.5 5.5]);
        legend('off');

        % For sliding window
        % labels = startp+round(decode_details.timing(1)/2)-startp(1);
        yline(0.5, 'Color', 'k');
        yticks();
        ylim([50 90]);
        ylabel('AUC')    
        set(gca, 'FontSize', 15);
        box off;
        
        % For linear mixed effect model
        AUC(ctr:ctr+nfolds*2-1) = [accls(1, :)'; accls(2, :)';];
        % Corpus
        comp(ctr:ctr+nfolds*2-1) = [ones(nfolds, 1)*double(s==1); ones(nfolds, 1)*double(s==2);];
        % Subject language background
        lang(ctr:ctr+nfolds*2-1) = [ones(nfolds, 1); ones(nfolds, 1)*2;];
        % Electrode subset
        type(ctr:ctr+nfolds*2-1) = ones(nfolds*2, 1)*t;
        ctr = ctr+nfolds*2;

        pos = [position(1)+(t-1)*3, position(2)+(t-1)*3];
        combo = [1, 2];
        
        for i = 1
            [p, ~] = ranksum(accls(combo(i, 1), :), accls(combo(i, 2), :));
            
            if ~isempty(getSigStr(p))
                plot([pos(i, 1)-0.02 pos(i, 2)+0.02], ...
                    [80+5*i 80+5*i], 'LineWidth', 1.5, 'Color', 'k');
                text(mean(pos(i, :))-0.25, 80+5*i+1, getSigStr(p, 2), 'FontSize', 15);
            end
        end
    end
    
    xticks([1.5 4.5])
    xlabel('Electrode set')
    xticklabels({'acoustic only', 'word boundary only'})
    legend({'Spanish mono', 'English mono'}, 'Location','best')
    yticks([50 90])
end

tbl = table();
tbl.AUC = AUC; 
tbl.comp = comp; 
tbl.type = type;
tbl.lang = lang;

lme = fitlme(tbl,'AUC~1+comp+type+lang');
disp(lme)

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* *elecs;
%% FIX: misclassification examples

% misclassified as a syllable but actually a word
% misclass = decode_details.yidx(decode_details.y==1); % decode_details.score{1} <0.03 & 
% subplot(2, 2, 1)
% cols = inferno(9);
% cols = cols(2:end, :);
% for i = 2:8
%     idx = misclass_tbl.wordLen==i;
%     scores = decode_details.score{ls}(decode_details.y==1);
%     [f,xi] = ksdensity(scores(idx), 'Bandwidth', 0.05);
%     plot(xi, f+i*0.7, 'Color', cols(i, :), 'LineWidth', ...
%         sum(~isnan(scores(idx)))/70); hold on;
%     disp(sum(~isnan(scores(idx))));
%     xlim([0 1]);
% 
%     ylabel('word length (# phones)');
%     xlabel('score');
%     set(gca, 'FontSize', 15);
% end
figure;
label ={'Spanish monolingual', 'English monolingual'};
for ls = 1:2
    prctiles = [prctile(decode_details.score{ls}(decode_details.y==1), 10), ...
        prctile(decode_details.score{ls}(decode_details.y==1), 90)];
    
    misclass = decode_details.yidx(decode_details.y==1 & ...
        decode_details.score{ls} <prctiles(1));
    misclass_tbl = Swrd(misclass, :);
    
    corrclass = decode_details.yidx(decode_details.y==1 & ...
        decode_details.score{ls} >prctiles(2));
    corrclass_tbl = Swrd(corrclass, :);
    % 
    % boxplot([misclass_tbl.wordLen, corrclass_tbl.wordLen]);
    % boxplot([misclass_tbl.precWordLen, corrclass_tbl.precWordLen]);
    % histogram(misclass_tbl.syll); hold on; histogram(corrclass_tbl.syll)
    subplot(3, 2, ls)
    colors = [0.6 0 0.6; 0.1 0.7 0.2];
    
    shadedErrorBar(-0.5:0.01:0.5, mean(cat(1, misclass_tbl.env{:})), ...
        nansem(cat(1, misclass_tbl.env{:})), {'Color', colors(1, :), ...
        'LineWidth', 1.5}, 0.5); hold on;
    
    shadedErrorBar(-0.5:0.01:0.5, mean(cat(1, corrclass_tbl.env{:})), ...
        nansem(cat(1, corrclass_tbl.env{:})), {'Color', colors(2, :), ...
        'LineWidth', 1.5}, 0.5); hold on;
    title(label{ls})
    box off;
    set(gca, 'FontSize', 13);
    ylabel('Amplitude');
    xlabel('Time (s)');
    xline(0);
    ylim([-0.25 1.3]);


    subplot(3, 2, ls+2);
    imagesc(-0.5:0.01:0.5, 1:80, median(cat(3, misclass_tbl.aud{:}), 3));
    set(gca, 'YDir', 'normal');
    title('Misclassified as syllable');
    xline(0)
    
    subplot(3, 2, ls+4);
    imagesc(-0.5:0.01:0.5, 1:80, median(cat(3, corrclass_tbl.aud{:}), 3));
    title('Correctly classified as word')
    xline(0)
    set(gca, 'YDir', 'normal');
    colormap(inferno)

end

figure;
histogram(corrclass_tbl.precPhn); hold on; histogram(misclass_tbl.precPhn)
figure; 
histogram(corrclass_tbl.onsPhn); hold on; histogram(misclass_tbl.onsPhn)

misclass = decode_details.yidx(decode_details.y==1);
misclass_tbl = Dwrd(misclass, :);

[~, idx] = mink(decode_details.score{2}(decode_details.y==1), 5);
arrayfun(@(i) join(misclass_tbl.phns{i}), idx)

%% Visualize word erps with overlaid decoding weights

Swrd = Dwrd;
lss = 1;
if startsWith(Swrd.sentid{1}, 's')
    corpus = 'DIMEX';
else
    corpus = 'TIMIT';
end

switch lss
    case 4
        filename =  [corpus '_word_decode_bilingual_600ms.mat']; 
    otherwise % ls = [1, 2]
        filename =  [corpus '_word_decode_monolingual_600ms.mat']; 
end

% Dwrd = loadDwrdMulti('dimex',  bef, aft, {'EC163', 'EC172', ...
%     'EC214', 'EC219', 'EC222', 'EC129', 'EC260', 'EC266'}, dimex_details);
% TDwrd = loadDwrdMulti('timit', bef, aft, {'EC172', 'EC214', ...
%     'EC212', 'EC222', 'EC129', 'EC260', 'EC266'}, timit_details);

decodeSID = {'EC163'};
% Dwrd = loadDwrdMulti('dimex',  bef, aft, decodeSID, dimex_details);
% TDwrd = loadDwrdMulti('timit', bef, aft, decodeSID, timit_details);

load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

% find top 3 weighted electrodes and look at word erps
f = figure; 
ctr = 1;
numel = 20;
for ls = lss
    weights = squeeze(median(decode_details.weights{ls}, 1));

    % Either show electrodes with highest decoding weights
    %[~, maxel] = maxk(sum(abs(weights)), numel);

    elidx = decode_details.elecs{ls};
    maxel = find(strcmp(decode_details.encoding.SID(elidx) , decodeSID{1}));
    maxel = maxel(1:numel);
    %maxel = [7, 8, 25];

    for e = maxel'
        SID = decode_details.encoding.SID{elidx(e)};
        %disp(['SID:' SID ', ' num2str(decode_details.encoding.el(elidx(e)))])
    
        dummy = struct(); 
        nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
            %| (Swrd.syll<2 & ~isnan(Swrd.syll));
        dummy.(SID).resp = cat(3, Swrd.(SID){~nanidx});
        
        dummy.wordOns = Swrd.wordOns (~nanidx);
        dummy.syllOns = ones(sum(~nanidx), 1);
    
        subplot(2, numel, ctr)
        addpath(genpath('shadederror'))
        plotWordErp(dummy, SID, decode_details.encoding.el(elidx(e)), ...
            [], f, [0 0 0 ;0.6 0.1 0.28], 1, 0.5); hold on;
        %[0.6 0 0.6;0.1 0.7 0.2]
        yticks(0);
        ylabel('HFA (z)');
        set(gca, 'FontSize', 13);
        yyaxis right
        tp = -0.21:0.01:0.36; %0.01:0.01:0.48
        plot(tp, squeeze(weights(:, e)), 'LineWidth', 2, ...
            'LineStyle', ':', 'Color', 'k');
        yticks(0);
        legend('off');
        xlim([-0.2, 0.4]);

        el = decode_details.encoding.el(elidx(e));
        [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
            dummy.(SID).resp(el, :, :), dummy.wordOns+1, [1, 2]);
        fthresh = finv(1-0.0001, df1, df2);  

        x = -0.5:0.01:0.5;
        yyaxis left
        scatter(x(fvals>fthresh), 0*ones(1, sum(fvals>fthresh)), 15, ...
            fvals(fvals>fthresh), 'filled');
        cm = colormap("gray");
        colormap(flipud(cm(1:200, :)))

        ctr = ctr + 1;
    end
    clear dummy
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;


%% Plot word ERPs (syllable versus word) for monolingual example (no weights)

% find top 3 weighted electrodes and look at word erps
f = figure; 

% SIDs = {'EC260', 'EC260', 'EC260', 'EC260', 'EC260', 'EC260', ...
%     'EC260', 'EC260'}; % 'EC266', 'EC266'
% els = [ 221, 192, 205, 206, 207, 208, 217, 214]; %236, 205 for EC260, 178, 175

SIDs = {'EC183'};% 'EC100', 'EC100', 'EC100','EC100'};
els = 71;%,135, 22, 70, 71, 150]; %236, 205 for EC260, 178, 175

% SIDs = {'EC260'};% 'EC100', 'EC100', 'EC100','EC100'};
% els = 221;%,135, 22, 70, 71, 150]; %236, 205 for EC260, 178, 175
% 
% % in figure
% % EC100: 21, EC183: 135
% 
% els = [55]; %236, 205 for EC260, 178, 175
% SIDs = repmat({'EC163'}, length(els), 1);


numel = length(els);
Swrds = {Dwrd, TDwrd};
for s = 1:2
    Swrd = Swrds{s};

    for ctr = 1:length(els)
        SID = SIDs{ctr};    
    
        dummy = struct();
        nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
            %| (Swrd.syll<2 & ~isnan(Swrd.syll));
        dummy.(SID).resp = cat(3, Swrd.(SID){~nanidx});
        
        dummy.wordOns = Swrd.wordOns (~nanidx);
        dummy.syllOns = ones(sum(~nanidx), 1);
    
        subplot(2, numel, ctr+(s-1)*numel)
        addpath(genpath('shadederror'))
        plotWordErp(dummy, SID, els(ctr), ...
            [], f, getColorsCrossComp(3), 1, 0.5); hold on;
%         yticks(0);
        ylabel('HFA (z)');
        set(gca, 'FontSize', 13);

        [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
            dummy.(SID).resp(els(ctr), :, :), dummy.wordOns+1, [1, 2]);
        fthresh = finv(1-0.0001, df1, df2);  
    
        x = -0.5:0.01:0.5;
        scatter(x(fvals>fthresh), 0*ones(1, sum(fvals>fthresh)), 25, ...
            fvals(fvals>fthresh), 'filled', 'HandleVisibility', 'off');
        cm = colormap("gray");
        colormap(flipud(cm(1:200, :)))
        ylim([-0.1 1]);
    end

    clear dummy
end

corpusnames = {'dimex', 'timit'};
modelname = {'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin'};
figure; 
modelfeatures = ['word onset', 'word length', 'surprise'];
for ctr = 1:length(els)
    for s = 1:2

        SID = SIDs{ctr};   
%         Swrd = Swrds{s};
        % read in model strf weights               
        [strf] = loadMultModelStrf(SID, modelname(s), corpusnames{s}, ...
            datapath, 1);  
        weights = strf{1}.meanStrf(13:end, :, els(ctr));
        
        subplot(1, 2, s);
        imagesc(weights(:, 1:40));
        colormap(inferno)
    
%         yticks(1:length(modelfeatures));
%         yticklabels(modelfeatures);
        xlim([0.5 40])
        xticks([1 40]);
        xticklabels({'0', '-0.4'});
        xlabel('Time (s)');
        set(gca, 'FontSize', 13);
        clear strf
        title(els(ctr))
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: MNI Brain for decoding weights

corpus = 'DIMEX';
filename =  [corpus '_word_decode_bilingual_600ms.mat']; 

load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

% colors
numbins = 8;
colors = getColorsCrossComp(1);
colors_sp = [linspace(1,colors(1, 1),numbins)', ...
    linspace(1,colors(1, 2),numbins)', ...
    linspace(1,colors(1, 3),numbins)'];
colors_en = [linspace(1,colors(2, 1),numbins)', ...
    linspace(1,colors(2, 2),numbins)', ...
    linspace(1,colors(2, 3),numbins)'];
% colors_b = [linspace(1,colors(4, 1),numbins)', ...
%     linspace(1,colors(4, 2),numbins)', ...
%     linspace(1,colors(4, 3),numbins)'];

% colors_b = [linspace(1,0.2,numbins)', ...
%     linspace(1,0.7,numbins)', ...
%     linspace(1,0.2,numbins)'];
% colors_b = [linspace(1,0.6,numbins)', ...
%     linspace(1,0.1,numbins)', ...
%     linspace(1,0.28,numbins)'];
colors_b = [1 1 1 ; flipud(fpurple(numbins))];

ctr = 1;
numsubj = 3; % 3 subject groups

desel=struct();
desel.cols = [colors_sp; colors_en; colors_b];
desel.labels = split(num2str(1:numbins*numsubj));
desel.conds = 1:numbins*numsubj; % change to accomodate subject groups
desel.sz = [10, repmat(45, 1, numbins*numsubj-1)];
for ls = [1:2 4]
    % initialize design electrode structure   

    if ls == 4
        filename =  [corpus '_word_decode_bilingual_600ms.mat'];
        load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    end

    weights = squeeze(sum(abs(mean(decode_details.weights{ls})), 2));
    [~, binedges] = discretize(weights, numbins);
    for e = 1:length(decode_details.elecs{ls})
        el = decode_details.elecs{ls}(e);
        SID = decode_details.encoding.SID{el};
        if ~isfield(desel, SID)
            desel.(SID).elid = [];
            desel.(SID).condition = [];
        end
        desel.(SID).elid = [desel.(SID).elid decode_details.encoding.el(el)];
        desel.(SID).condition = [desel.(SID).condition discretize(weights(e), ...
            binedges)+(ctr-1)*numbins];
    end
    ctr = ctr+1;
end

names = fieldnames(desel);
sids = {'EC260', 'EC266'}; % ,%names(startsWith(names, 'EC'));
% desel.EC260.selid = [221, 192]; % 221 192
% desel.EC266.selid = [177, 164];

[native_plot] = plotNativeElec(sids, desel, 1);

%[mni_lh] = plotMNIElec(sids, desel, 'lh', 0);
%print(fullfile('mniBrain_LH_decodeweights.jpg'), '-djpeg', '-painters', '-r600')
%[mni_rh] = plotMNIElec(sids, desel, 'rh', 0);
%print(fullfile('mniBrain_RH_decodeweights.jpg'), '-djpeg', '-painters', '-r600')

% scatter(1:8, ones(10, 1), 55,1:10);
ctr = 1;
for ls = [1:2 4]
    figure;
    colormap(desel.cols((1:numbins)+(ctr-1)*numbins, :));
    cbh = colorbar;
    cbh.Ticks = [0 1];
    cbh.TickLabels = {'min', 'max'};
    axis off;
    set(gca, 'FontSize', 15);
    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd ;

%% comparing decoding weights to encoding uv

corpora = {'DIMEX', 'TIMIT'};
fields = {'sp_uv_all', 'eng_uv_all'};
r = nan(2, 5);
p = nan(2, 5);

figure;
for corp = 1:2
    corpus = corpora{corp};
    filename =  [corpus '_word_decode_monolingual_600ms.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
    for ls = corp % only show the case where its native
        weights = squeeze(mean(decode_details.weights{ls}, 1));
        avg_weights = mean(abs(decode_details.weights{ls}), [1, 3]);
        [~, max_tp] = max(avg_weights);

        % extract decoding weights
        tbl = decode_details.encoding(decode_details.elecs{ls}, :);
        tbl.mdl_weight = mean(abs(weights(max_tp-5:min(max_tp+5, 58), :)))';

        % remove electrodes if they do not contribute to the model (?)
        encode_tbl = innerjoin(tbl, wordsurp_encoding, 'Keys', {'SID', 'el'});
        featureOrd = {'onset', 'peakrate', 'formant', 'consonant', ...
            'word+surp', 'word', 'surp'};
        
        for index = 1:5
            idx = encode_tbl.(fields{corp})(:, index)>0;
            x = encode_tbl.mdl_weight(idx);
            y = encode_tbl.(fields{corp})(idx, index);
        
            % correlate the encoding table unique variance and the decoding
            % model weight
            [r(corp, index), p(corp, index)] = corr(x, y, 'Type', 'Spearman'); 

            subplot(2, 5, index+(corp-1)*5);
            scatter(x, y, 10, 'k', 'filled');
            l = lsline;
            l.LineWidth = 1.8;
            yticks([]);
            xticks([]);
            title(getSigStr(p(corp, index), 2));
            xlabel('model weight');
            ylabel([featureOrd{index} '\Delta R^2 ']);
        end
%         subplot(2, 5, 1:5);
%         bar(r(corp, :), 'FaceAlpha', 0.3, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');       
    end
%     sgtitle(corpora{corp})
end

figure;
% compare TIMIT and DIMEx
bar(r', 'grouped', 'FaceAlpha', 0.6, 'EdgeColor', 'none');
for corp = 1:2
    for i = 1:5
        text(i+0.2*(corp-1), r(corp, i)+0.1, ...
            getSigStr(p(corp, i), 2), 'FontSize', 13);
    end
end

ylim([-0.1 0.6]);
box off;
set(gca, 'FontSize', 15, 'Xtick', 1:5, 'Xticklabel', ...
    featureOrd(1:5), 'Ytick', [0 0.5]);
ylabel('Spearman corr')
legend({'Spanish speech', 'English speech'})

%

figure;
% compare TIMIT and DIMEx
bar(r, 'FaceAlpha', 0.6, 'EdgeColor', 'none');
for corp = 1:2
    for i = 1:5
        text((corp-1)+0.2*(i-1), r(corp, i)+0.1, ...
            getSigStr(p(corp, i), 2), 'FontSize', 13);
    end
end

ylim([-0.1 0.6]);
box off;
set(gca, 'FontSize', 15, 'Ytick', [0 0.5]);
ylabel('Spearman corr');
xticklabels({'Spanish', 'English'});
xlabel('Language')
legend(featureOrd)


%% comparing decoding weights to encoding uv v2

corpora = {'DIMEX', 'TIMIT'};
fields = {'sp_uv_all', 'eng_uv_all'};
colors = [0, 0, 1; 1, 0, 0];
r = nan(2, 5);
p = nan(2, 5);

figure;
for corp = 1:2
    corpus = corpora{corp};
    filename =  [corpus '_word_decode_monolingual_600ms.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
    for ls = corp % only show the case where its native
        weights = squeeze(mean(decode_details.weights{ls}, 1));
        avg_weights = mean(abs(decode_details.weights{ls}), [1, 3]);
        [~, max_tp] = max(avg_weights);

        % extract decoding weights
        tbl = decode_details.encoding(decode_details.elecs{ls}, :);
        tbl.mdl_weight = mean(abs(weights(max_tp-5:min(max_tp+5, 58), :)))';

        % remove electrodes if they do not contribute to the model (?)
        encode_tbl = innerjoin(tbl, wordsurp_encoding, 'Keys', {'SID', 'el'});
        featureOrd = {'onset', 'peakrate', 'formant', 'consonant', ...
            'word+surp', 'word', 'surp'};
        
        ctr=1;
        for index = [1:4 7 6]
            idx = encode_tbl.(fields{corp})(:, index)>0;
            x = encode_tbl.mdl_weight(idx);
            y = encode_tbl.(fields{corp})(idx, index);
        
            % correlate the encoding table unique variance and the decoding
            % model weight
            [r(corp, index), p(corp, index)] = corr(x, y, 'Type', 'Spearman'); 

            subplot(2, 3, ctr);

            scatter(x, y, 10, colors(corp, :), 'filled', ...
                'MarkerFaceAlpha', 0.3); hold on;
            
            % plot least squares line
            coef1 = polyfit(x,y,1);
            y1 = polyval(coef1,x);           
            
            if p(corp, index)<0.1
                plot(x, y1, 'Color', colors(corp, :), 'LineWidth', 2);
                text(0.01, 0.01, getSigStr(p(corp, index), 2)) 
            else
                plot(x, y1, 'Color',[0.4 0.4 0.4], 'LineWidth', 1.5, ...
                    'LineStyle',':');
            end
%             yticks([]);
%             xticks([]);
            xlabel('model weight');
            ylabel('\Delta R^2 ');
            title(featureOrd{index})
            ctr=ctr+1;
        end     
    end
end

figure;
% compare TIMIT and DIMEx
bar(r', 'grouped', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
for corp = 1:2
    for i = 1:5
        text(i+0.2*(corp-1), r(corp, i)+0.1, ...
            getSigStr(p(corp, i), 2), 'FontSize', 13);
    end
end

ylim([-0.1 0.6]);
box off;
set(gca, 'FontSize', 15, 'Xtick', 1:5, 'Xticklabel', ...
    featureOrd(1:5), 'Ytick', [0 0.5]);
ylabel('Spearman corr')
legend({'Spanish speech', 'English speech'})

%% comparing unique variance across features

fields = {'sp_uv_all', 'eng_uv_all'};
labels = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
numlabels = length(labels);
ls = [1, 2];

figure;
for l = ls
    % correlation matrices
    bil_corr = nan(numlabels, numlabels);
    bil_pval = nan(numlabels, numlabels);
    for i = 1:numlabels
        for j = 1:numlabels
            x = wordsurp_encoding.(fields{l})(ismember(wordsurp_encoding.ls, l), i);
            y = wordsurp_encoding.(fields{l})(ismember(wordsurp_encoding.ls, l), j);
            idx = x>0 & y>0;
            [bil_corr(i, j), bil_pval(i, j)]  = corr(x(idx), y(idx), ...
                'Type', 'Spearman');
        end
    end    
    
%     subplot(1, 2, 1);
%     imagesc(mono_corr);
%     xticks(1:5);
%     yticks(1:5);
%     xticklabels(labels)
%     yticklabels(labels)
%     caxis([-1 1]);
%     colormap(balanced);
%     
    subplot(1, 2, l);
    imagesc(bil_corr);
    for i = 1:5
        for j = 1:5
            if bil_pval(i, j)>0.05
                text(i-0.2, j, num2str(bil_corr(i, j), 2), 'FontSize', 15);
            else
                text(i-0.2, j, num2str(bil_corr(i, j), 2), ...
                    'FontSize', 15, 'FontWeight', 'bold', 'Color', 'w');
            end
        end
    end
    xticks(1:5);
    yticks(1:5);
    xticklabels(labels)
    yticklabels(labels)
    caxis([0 1]);
    colormap(flipud(gray));
    set(gca, 'FontSize', 15);
    box off;
end


%% ------------------------- single subjects ------------------------------
%% Logistic on neural data for single subjects

Swrd = TDwrd;

% expanding window
% timing = [2, 5, 10, 20];
% startp = repmat(51, 4, 1);
% timelabel = '20-200ms';

% sliding window
% startp = 50:5:90;`
% timing = repmat(4, 1, length(startp));
% timelabel = '50ms-sliding';

% % single window
% startp = 51;
% timing = 47;
% timelabel = '500ms';

% single window before onset
startp = 31;
timing = 57;
timelabel = '600ms';

% sliding window
% startp = 30:5:90;
% timing = ones(length(startp), 1)*5;
% timelabel = '5ms-sliding';

% subset to only the monolinguals
subj = 'monolingual';
mono_encoding = word_encoding(ismember(word_encoding.ls, [1, 2]), :);

tps = startp(1):startp(end)+max(timing);
nreps = 20;

uSIDs = [unique(mono_encoding.SID)];
tic
varnames = {'SID', 'elecs', 'ls', 'acc', 'weights', 'auc'}; 
decode_tbl = array2table(zeros(0, 6), 'VariableNames', varnames);

for i = 1:length(uSIDs)
    disp(['loading subject ' num2str(i) ' of ' num2str(length(uSIDs))])
    SID = uSIDs{i};
    
    elidx = find(strcmp(mono_encoding.SID, SID));
    els = mono_encoding.el(elidx);
    ls = mono_encoding.ls(elidx(1));

    % find all relevant neural responses
    tmpidx = cellfun(@(x) ~isempty(x), [Swrd.(SID)]);
    X_tmp = nan(size(Swrd.(SID){find(tmpidx, 1)}, 1), length(tps), ...
        height(Swrd));
    
    % include all non-nan responses
    tmp = cat(3, Swrd.(SID){:});
    X_tmp(:, :, tmpidx) = tmp(:, tps, :);
    
    % remove all sentence onset responses
    idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));
    X_neural = X_tmp(els, :, idx);
    clear X_tmp

    AUC = nan(length(timing), nreps);
    acc = nan(nreps, length(timing));
    scores = nan(sum(idx), length(timing));
    weights = cell(1, length(timing));
    elecs = [];
    for t = 1:length(timing)
        start = startp(t)-startp(1)+1;
        tps = start:start+timing(t);

        X = X_neural(:, tps, :);
        mintrl = 1905;

        y = Swrd.wordOns(idx)>0;
        [X_mod, nanrow, nancol, y] = makeDataMatrix(X, y, ...
            ones(length(els), 1), mintrl);
    
        if ~isempty(X_mod)
            % find true electrodes retained after de-nanning
            elshape = reshape(repmat(elidx, 1, length(tps)), ...
                length(elidx)*length(tps), []);
            elshape(nanrow) = [];
            elecs = {unique(elshape)};
    
            disp(['-------------------' SID '----------------------'])
            disp(['Electrodes in ls ' num2str(ls) ...
                ', t' num2str(timing(t)) ': ' num2str(size(X_mod, 1)/length(tps))])           
            disp(['trials: ' num2str(size(X_mod, 2))]);
            
            tic
            % computes bootstrapped accuracy across 20 repeats of 80-20 train
            % test splits, with pca computed
            % weights are significant model components x electrode
            [fp, tp, AUC(t, :), ~, tmp_scores, acc(:, t), weights{t}] = ...
                logistic(X_mod', y, 1, [], tps);
            toc
            scores(~nancol) = tmp_scores;
            disp(['mean AUC: ' num2str(mean(AUC(t, :)))])
        else
            warning([SID ' missing']);
        end
    end

    if ~isempty(X_mod)
        tmp = table({SID}, elecs, ls, {acc}, {weights}, {AUC}, 'VariableNames', ...
                varnames);
        decode_tbl = [decode_tbl; tmp];
    end
    clear tmp elecs ls acc weights AUC X_neural
end
toc
 
decode_details = struct();
decode_details.timing = timing;
decode_details.startp = startp;
decode_details.tbl = decode_tbl;
decode_details.encoding = word_encoding; % for threshold 

if startsWith(Swrd.sentid{1}, 's')
    corpus = 'DIMEX';
else
    corpus = 'TIMIT';
end
    
filename = [corpus '_word_decode_' subj '_' timelabel '_bysubj.mat']; % dimex filename
save([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: decoding accuracy across corpora per subject

timelabel = '600ms';
subj = 'bilingual'; % monolingial

% for bilinguals, display proficiency
prof.DIMEX = {5, 5, 5, 3, 3, 5, 5, 5}';
prof.TIMIT = {5, 5, 0, 5, 5, 4, 4, 3}'; 

for corp = {'TIMIT', 'DIMEX'}
    corpus = corp{1};
    load([datapath ...
        'ecog_decode/wordOnset/' corpus '_word_decode_' subj '_' timelabel '_bysubj.mat'], ...
        'decode_details');
    figure;
    [v, c] = sort(arrayfun(@(x) median(decode_details.tbl.auc(x, :), 2), ...
        1:height(decode_details.tbl)));
    tmp = [decode_details.tbl.auc]';

    if strcmpi(subj, 'bilingual')
        % get color map
        spec = spectral(8);
        spec = [0 0 0; spec([1, 1, 8, 6, 3], :)];

        ctr=1;
        for i = c
            boxplot(tmp(:, i), decode_details.tbl.SID(i), 'PlotStyle','compact', ...
            'Colors', spec(prof.(corpus){i}+1, :  ), 'Position', ctr); hold on;
            ctr=ctr+1;
        end
    else
        h = boxplot(tmp(:, c), decode_details.tbl.SID(c), 'ColorGroup', ...
            decode_details.tbl.ls(c), 'PlotStyle','compact', ...
            'Colors', getColorsCrossComp(1), 'HandleVisibility', 'on');
    end

    numel = arrayfun(@(x) length(decode_details.tbl.elecs{x}), ...
        1:height(decode_details.tbl));
    
    ctr=1;
    for i = c
        text(ctr-0.05, median(decode_details.tbl.auc(i, :), 2)+0.08, ...
            num2str(numel(i)));
        
%         if strcmpi(subj, 'bilingual')
%             text(ctr-0.05, median(decode_details.tbl.auc(i, :), 2)+0.1, ...
%                 num2str(prof.(corpus){i}))
%         end
        ctr=ctr+1;
    end
    yline(0.5)
    ylim([0.4 0.78]);
    yticks(0.4:0.1:0.8);
    box off;
    set(gca, 'FontSize', 15);
    ylabel('AUC');
    xlabel('Subject');
    sgtitle(corpus)
end

tmp = load([datapath ...
        'ecog_decode/wordOnset/DIMEX_word_decode_' subj '_' timelabel '_bysubj.mat'], ...
        'decode_details');
dimex=tmp.decode_details.tbl;
tmp = load([datapath ...
        'ecog_decode/wordOnset/TIMIT_word_decode_' subj '_' timelabel '_bysubj.mat'], ...
        'decode_details');
timit=tmp.decode_details.tbl;

joined = join(timit, dimex, 'Keys', 'SID');
cols = getColorsCrossComp(1);
figure;
for i = 1:height(joined)
    [h, p] = ttest2(joined.auc_timit(i, :),  joined.auc_dimex(i, :));

    if p<0.01
        style = "-";
    elseif p<0.05
        style = "--";
    else
        style = ":";
    end
    disp([joined.SID(i) num2str(p)])

    if joined.ls_dimex(i)==1
        subplot(1, 2, 1);
%         boxplot([joined.auc_timit(i, :); joined.auc_dimex(i, :)]', 'PlotStyle','compact'); 
        
        plot([1, 2], median([joined.auc_timit(i, :); joined.auc_dimex(i, :)], 2), ...
            'o-', 'LineWidth', 1.5, 'Color', cols(joined.ls_dimex(i), :), 'LineStyle', style);
        text(2, median(joined.auc_dimex(i, :)), joined.SID(i));
        title('Spanish subjects')
        ylim([0.45 0.75])
        hold on;
    else
        subplot(1, 2, 2);
%         boxplot([joined.auc_timit(i, :); joined.auc_dimex(i, :)]', 'PlotStyle','compact');
        plot([1, 2], median([joined.auc_timit(i, :); joined.auc_dimex(i, :)], 2), 'o-', ...
                'LineWidth', 1.5, 'Color', cols(joined.ls_dimex(i), :), 'LineStyle', style);
        text(2, median(joined.auc_dimex(i, :)), joined.SID(i));
        title('English subjects')
        hold on;
        ylim([0.45 0.75])
    end
    ylabel('AUC');
    xlabel('Language presented');
    xticks([1, 2]);
    xticklabels({'English', 'Spanish'});
    box off;
    set(gca, 'FontSize', 15);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;



%% vis: decoding weights across corpora per subject

timelabel = '600ms';
a = load([datapath ...
    'ecog_decode/wordOnset/TIMIT_word_decode_monolingual_' timelabel '_bysubj.mat'], ...
    'decode_details');
timit_word = a.decode_details.tbl;

b = load([datapath ...
    'ecog_decode/wordOnset/DIMEX_word_decode_monolingual_' timelabel '_bysubj.mat'], ...
    'decode_details');
dimex_word = b.decode_details.tbl;

clear a b

% SIDs = {'EC100', 'EC172', 'EC183', 'EC105', 'EC163', 'EC252', 'EC222'};
SIDs = unique(timit_word.SID);
for s = SIDs'
    SID = s{1};

    didx = strcmp(dimex_word.SID, SID);
    tidx = strcmp(timit_word.SID, SID);
    shared = intersect(dimex_word.elecs{didx},timit_word.elecs{tidx});
    dels = find(ismember(dimex_word.elecs{didx}, shared));
    tels = find(ismember(timit_word.elecs{tidx}, shared));
    
    figure;
    subplot(1, 3, 1)
    
    % sum of absolute weights in time
    dimex_weights = sum(abs(squeeze(mean(dimex_word.weights{didx}(:, :, dels)))));
    
    % sum of absolute weights in time
    timit_weights = sum(abs(squeeze(mean(timit_word.weights{tidx}(:, :, tels)))));
    
    scatter(dimex_weights, timit_weights, 40, 'filled', 'k');
    lsline
    [rho, pval] =corr(dimex_weights', timit_weights', 'Type', 'Spearman'); 
    ylabel('timit abs weight'); 
    xlabel('dimex abs weight');
    title([num2str(rho, 2) ', ' num2str(pval, 2)]);
    set(gca, 'FontSize', 15);

    subplot(1, 3, 2)
    [~, idx] = sort(squeeze(sum(abs(mean(dimex_word.weights{didx}(:, :, dels))))), 'descend');
    imagesc(squeeze(mean(dimex_word.weights{didx}(:, :, dels(idx))))');
    title('DIMEX weights');
    ylabel('Electrodes'); 
    xlabel('Time (10ms)'); 
    colormap(vik); 
    set(gca, 'FontSize', 15);
    
    subplot(1, 3, 3)
    [~, idx] = sort(squeeze(sum(abs(mean(timit_word.weights{tidx}(:, :, tels))))), 'descend');
    imagesc(squeeze(mean(timit_word.weights{tidx}(:, :, tels(idx))))');
    title('TIMIT weights');
    ylabel('Electrodes'); 
    xlabel('Time (10ms)'); 
    colormap(vik); 
    set(gca, 'FontSize', 15);

    sgtitle([SID ', d' num2str(median(dimex_word.auc(didx, :))) ':t' ...
        num2str(median(timit_word.auc(tidx, :)))]);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;


%% FIX: correlation of weights to word unique variance, peak rate unique
% variance, phonetic feature unique variance, and base rsq
%figure;
data = cell(2, 1); % comparing weights of model to specific continuous ...
% unique variance variables
for ls = 1:2
    for t = 1:4
        %subplot(2, 4, (ls-1)*4+t)

        % remove all components that were not weighted significantly
        % differently from 0
        weights = abs(comp_weighted{ls, t});
        weights(weights==0)=NaN;
        weights(sum(isnan(weights), 2)>0, :)=[];

        data{ls}(t, :) = sum(weights)';

%         imagesc(weights)
%         scatter(weights(1, :), mono_encoding.uv_wordOns(elecs{ls, t}, 1), ...
%             35, 'filled'); hold on;
%         lsline;
    end
    % change to match the decoding language
    data{ls}(t+1, :) = mono_encoding.uv_wordOns(elecs{ls, t}, 1);
    data{ls}(t+2, :) = mono_encoding.uv_peakRate(elecs{ls, t}, 1);
end

figure; 
varnames = [split(num2str(timing*10)); {'word \Delta R^2'; 'peakrate \Delta R^2'}];
[~,PValue] = corrplot(data{1}', 'Varnames', varnames);
arrayfun(@(x) delete(subplot(6, 6, x)), find(PValue>0.05 & PValue<1));

figure; 
[~,PValue] = corrplot(data{2}', 'Varnames', varnames);
arrayfun(@(x) delete(subplot(6, 6, x)), find(PValue>0.05 & PValue<1));


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* X_neural;



%% FIX:
corpus = 'TIMIT';

labels = {'word + surprise', 'peak rate', 'phonetic features'};
subjs = {'Spanish', 'English','Mandarin', 'Bilingual'};
colors = getColorsCrossComp(1);

f = figure;
featuvs = cell(4, 2);
for ls = [1:2 4]

    if ismember(ls, [1, 2])
        subj = 'monolingual';
    else
        subj = 'bilingual';
    end

    load([datapath ...
        'ecog_decode/wordOnset/' corpus '_word_decode_' subj '_500ms.mat'], ...
        'decode_details');
    % find highest weighted electrode
    weights = decode_details.weights{ls};
    
    p = nan(1, size(weights, 3));
    for e = 1:size(weights, 3)
        [~, p(e)] = ttest(squeeze(mean(weights(:, :, e), 2)));
    end
    sigw = p<0.001; 

    weights = squeeze(sum(abs(mean(weights(:, :, sigw))), 2));
    els = decode_details.encoding.el(decode_details.elecs{ls}(sigw));
    sids = decode_details.encoding.SID(decode_details.elecs{ls}(sigw));

    if strcmp(corpus, 'TIMIT')
        featuv  = arrayfun(@(x)  word_encoding.uv_timit(strcmp(word_encoding.SID, sids(x)) & ...
                word_encoding.el==els(x)), 1:length(els));
    else
        featuv  = arrayfun(@(x)  word_encoding.uv_dimex(strcmp(word_encoding.SID, sids(x)) & ...
                word_encoding.el==els(x)), 1:length(els));
    end

    for feat = 1:3
        set(0, 'CurrentFigure', f);
%         subplot(1, 3, feat);        
        scatter(rescale(weights), cellfun(@(x) x(feat), featuv), ...
            35, colors(ls, :), 'filled', 'MarkerFaceAlpha', 0.3);
        [r, pv] = corr(weights, cellfun(@(x) x(feat), featuv)', 'type', 'Spearman');
        hold on;
        disp(['ls ' num2str(ls) ': ' num2str(pv) ', ' num2str(r)]);
        ylim([-0.005 0.05]);
        yticks(0:0.02:0.05);
        xlim([0 1]);
        ylabel('\Delta R^2');

        set(gca, 'FontSize', 15);
        title(labels(feat))
        %set(gca, 'YScale', 'log');
    end

    featuv = cat(1, featuv{:});
    featuv(featuv<0) = 0;

    idx = find(weights>prctile(weights, 50));
    [~, sorted] = sort(weights(idx), 'descend');
%     figure; bar(featuv(idx(sorted), 1:3)./ sum(featuv(idx(sorted), 1:3), 2));

    eva = evalclusters(featuv(:, 1:3), 'kmeans', 'CalinskiHarabasz', 'klist', 1:3);
    
    numel = 3;
    widx = idx(sorted);
    cols = crest(numel);
    i = 1; 
    featuvs{ls, 1} = featuv(:, 1:3);
    featuvs{ls, 2} = weights;

    figure; 
    p = pie(mean(featuv(:, 1:3).*weights, 1)*2000, [1, 1, 1]);
    for j = 1:3
        p((j-1)*2+1).FaceColor = cols(j, :); 
        p((j-1)*2+1).EdgeColor = 'none';
    end 
    legend({'word + surprise', 'peakRate', 'phonetic feature'});
    sgtitle(subjs{ls});  
    set(gca, 'FontSize', 15);

%     while ctr < numel+1 && i < length(widx)
%         subplot(1, numel, ctr); 
%              
%         % show pie charts for electrodes with more than one feature explaining variance
%         if sum(featuv(widx(i), 1:3)>0)>1
%             p = pie(featuv(widx(i), 1:3)./sum(featuv(widx(i), 1:3)), ...
%                 [1, 1, 1], {'', '', ''}); 
% 
%             for j = 1:3
%                 p((j-1)*2+1).FaceColor = cols(j, :); 
%                 p((j-1)*2+1).EdgeColor = 'none';
%             end            
% %             subplot(2, numel, ctr+numel)
% %             imagesc(squeeze(decode_details.weights{ls}(:, :, widx(i))));   
% %             colormap(vik);
%             ctr = ctr+1;            
%         end 
%         i = i+1;
% 
%         if ctr==numel
%             legend({'word + surprise', 'peakRate', 'phonetic feature'});
%             set(gca, 'FontSize', 15);     
%         end
%         sgtitle(subjs{ls});   
% 
% 
% 
%     end
    % show unique variance plots for dimex and timit
%     figure;
%     subplot(1, 2, 1);
%     pie(word_encoding.uv_timit{el}*100, labels);
% 

%     subplot(1, 2, 2);
%     pie(word_encoding.uv_dimex{el}*100, labels);
%     sgtitle([word_encoding.SID{el} ' el: ' num2str(word_encoding.el(el)) ...
%         ' ' num2str(word_encoding.ls(el))])
%     colormap(brewermap(5, 'Set2'));
end

figure;
b = [];
for ls = [1, 2, 4]
    normw = featuvs{ls, 2}./norm(featuvs{ls, 2});
    tmp = squeeze(cat(1, featuvs{ls, 1})).*normw;
    nofeat = all(tmp==0, 2);
    b = [b; mean(tmp(~nofeat, :))];

    
end
legend({'word + surprise', 'peakRate', 'phonetic feature'});
ba = bar(b);

cm = gray(5);
cm = cm([1 3 4], :);
for k = 1:size(b,1)
    ba(k).FaceColor = cm(k, :);
end

xticks(1:3);
xticklabels(subjs([1, 2, 4]));

%% FIX: erps word boundary
rmpath('util/shadederror');
addpath(genpath('util'))
SID = 'EC222'; % 'EC222'
el = 119; %119

% EC172, el 122 -- Dip is Spanish word onset?
% EC100, el 21, 118 -- peakrate
% EC235, el 249, 250 -- English word onset?
% EC222, el 101 -- English peakRate?

ctr = 1;
dummy = cell(1, 2);
for tmp = {Dwrd, TDwrd}
    swrd = tmp{1};
    dummy{ctr} = struct();
    nanidx = cellfun(@(x) isempty(x), swrd.(SID)) | ...
        (swrd.syll<2 & ~isnan(swrd.syll));
    dummy{ctr}.(SID).resp = cat(3, swrd.(SID){~nanidx});
    
    dummy{ctr}.wordOns = swrd.wordOns(~nanidx);
    dummy{ctr}.precWordLen = swrd.precWordLen(~nanidx);
    dummy{ctr}.syllOns = ones(sum(~nanidx), 1);
    ctr = ctr + 1;
end

% word responses
f = figure;
ax(1) = subplot(1, 2, 1);
cols = [0.6 0 0.6; 0.1 0.5 0.05;]; %  0.6 0.2 0.6
plotWordErp(dummy{2}, SID, el, 'timit', f, cols, 1, bef./100);
xlabel('Time (s)');
set(gca, 'FontSize', 15);
ax(2) = subplot(1, 2, 2);
plotWordErp(dummy{1}, SID, el, 'dimex', f, cols, 1, bef./100);
xlabel('Time (s)');
set(gca, 'FontSize', 15);
%linkaxes(ax, 'y');

% word responses by previous word length
f = figure;
ax(1) = subplot(1, 2, 1);
cols = brewermap(6, 'Spectral'); 
plotWordErp(dummy{2}, SID, el, 'timit', f, cols, 2, bef./100);
ylabel('HGA (z-score)');
xlabel('Time (s)');
set(gca, 'FontSize', 15);
ax(2) = subplot(1, 2, 2);
cols =  brewermap(6, 'Spectral');
addpath(genpath('util'))
plotWordErp(dummy{1}, SID, el, 'dimex', f, cols, 2, bef./100);
xlabel('Time (s)');
set(gca, 'FontSize', 15);
linkaxes(ax, 'y');

% word responses by vowel vs. consonant as end/beginning

% look at STRF weights for electrode

% SID = 'EC222'; % 'EC222'
% el = 87;

modelnames={'onset_aud_maxDtL_wordOns_wordL'};  
corpusStrf = loadMultModelStrf(SID, modelnames, 'timit', datapath, 1);
figure;
% imagesc for TIMIT
subplot(2, 3, [1 4]);
x =  0:-0.01:-0.6;
imagesc(squeeze(corpusStrf{1}.meanStrf(:, :, el)), 'Xdata', x);
yline(1.5, 'LineWidth', 2,'Color','k');
yline(83.5, 'LineWidth', 2,'Color','k');
xlabel('Time (-s)');
ylabel('predictor');
yticks([1 40 83]);
yticklabels({'onset', 'aud', 'word'});
ytickangle(45)
set(gca, 'FontSize', 15);

% beta weights for corpus STRFs
subplot(2, 3, 3);
plot(x, squeeze(corpusStrf{1}.meanStrf(82, :, el)), ...
    'LineWidth', 2, 'Color', 'k'); hold on;
plot(x, squeeze(corpusStrf{1}.meanStrf(83, :, el)), ...
    'LineWidth', 2, 'Color', 'r');
ylabel('beta weight');
yticks([]);
xlabel('Time (-s)');
set(gca, 'FontSize', 15);
box off;

corpusStrf = loadMultModelStrf(SID, modelnames, 'dimex', datapath, 1);
subplot(2, 3, [2, 5]);
imagesc(squeeze(corpusStrf{1}.meanStrf(:, :, el)), 'Xdata', x);
yline(1.5, 'LineWidth', 2,'Color','k');
yline(83.5, 'LineWidth', 2,'Color','k');
xlabel('Time (-s)');
ylabel('predictor');
yticks([])
set(gca, 'FontSize', 15);

subplot(2, 3, 6);
plot(x, squeeze(corpusStrf{1}.meanStrf(82, :, el)), ...
    'LineWidth', 2, 'Color', 'k'); hold on;
plot(x, squeeze(corpusStrf{1}.meanStrf(83, :, el)), ...
    'LineWidth', 2, 'Color', 'r');
legend({'peakRate', 'word'});
ylabel('beta weight');
yticks([]);

xlabel('Time (-s)');
set(gca, 'FontSize', 15);
box off;

clearvars -except *all subj *vow* *cons* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% ------------ CONTROL: NEURAL CLASSIFICATION OF SONORITY ----------------
%% Logistic on neural data for different times around onset
% requires you to load word responses for all SIDs
Swrd = Dwrd;
% english/spanish - [1, 2], mandarin - 3, bilingual - 4
ls = [4];% [1, 2]; 
subtype = ''; % using different model, aud + word/wordOnset _word
if startsWith(Swrd.sentid{1}, 's')
    corpus = 'DIMEX';
    corpus_details = dimex_details;
    column = 2;
else
    corpus = 'TIMIT';
    corpus_details = timit_details;
    column = 1;
end
    
% Swrd = Dwrd;
% corpus = 'DIMEX';

idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));

% single window
% startp = 51;
% timing = 47;
% timelabel = '500ms';

% single window before onset
startp = 31;
timing = 57;
timelabel = '600ms';

% subset to only the monolinguals
% mono_encoding = word_encoding(word_encoding.ls<3, :);

encoding_tbl = word_encoding;
word_uv = nan(height(encoding_tbl ), 1);
for i = 1:height(encoding_tbl )
    if ~isempty(encoding_tbl .(['uv_' lower(corpus)]){i})
        % true word unique variance versus word+surprise (1 or 4)
        word_uv(i) = encoding_tbl.(['uv_' lower(corpus)]){i}(4);
    end
end
wordAud_uv = encoding_tbl.uv_word(:, column);

switch subtype
% was using surprisal+word model
    case '_acs' % only "acoustic" electrodes
        tblidx = ismember(encoding_tbl.ls,ls) & (word_uv<0 & ~isnan(word_uv));
        mono_encoding = encoding_tbl (tblidx, :);
    case '_acsTop'
        thresh = 0.015;
        % threshold on acoustic feature rsq
        tblidx = ismember(encoding_tbl.ls,ls) & ...
            (word_uv<0 & ~isnan(word_uv)) & ...
            encoding_tbl.(['acsfeat_rsq_' lower(corpus)])>thresh;
        mono_encoding = encoding_tbl (tblidx, :);
    case '' % all
        mono_encoding = word_encoding(ismember(word_encoding.ls,ls), :);
end
clear word_uv

tps = startp(1):startp(end)+max(timing);
X_neural = nan(height(mono_encoding), length(tps), sum(idx));

% aggregation of neural data
disp('Loading all neural data....');
uSIDs = [unique(mono_encoding.SID)];
tic
for i = 1:length(uSIDs)
    disp(['loading subject ' num2str(i) ' of ' num2str(length(uSIDs))])
    SID = uSIDs{i};
    elidx = find(strcmp(mono_encoding.SID, SID));
    els = mono_encoding.el(elidx);

    % find all relevant neural responses
    tmpidx = cellfun(@(x) ~isempty(x), [Swrd.(SID)]);
    X_tmp = nan(size(Swrd.(SID){find(tmpidx, 1)}, 1), length(tps), ...
        height(Swrd));
    
    % include all non-nan responses
    tmp = cat(3, Swrd.(SID){:});
    X_tmp(:, :, tmpidx) = tmp(:, tps, :);
    % remove all sentence onset responses

    % make sure all electrodes are under maximum electrode number
    els = els(els<size(X_tmp, 1));
    elidx = elidx(els<size(X_tmp, 1));

    X_neural(elidx, :, :) = X_tmp(els, :, idx);
    clear X_tmp
end
toc
disp('Loading neural data complete!');

% logistic regression with 5-fold cross-validation
idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));
if startsWith(Swrd.sentid{1}, 's')
    cols = brewermap(length(startp)+1, 'YlGnBu'); % dimex color map
else
    cols = brewermap(length(startp)+1, 'YlOrRd'); % timit color map
end

cols = cols(2:end, :);
nfolds = 20;
subjs = getSubjectLabels(ls);
if ls == 3, nfolds = 5; end

AUC = nan(length(subjs), length(timing), nfolds);
acc = cell(length(subjs), 1);
% electrode weighting through the logistic regression and PCA 
comp_weighted = cell(length(subjs), length(timing));
elecs = cell(length(subjs), length(timing));
scores = cell(length(subjs), length(timing));
yidx = cell(length(subjs), 1);

f = figure;
for ls = ls % [1 , 2]
%     subplot(1, 2, ls)
    for t = 1:length(timing)
        % assumes that X_neural starts at startp(1)
        start = startp(t)-startp(1)+1;
        tps = start:start+timing(t);

        elidx = mono_encoding.ls==ls;
        elshape = reshape(repmat(find(elidx), 1, length(tps)), ...
        sum(elidx)*length(tps), []);

        X = X_neural(elidx, tps, :);        
    
        % for onset phone sonority decoding
        onsPhns = arrayfun(@(x) corpus_details.phnnames(x), [Swrd.onsPhn]);
        onsSonor = cellfun(@(x) ismember(x, corpus_details.features.sonorant), onsPhns);
        y = onsSonor;

        if ls == 3
            mintrl = 500;
        else
            mintrl = 2000;
        end            

        [X, nanrow, nancol, y] = makeDataMatrix(X, y, ...
            mono_encoding.SID(elidx), mintrl);
        tmp = find(idx);
        yidx{ls} = tmp(~nancol);
        elshape(nanrow) = [];
        elecs(ls, t) = {unique(elshape)};

        % how many electrodes and subjects are included in the final matrix
        disp(['Electrodes in ls ' num2str(ls) ...
            ', t' num2str(timing(t)) ': ' num2str(size(X, 1)/length(tps))])
        disp(['Unique subjects included: ' ...
            num2str(length(unique(mono_encoding.SID(elecs{ls, t}))))...
            ', trials: ' num2str(size(X, 2))])
 
        % computes bootstrapped accuracy across 50 repeats of 80-20 train
        % test splits, with pca computed
        % weights are significant model components x electrode
        [fp, tp, AUC(ls, t, :), ~, tmp_score, acc{ls}(:, t), weights] = ...
            logistic(X', y, 1, [], tps, nfolds);
        scores{ls, t} = nan(sum(idx), 1);
        scores{ls, t}(~nancol) = tmp_score;
        comp_weighted(ls, t) = {weights};

        % ROC curve
        set(0, 'CurrentFigure', f)
        plot(fp, tp, 'LineWidth', 2, 'Color', cols(t, :)); hold on;   
    end

    h = refline(1, 0);
    h.LineStyle = '--';
    h.Color = 'k';
    xlabel('False positive rate') 
    ylabel('True positive rate')
    
    box off; set(gca, 'FontSize', 15);  
    %title(['Subject Group: ' subjs{ls}])
    %legend(split(num2str(timing*10)));
end

decode_details = struct();
decode_details.acc = acc; % accuracy for both subject group types
decode_details.timing = timing;
decode_details.startp = startp;
decode_details.weights = comp_weighted;
decode_details.AUC = AUC;
decode_details.elecs = elecs;
decode_details.encoding = mono_encoding; % for threshold 
decode_details.score = scores;
decode_details.y = Swrd.wordOns(idx)>0;
decode_details.yidx = yidx;

if isscalar(ls) && ls == 3
    subj = 'mandarinmono';
elseif isscalar(ls) && ls == 4
    subj = 'bilingual';
else
    subj = 'monolingual';
end

filename = [corpus '_onsSonority_decode_' subj '_' timelabel subtype '.mat']; % dimex filename
save([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd* X_neural;



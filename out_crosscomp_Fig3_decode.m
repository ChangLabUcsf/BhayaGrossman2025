% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;
SIDs = [sSIDs eSIDs bSIDs];

% corpus details
% tps = 50:55;

% selected electrodes
timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

% before and after word time points
bef=50;
aft=50;
%
% loading in word-level subject data % eSIDs, sSIDs
TDwrd = loadDwrdMulti('timit', bef, aft, [eSIDs sSIDs], timit_details);
Dwrd = loadDwrdMulti('dimex',  bef, aft, [eSIDs sSIDs], dimex_details);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% Load in the electrodes to use (all speech responsive
% 
% % load in TRF models   
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

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% Running the logistic classification on neural data for different times around onset

% requires you to load word responses for all SIDs [acoustic decoding]
Swrd = Dwrd;
window = 20; % 200ms window around boundary
Swrd.ambiguity = getAAI(window, Swrd, 'aud');

% english/spanish - [1, 2], mandarin - 3, bilingual - 4
ls = [1, 2]; %[1, 2];% [1, 2]; 
subtype = ''; % using different model, aud + word/wordOnset _word
trialset = ''; % all trials, high ambiguity or low ambiguity _highAAI
maxsid = 3; % maximum number of subjects to include

if startsWith(Swrd.sentid{1}, 's')
    corpus = 'DIMEX';
    lowAAIthresh = prctile(Swrd.ambiguity, 20);
    highAAIthresh = prctile(Swrd.ambiguity, 80);
    column = 2;
else
    corpus = 'TIMIT';
    lowAAIthresh = prctile(Swrd.ambiguity, 20);
    highAAIthresh = prctile(Swrd.ambiguity, 80);
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
% startp = 31;
% timing = 57;
% timelabel = '600ms';

% sliding window
startp = 30:5:90;
timing = ones(length(startp), 1)*10;
timelabel = '10ms-sliding';

% sliding window
% startp = 30:5:90;
% timing = ones(length(startp), 1)*5;
% timelabel = '5ms-sliding';

% subset to only the monolinguals
% mono_encoding = word_encoding(word_encoding.ls<3, :);

encoding_tbl = word_encoding;
% for when you want to use only elecs with uv for word features
% word_uv = nan(height(encoding_tbl ), 1);
% for i = 1:height(encoding_tbl )
%     if ~isempty(encoding_tbl .(['uv_' lower(corpus)]){i})
%         % true word unique variance versus word+surprise (1 or 4)
%         word_uv(i) = encoding_tbl.(['uv_' lower(corpus)]){i}(4);
%     end
% end
% wordAud_uv = encoding_tbl.uv_word(:, column);

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
yidx = cell(length(subjs), length(timing));
ys = cell(length(subjs), length(timing));
Mdls = cell(length(subjs), length(timing));
pcaXs = cell(length(subjs), length(timing));

f = figure;
for ls = ls 

    for t = 1:length(timing)
        % reset idx for each subject group
        idx = Swrd.sentOns<1 & (Swrd.syll>1 | isnan(Swrd.syll));
        % subset to both trialsets to ensure same subjects etc are used across sets
        if any(strcmp(trialset, {'_highAAI', '_lowAAI'}))
            idx = idx & (Swrd.ambiguity>highAAIthresh ...
                | Swrd.ambiguity<lowAAIthresh);
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
            mono_encoding.SID(elidx), mintrl, maxsid, 1);
        tmp = find(idx);
        yidx{ls, t} = tmp(~nancol); % this won't be correct in AAI cases
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
        ys{ls, t} = y;

         % Display how many electrodes and subjects are included in the final matrix
        disp(['Electrodes in ls ' num2str(ls) ...
            ', t' num2str(timing(t)) ': ' num2str(size(X, 1)/length(tps))])
        disp(['Unique subjects included: ' ...
            num2str(length(unique(mono_encoding.SID(elecs{ls, t}))))...
            ', trials: ' num2str(size(X, 2))])
 
        % Computes bootstrapped accuracy across 50 repeats of 80-20 train
        % test splits, with pca computed
        % Weights are significant model components x electrode
        [fp, tp, AUC(ls, t, :), pcaX, tmp_score, acc{ls}(:, t), weights, ...
            ~, mdl] = logistic(X', y, 1, [], tps, nfolds);
        pcaXs(ls, t) = {pcaX}; % for running the mdl on different AAI bins
        scores{ls, t} = nan(sum(idx), 1);
        Mdls(ls, t) = {mdl};
        clear mdl
        if ~strcmp(trialset, '_highAAI') && ~strcmp(trialset, '_lowAAI')
            assert(length(tmp_score) == length(y));
            %scores{ls, t}(~nancol) = tmp_score;
            scores{ls, t} = tmp_score;
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
% should match the scores
% decode_details.y = Swrd.wordOns(idx)>0;
decode_details.ys = ys;
decode_details.yidx = yidx;
decode_details.Mdl = Mdls;
decode_details.pcaX = pcaXs;

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
varnames = {'SID', 'elecs', 'ls', 'numtrials', 'acc', 'weights', 'auc'}; 
decode_tbl = array2table(zeros(0, 7), 'VariableNames', varnames);

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
        mintrl = 1900;

        y = Swrd.wordOns(idx)>0;
        [X_mod, nanrow, nancol, y] = makeDataMatrix(X, y, ...
            ones(length(els), 1), mintrl, 1, 1);
    
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
        tmp = table({SID}, elecs, ls, size(X_mod, 2), {acc}, ...
            {weights}, {AUC}, 'VariableNames', ...
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


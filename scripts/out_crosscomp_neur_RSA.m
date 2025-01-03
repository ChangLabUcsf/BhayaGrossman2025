
% Ilina Bhaya-Grossman
% 08/02/2023
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('loadData'))

addpath brewer
addpath(genpath('util'))
zFolder = 'block_z'; % 'block_z'
[datapath, ~] = setDatapath;
addpath(genpath(datapath))

bef=10;
aft=50;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
% sSIDs = [sSIDs, {'EC225', 'EC152'}];
% eSIDs = [eSIDs, {'EC242', 'EC142'}];
% SIDs = [sSIDs, eSIDs, {'HS11', 'HS9', 'HS10'}];
SIDs = [sSIDs];
subjectgroup = 'spanishmono';
phngroup = 'all'; % sonorant, obstruent, all 

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
% all the spectral slices to correlate with neural RDMs
tps = [50:60; 60:70; 70:80; 80:90];

exclude_phns = {'b_c', 'd_c', 'g_c', 'p_c', 't_c', 'tS_c', 'k_c', 'bcl', 'dcl', ...
'gcl', 'pcl', 'tcl', 'kcl', 'epi', 'h#', 'pau', 'epi', 'ax-h', ...
'-B','-D', 'Z','n~','tS', ... % Spanish phones with less than 100 trials
'ch', 'el', 'em', 'en', 'eng', 'jh', 'nx', 'oy', 'th', 'uh', 'uw', 'y', 'zh'}; % English phones

% exclude sonorant or obstruent phones from the RDMs
if strcmp(phngroup, 'obstruent')
    exclude_phns = [exclude_phns, timit_details.features.sonorant, dimex_details.features.sonorant];
elseif strcmp(phngroup, 'sonorant')
    exclude_phns = [exclude_phns, timit_details.features.obstruent, dimex_details.features.obstruent];
end

% exclude closure phns
include_phns = dimex_details.phnnames(~ismember(dimex_details.phnnames, exclude_phns));
Dcons = loadDDcons('dimex', bef, aft, SIDs, [], [], ...
    dimex_details.sentdet, include_phns);
Dcons.name = 'dimex';

include_phns = timit_details.phnnames(~ismember(timit_details.phnnames, exclude_phns));
TDcons = loadDDcons('timit', bef, aft, SIDs, [], [], ...
    timit_details.sentdet, include_phns);
TDcons.name = 'timit';

clearvars -except Dcons TDcons dimex_details timit_details tps SIDs subjectgroup ...
    datapath dpath zFolder bef aft exclude_phns RDMs neural_RDMs nancol phngroup

%% -------------------------- Set RDM details --------------------------

% remove all phones for which there are less than 100 trials
% already removed in loadDDcons
% [Dcons, removedPhns] = removePhns(Dcons, 100);
% [TDcons, removedPhns] = removePhns(TDcons, 100);

% randomly select 3000 phoneme trials from each corpora
rng(1);
idx = randperm(size(Dcons.aud, 3), 1000);
Dcons.randidx = idx;

idx = randperm(size(TDcons.aud, 3), 1000);
TDcons.randidx = idx;

% convert phns to phonetic features
Dcons = phn2feat(Dcons, dimex_details);
TDcons = phn2feat(TDcons, timit_details);

clearvars -except Dcons TDcons dimex_details timit_details tps neurtps SIDs subjectgroup ...
    datapath dpath zFolder bef aft exclude_phns RDMs neural_RDMs nancol phngroup

%% -------------------------- Generate acoustic RDMs --------------------------

%% RDMs from spectrograms, phoneme label, and phonetic features

Scons = {Dcons, TDcons};
% corpus x (spectrogram, phoneme label, phonetic features)
RDMs = cell(length(Scons), size(tps, 1)+2);

for s = 1:length(Scons)
    Scon = Scons{s};

    % sort trial by phoneme label 
    [~, ~, phnidx] = unique(Scon.phn(Scon.randidx));

    for tprow = 1:size(tps, 1) 
        % use preselected trials
        aud = Scon.aud(:, tps(tprow, :), Scon.randidx);
        
        % (freq x time x trial) reshape to (freq*time x trial)
        aud = reshape(aud, size(aud, 1)*size(aud, 2), size(aud, 3));

        % pca on spectrogram
        [~, score, ~, ~, explained, ~] = pca(aud');

        % compute RDM from the components that explain 90% of the variance
        n = find(cumsum(explained) > 90, 1);
        disp(['Number of components that explain 90% of the variance: ', num2str(n)])
        aud = score(:, 1:n);

        % compute RDM for spectrogram
        RDM = squareform(pdist(aud, 'correlation'));
        RDMs{s, tprow} = RDM;
        clear RDM
    end

    % compute RDM for phonetic features
    phnfeat = Scon.phnfeat(:, Scon.randidx);
    RDM = squareform(pdist(phnfeat', 'hamming'));
    RDMs{s, size(tps, 1)+1} = RDM;

     % compute RDM for phoneme label (1 if different, 0 if same)
    RDM = phnidx ~= phnidx';
    RDMs{s, size(tps, 1)+2} = double(RDM);
    clear RDM
end

clearvars -except Dcons TDcons dimex_details timit_details neurtps tps SIDs subjectgroup ...
    datapath dpath zFolder bef aft exclude_phns RDMs neural_RDMs nancol phngroup

%% Plot acoustic RDMs

Scons = {Dcons, TDcons};

figure;
ctr = 1;
for s = 1:length(Scons)
    Scon = Scons{s};

    % sort trial by phoneme label 
    [phn, ~, phnidx] = unique(Scon.phn(Scon.randidx));
    [phntype, sortidx] = sort(phnidx);

    for type = 1:length(RDMs(1, :))
        % plot RDMs
        subplot(2, length(RDMs(1, :)), ctr)
        imagesc(RDMs{s, type}(sortidx, sortidx))
        
        title(Scon.name);

        % plot lines on RDM to separate phonemes
        hold on;
        startphn = zeros(1, length(phn));
        numphn = zeros(1, length(phn));
        for i = 1:length(phn)
            startphn(i) = find(phntype == i, 1);
            numphn(i) = sum(phntype == i);
            xline(startphn(i), 'k', 'LineWidth', 0.8)
            yline(startphn(i), 'k', 'LineWidth', 0.8)
        end

        % set the tick labels to be the phoneme labels
        set(gca, 'FontSize', 13, 'Xtick', startphn+numphn/2, 'XtickLabel', phn, ...
            'XTickLabelRotation', 90, 'Ytick', startphn+numphn/2, 'YtickLabel', phn);
        ctr = ctr + 1;
        colormap(mako)
    end
end

clearvars -except Dcons TDcons dimex_details timit_details neurtps tps SIDs subjectgroup ...
    datapath dpath zFolder bef aft exclude_phns RDMs neural_RDMs nancol phngroup 

%% -------------------------- Generate neural RDMs --------------------------

%% RDMs from neural data
Scon = Dcons;
Scons = {Dcons};

% read in relevant neural data / load electrode set
load(['select_elec/out_elecs_speechtypeftest_bychan_' Scon.name '_all.mat'], 'allidx')

% aggregate neural data across all subjects
resp_all = [];
sids = [];
neurtps = bef+10:aft+bef;
for s = 1:length(SIDs)
    SID = SIDs{s};

    % if subject has no speech responsive sites, skip
    if ~isfield(allidx, SID)
        continue
    end

    resp_all = cat(1, resp_all, ...
        Scon.(SID).resp(allidx.(SID), neurtps, Scon.randidx));
    sids = cat(2, sids, repmat(str2double(SID(3:end)), ...
        1, length(allidx.(SID))));
end

% compute one RDM per time point
neural_RDMs = cell(1, size(resp_all, 2));

% find nan trials
nancol = zeros(size(resp_all, 2), size(resp_all, 3));
for t = 1:size(resp_all, 2)
    [X, ~, nancol(t, :)] = makeDataMatrix(resp_all(:, t, :), ...
        ones(size(resp_all, 3), 1), sids, 500);

    % compute pca on neural data
    [~, score, ~, ~, explained, ~] = pca(X');
    n = find(cumsum(explained) > 80, 1);
    disp(['Number of components that explain 80% of the variance: ', num2str(n)])
    RDM = squareform(pdist(score(:, 1:n), 'correlation'));
    % RDM = 1-corr(score(:, 1:n)');

    % compute RDM
    % RDM = 1 - corr(squeeze(resp_all(:, t, :)), 'rows', 'pairwise');
    % RDM = squareform(pdist(squeeze(resp_all(:, t, :))',...
    % 'correlation'));
    neural_RDMs{t} = RDM;
end

assert(sum(diff(nancol), "all")==0, "selected trials not the same across time points")
nancol = nancol(1, :);

% Plot neural RDMs

figure;
% sort trial by phoneme label 
% find non-nan trials in sorted order
[phn, ~, phnidx] = unique(Scon.phn(Scon.randidx & ~nancol));
[phntype, sortidx] = sort(phnidx);

ctr=1;
for tp = 20:30
    subplot(1, 11, ctr)
    imagesc(neural_RDMs{tp}(sortidx, sortidx));
    title(['Time point ' num2str(tp)])
    colormap(mako)

    % plot lines on RDM to separate phonemes
    hold on;
    startphn = zeros(1, length(phn));
    numphn = zeros(1, length(phn));
    for i = 1:length(phn)
        startphn(i) = find(phntype == i, 1);
        numphn(i) = sum(phntype == i);
        xline(startphn(i), 'k', 'LineWidth', 0.8)
        yline(startphn(i), 'k', 'LineWidth', 0.8)
    end

    % set the tick labels to be the phoneme labels
    set(gca, 'FontSize', 13, 'Xtick', startphn+numphn/2, 'XtickLabel', phn, ...
        'XTickLabelRotation', 90, 'Ytick', startphn+numphn/2, 'YtickLabel', phn);
    ctr = ctr + 1;
end

% -------------------------- Compute partial correlation --------------------------

% compute partial correlation between neural and acoustic RDMs

partialrho = cell(length(Scons), 1);
partialpv = cell(length(Scons), 1);
for s = 1:length(Scons)
    partialrho{s} = nan(length(neural_RDMs), length(RDMs(1, :)));
    partialpv{s} = nan(length(neural_RDMs), length(RDMs(1, :)));
    for tp = 1:length(neural_RDMs)
        % compute partial correlation between neural and acoustic RDMs

        % vectorize the upper triangle of the RDMs for partial correlation
        neuralRDM = squareform(neural_RDMs{tp});
        acousticRDMs = cellfun(@(x) squareform(x(~nancol, ~nancol)), ...
            RDMs(s, :), 'UniformOutput', false);

        % partial correlation between neural and each acoustic RDMs     
        data = cat(2, neuralRDM, acousticRDMs{:});
        [rho, pval] = partialcorr(cat(1, neuralRDM, ...
            acousticRDMs{1:end})', 'Type', 'Spearman');
        
        % fill in partial correlation
        partialrho{s}(tp, :) = rho(1, 2:end);
        partialpv{s}(tp, :) = pval(1, 2:end);
    end
end

% plot partial correlation
figure;
yticklabels = [cellfun(@(x) ['spectro-' num2str(tps(x, 1)) 'ms'], ...
    num2cell(1:size(tps, 1)), 'UniformOutput', false) {'phonetic feature', 'phoneme'}];

subplot(1, 2, 1)
imagesc(partialrho{1}');
set(gca, 'FontSize', 13, 'Xtick', 1:length(neural_RDMs), ...
    'XtickLabel', neurtps, 'Ytick', 1:length(RDMs(1, :)), 'YtickLabel', yticklabels);

subplot(1, 2, 2)
imagesc(partialpv{1}');
set(gca, 'FontSize', 13, 'Xtick', 1:length(neural_RDMs), ...
    'XtickLabel', neurtps, 'Ytick', 1:length(RDMs(1, :)), 'YtickLabel', yticklabels);
caxis([0 0.05])

save(['partialcorr/partialcorr_' Scons{1}.name '_' subjectgroup '_' phngroup '.mat'], 'partialrho', 'partialpv', ...
    'neural_RDMs', 'RDMs', 'nancol');
    
clearvars -except Dcons TDcons dimex_details timit_details neurtps tps SIDs subjectgroup ...
    datapath dpath zFolder bef aft exclude_phns RDMs neural_RDMs nancol phngroup

%% -------------------------- Plot partial correlation --------------------------

% load in partial correlation for same language, across subjects
corpus = {'dimex', 'timit'};
groups = {'spanishmono', 'englishmono'};
Scons = {Dcons, TDcons};
rho = cell(length(Scons), 1);

for s = 1:length(Scons)
    figure;
    for subj = 1:2 
        % if the file exists, load it
        filename = ['partialcorr/partialcorr_' corpus{s} '_' groups{subj} '_' phngroup '.mat'];
        if exist(filename, 'file') 
            load(filename, 'partialrho', 'partialpv', ...
                'neural_RDMs', 'RDMs', 'nancol');
        else
            continue
        end
        
        % plot partial correlation
        subplot(1, 2, subj)
        imagesc(smoothdata(partialrho{1}', 2, 'gaussian'));
        rho{subj} = partialrho{1};
        title(groups{subj});

        yticklabels = [cellfun(@(x) ['spectro-' num2str(tps(x, 1)) 'ms'], ...
            num2cell(1:size(tps, 1)), 'UniformOutput', false) {'phonetic feature', 'phoneme'}];
        set(gca, 'FontSize', 13, 'Xtick', 1:length(neural_RDMs), ...
            'XtickLabel', neurtps, 'Ytick', 1:length(RDMs(1, :)), 'YtickLabel', yticklabels);

    end
    sgtitle(Scons{s}.name);

    % plot partial correlation overlaid for english / spanish
    feats = 5:6;
    featnames = {'phonetic feature', 'phoneme'};
    colors = getColorsCrossComp(1);
    
    ctr = 1;
    figure;
    for feat = feats
        subplot(1, length(feats), ctr)
        for subj = 1:2
            plot(smoothdata(rho{subj}(:, feat)'), 'LineWidth', 2, 'Color', colors(subj, :)); hold on;
        end
        legend(groups);
        title(featnames{ctr});
        ctr = ctr + 1;
        set(gca, 'FontSize', 13);
    end
    sgtitle(Scons{s}.name);
end

clearvars -except Dcons TDcons dimex_details timit_details tps neurtps SIDs subjectgroup ...
    datapath dpath zFolder bef aft exclude_phns RDMs neural_RDMs nancol phngroup

%% -------------------------- Functions --------------------------

function [Scon, removedPhns] = removePhns(Scon, mintrials)
    % remove all phones for which there are less than 10 trials
    phn = unique(Scon.phn);
    removeidx = false(1, length(Scon.phn));
    for i = 1:length(phn)
        if sum(ismember(Scon.phn, phn{i})) < mintrials
            % add all trials to the removeidx
            removeidx = removeidx | ismember(Scon.phn, phn{i});            
        end
    end
    removedPhns = unique(Scon.phn(removeidx));

    % go through each field and remove the trials
    fields = fieldnames(Scon);
    for f = 1:length(fields)
        if size(Scon.(fields{f}), 3) == length(Scon.phn)
            Scon.(fields{f})(:, :, removeidx) = [];
        elseif size(Scon.(fields{f}), 2) == length(Scon.phn)
            Scon.(fields{f})(:, removeidx) = [];
        elseif startsWith(fields{f}, 'EC')
            Scon.(fields{f}).resp(:, :, removeidx) = [];
            if isfield(Scon.(fields{f}), 'predResp')
                Scon.(fields{f}) = rmfield(Scon.(fields{f}), 'predResp');
            end
        end
    end
end

function [Scon] = phn2feat(Scon, corpus_details)
    % use corpus details to convert phns to phonetic features
    % get phonetic features for each phoneme
    phnfeat = nan(length(corpus_details.features.names), length(Scon.phn));

    % Note: voicing label is incorrect for dimex
    for i = 1:length(Scon.phn)
        pidx = find(ismember(corpus_details.phnnames, Scon.phn{i}));
        phnfeat(:, i) = corpus_details.features.mat(:, pidx);
    end
    Scon.phnfeat = phnfeat;
end

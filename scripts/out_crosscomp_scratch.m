%% CROSSCOMP SCRATCH

%% ----------------------- FIX ----------------------------------------------
%% Plot single native brain 

load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat', 'fvals');

% initialize design electrode structure
desel=struct();
desel.conds = 1;
desel.sz = 30;
desel.cols = [0 0 0];
desel.labels = {'none'};
for s={'EC221'}
    SID = s{1};
    desel.(SID).elid = 1:length(fvals.(SID));
    desel.(SID).condition = ones(1, length(fvals.(SID)));
    plotNativeElec(s, desel);
end

clearvars -except *all subj *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *D*;

%% TO CLEAN: TRF electrode encoding

allidx_timit = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat', 'allidx');
allidx_dimex = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat', 'allidx');
allidx_asccd = load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat', 'allidx');

% load in TRF models
% with pitch change feature
modelnames={'onset_phnfeatonset_maxDtL', ...
    'phnfeatonset_maxDtL', ...
    'onset_phnfeatonset', ...
    'onset_maxDtL', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0Bin_relPitchBin_F0ChangeBin', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0Bin', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_relPitchBin_F0ChangeBin', ...
    'aud', ...
    'onset', ...
    'onset_maxDtL_vowelOnset_aud', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ...
    'onset_maxDtL_formantMedOnset', ...
    'onset_phnfeatConsOnset_maxDtL'};

% modelnames={'onset_phnfeatonset_maxDtL', ...
%     'phnfeatonset_maxDtL', ...
%     'onset_phnfeatonset', ...
%     'onset_maxDtL', ...
%     'onset_maxDtL_maxDtLOnset_vowelOnset_F0_relPitchBin', ...
%     'onset_maxDtL_maxDtLOnset_vowelOnset_F0', ...
%     'onset_maxDtL_maxDtLOnset_vowelOnset_relPitchBin'};

% determine unique variance per feature and primary encoding
varnames = {'SID', 'el', 'ls', 'base_rsq', 'spect_rsq', 'basepitch_rsq', ...
    'uv_onset', 'uv_pr', 'uv_phn', 'uv_f0', 'uv_relf0', 'uv_formant', 'uv_consonant'};
encodings =  array2table(zeros(0,length(varnames)), 'VariableNames', varnames);

[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
SIDs = [sSIDs, eSIDs, mSIDs];
SIDs(ismember(SIDs, {'EC116'}))=[]; % for now remove EC116 % , 'HS11'
for s = SIDs
    SID = s{1}; 
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs}));
    corpusStrf=loadMultModelStrf(SID, modelnames, 'timit', datapath, 1);
    
    % models without pitch
    base = (corpusStrf{1}.meanTestR.^2)';
    spect = (corpusStrf{8}.meanTestR.^2)';
    onset = (corpusStrf{9}.meanTestR.^2)';
    % from Yulia's paper
    uvOnset = onset - spect;
    uvPhnfeat = base - spect; 
    
    % models with pitch
    base_wtpitch = (corpusStrf{5}.meanTestR.^2)';    
    uvRelF0 = base_wtpitch - corpusStrf{6}.meanTestR.^2';
    uvF0 = base_wtpitch - corpusStrf{7}.meanTestR.^2';

    % formant and consonant models
    uvFormant = corpusStrf{11}.meanTestR.^2' - corpusStrf{13}.meanTestR.^2';
    uvConsonant = corpusStrf{11}.meanTestR.^2' - corpusStrf{12}.meanTestR.^2';
    uvPeakrate = corpusStrf{1}.meanTestR.^2' - corpusStrf{3}.meanTestR.^2';
    
    % find which electrodes meet threshold rsq for the basemodel
    %els = find(base > 0.01);
    %els = find(corpusStrf{10}.meanTestR.^2' > 0.025);
    if isfield(allidx_timit.allidx, SID) && isfield(allidx_dimex.allidx, SID)
        els = union(allidx_timit.allidx.(SID), allidx_dimex.allidx.(SID));
    elseif isfield(allidx_timit.allidx, SID) && isfield(allidx_asccd.allidx, SID)
        els = union(allidx_timit.allidx.(SID), allidx_asccd.allidx.(SID));
    elseif isfield(allidx_timit.allidx, SID)
        els = allidx_timit.allidx.(SID);
    else
        els = allidx_dimex.allidx.(SID);
    end
    
    sids = repmat({SID}, length(els), 1);
    lss = repmat(ls, length(els), 1);
    
    %uvPeakrate(els)
    tmp = table(sids, els, lss, base(els), spect(els), base_wtpitch(els), uvOnset(els), ...
            uvPeakrate(els),  uvPhnfeat(els), uvF0(els),  uvRelF0(els), ...
            uvFormant(els), uvConsonant(els), ...
            'VariableNames', varnames);
    encodings = [encodings; tmp];
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;

%% TO CLEAN: Discretize encodings
% at least one unique variance > 0.01 as a threshold
% pitch features calculated separately below
uvs = [encodings.uv_onset encodings.uv_pr encodings.uv_phn ...
    encodings.uv_relf0 encodings.uv_f0]';

% initial preferences, primary encoding
[~, ~] = max(uvs);
pref = zeros(height(encodings), 1);

% onsets
pref(encodings.uv_onset>0) = 1;
% peakrate and features
pref(encodings.uv_onset<0 & encodings.uv_phn>0) = 3;

% unique variance for pitch features calculated differently
thresh = 0.0065;
% if uv is greater than 0, primary feature is pitch (both or rel or abs)
pref(encodings.uv_f0 < 0 & encodings.uv_relf0 > thresh) = 4; % relative
pref(encodings.uv_f0 > thresh & encodings.uv_relf0 < 0) = 5; % absolute
pref(encodings.uv_f0 > thresh & encodings.uv_relf0 > thresh) = 6; % both

% 1 - onset, 2- peakrate, 3 - phonetic features, 4 - rel, 5  - abs, 6 - both 
% pref(~any(uvs > 0.02)) = 0;

encodings.primary = pref;

% visualization: plot distributions
figure;
colors = [0.45 0.67 0.18; ... % onset
%           0.7 0.3 0.7; ... % peakrate
          0.8 0.5 0.7; ... % phonetic feature
          0 0.44 0.74; ... % relative pitch
          0.8 0.2 0.2; ... % absolute pitch
          0.92 0.69 0.12];    % both
idx = encodings.primary~=0;
% ls 1 = spanish, ls 2 = english
count = [histcounts(encodings.primary(encodings.ls==1 & idx), 6); ...
    histcounts(encodings.primary(encodings.ls==2 & idx), 6); ...
    histcounts(encodings.primary(encodings.ls==3 & idx), 6)];
% remove peakrate as its own category
count = [count(:, 1) count(:, 2)+count(:, 3) count(:, 4:6)];
percent = (count./sum(count, 2))*100;
bar(percent, 'hist', 'EdgeColor', 'none');

% plotting details
xticklabels({'Spanish', 'English', 'Mandarin'});
yticks(0:25:50);
ylabel('Electrode Count (%)');
colormap(colors);
set(gca, 'FontSize', 13);
labels = {'onset', 'phonetic features + peak rate', ...
    'relative pitch', 'absolute pitch', 'both'};
legend(labels);
box off;

% alternate visualization
figure;
xl = [-0.02 0.1];
yl = [0 0.21];
idx = encodings.base_rsq>0.15;
colors = brewermap(3, 'Dark2');
subplot(1, 3, 1);
scatter(encodings.uv_onset(idx), max(encodings.uv_pr(idx), ...
    encodings.uv_phn(idx)), 45, encodings.ls(idx), 'filled', ...
    'MarkerFaceAlpha', 0.4);
ylabel('peak rate + phonetic'); xlabel('onset');
ylim(yl);
xlim(xl);
h = refline(1, 0); h.Color = 'k';
vertline(0); horzline(0);
yticks(0:0.1:0.2);
set(gca, 'FontSize', 15);

subplot(1, 3, 2);
scatter(max(encodings.uv_f0(idx), encodings.uv_relf0(idx)), ...
    encodings.uv_onset(idx), 45, encodings.ls(idx), 'filled', ...
    'MarkerFaceAlpha', 0.4);
ylabel('onset'); xlabel(' absolute + rel pitch');
ylim(yl);
xlim(xl);
h = refline(1, 0); h.Color = 'k';
vertline(0); horzline(0);
yticks(0:0.1:0.2);
set(gca, 'FontSize', 15);

subplot(1, 3, 3);
scatter(max(encodings.uv_f0(idx), encodings.uv_relf0(idx)), ...
    max(encodings.uv_pr(idx), encodings.uv_phn(idx)), 45, ...
    encodings.ls(idx), 'filled','MarkerFaceAlpha', 0.4);
ylabel('peak rate + phonetic'); xlabel(' absolute + rel pitch');
ylim(yl);
xlim(xl);
h = refline(1, 0); h.Color = 'k';
vertline(0); horzline(0);
colormap(colors);
yticks(0:0.1:0.2);
set(gca, 'FontSize', 15);

%legends('spanish', 'english', 'mandarin');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;


%% Discretize encodings v2
% at least one unique variance > 0.01 as a threshold
% pitch features calculated separately below
uvs = [encodings.uv_onset encodings.uv_pr encodings.uv_phn ...
    encodings.uv_relf0 encodings.uv_f0]';

% initial preferences, primary encoding
[~, ~] = max(uvs);
pref = zeros(height(encodings), 1);
thresh = 0.005;

% onsets
pref(encodings.uv_onset>0) = 1;
% peakrate and features
%pref(encodings.uv_onset<0 & encodings.uv_phn>0) = 3;
pref(encodings.uv_onset<0 & encodings.uv_pr>0)=2;

% unique variance for pitch features calculated differently
% if uv is greater than 0, primary feature is pitch (both or rel or abs)
% pref(encodings.uv_f0 < 0 & encodings.uv_relf0 > thresh) = 4; % relative
% pref(encodings.uv_f0 > thresh & encodings.uv_relf0 < 0) = 5; % absolute
% pref(encodings.uv_f0 > thresh & encodings.uv_relf0 > thresh) = 6; % both
pref(encodings.uv_f0 > thresh | encodings.uv_relf0 > thresh)=3;

% if consonant feature is greater than the threshold
pref(encodings.uv_formant > thresh & ...
    (encodings.uv_formant>encodings.uv_f0 | ...
    encodings.uv_formant>encodings.uv_relf0)) = 5;
pref(encodings.uv_consonant>thresh) = 4; % consonant


% 1 - onset, 2- peakrate, 3 - pitch, 4 - consonant, 5  - formant
% pref(~any(uvs > 0.02)) = 0;
encodings.primary = pref;

% visualization: plot distributions
figure;
colors = [0.45 0.67 0.18; ... % onset
%           0.7 0.3 0.7; ... % peakrate
          0.8 0.5 0.7; ... % phonetic feature
          0 0.44 0.74; ... % relative pitch
          0.8 0.35 0.15; ... % absolute pitch
          0.92 0.69 0.12];    % both
idx =  encodings.primary~=0;
% ls 1 = spanish, ls 2 = english
count = [histcounts(encodings.primary(encodings.ls==1 & idx), 0.5:5.5); ...
    histcounts(encodings.primary(encodings.ls==2 & idx), 0.5:5.5); ...
    histcounts(encodings.primary(encodings.ls==3 & idx), 0.5:5.5)];
percent = (count./sum(count, 2))*100;
bar(percent, 'hist', 'EdgeColor', 'none');

% plotting details
xticklabels({'Spanish', 'English', 'Mandarin'});
yticks(0:25:50);
ylabel('Electrode Count (%)');
colormap(colors);
set(gca, 'FontSize', 13);
labels = {'onset', 'peak rate', ...
    'relative/absolute pitch', 'consonant', 'formant'};
legend(labels);
box off;
% 1 - onset, 2- peakrate, 3 - pitch, 4 - consonant, 5  - formant

% alternate visualization
figure;
xl = [-0.02 0.1];
yl = [0 0.21];
idx = encodings.base_rsq>0.15;
colors = brewermap(3, 'Dark2');
subplot(1, 3, 1);
scatter(encodings.uv_onset(idx), max(encodings.uv_pr(idx), ...
    encodings.uv_phn(idx)), 45, encodings.ls(idx), 'filled', ...
    'MarkerFaceAlpha', 0.4);
ylabel('peak rate + phonetic'); xlabel('onset');
ylim(yl);
xlim(xl);
h = refline(1, 0); h.Color = 'k';
vertline(0); horzline(0);
yticks(0:0.1:0.2);
set(gca, 'FontSize', 15);

subplot(1, 3, 2);
scatter(max(encodings.uv_f0(idx), encodings.uv_relf0(idx)), ...
    encodings.uv_onset(idx), 45, encodings.ls(idx), 'filled', ...
    'MarkerFaceAlpha', 0.4);
ylabel('onset'); xlabel(' absolute + rel pitch');
ylim(yl);
xlim(xl);
h = refline(1, 0); h.Color = 'k';
vertline(0); horzline(0);
yticks(0:0.1:0.2);
set(gca, 'FontSize', 15);

subplot(1, 3, 3);
scatter(max(encodings.uv_f0(idx), encodings.uv_relf0(idx)), ...
    max(encodings.uv_pr(idx), encodings.uv_phn(idx)), 45, ...
    encodings.ls(idx), 'filled','MarkerFaceAlpha', 0.4);
ylabel('peak rate + phonetic'); xlabel(' absolute + rel pitch');
ylim(yl);
xlim(xl);
h = refline(1, 0); h.Color = 'k';
vertline(0); horzline(0);
colormap(colors);
yticks(0:0.1:0.2);
set(gca, 'FontSize', 15);

%legends('spanish', 'english', 'mandarin');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* encodings labels;

%% Plot discretized encoding by subject

colors = [0.45 0.67 0.18; ... % onset
%           0.7 0.3 0.7; ... % peakrate
          0.8 0.5 0.7; ... % phonetic feature
          0 0.44 0.74; ... % relative pitch
          0.8 0.35 0.15; ... % absolute pitch
          0.92 0.69 0.12];    % both
figure;
for ls = [1, 2, 3]
    percent = nan(5, length(unique(encodings.SID(encodings.ls==ls))));
    ctr = 1;
    for s = unique(encodings.SID(encodings.ls==ls))'
        SID = s{1};
        idx = strcmp(encodings.SID, SID);
        
        % get proportion of electrodes for each primary encoding
        count = histcounts(encodings.primary(idx), 0.5:5.5);
        percent(:, ctr) = (count./sum(count))*100;

        ctr = ctr + 1;
    end

    % plot bar graph with error bars for this group of subjects
    idx = [1:2 4:5];
    x = (1:length(idx)) + (ls-2)*7;
    h = bar(x, mean(percent(idx, :), 2)); hold on;

    % have each bar be a different color
    for i = 1:length(idx)
        h.FaceColor = 'flat';
        h.CData(i, :) = colors(idx(i), :);

        % remove line color
        h.EdgeColor = 'none';
    end

    err = nansem(percent(idx, :), 2);
    er = errorbar(x,mean(percent(idx, :), 2),err);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.CapSize = 0;
    er.HandleVisibility = 'off';

 end

 labels = {'onset', 'peak rate', ...
    'relative/absolute pitch', 'consonant', 'formant'};
legend(labels);

%% TO CLEAN: plot on MNI brain
% initialize design electrode structure
desel=struct();
desel.conds = 0:5; % or 0:6
desel.sz = [5 repmat(60, 1, length(desel.conds)-1)];
colors = [0.45 0.67 0.18; ... % onset
%           0.7 0.3 0.7; ... % peakrate
          0.8 0.5 0.7; ... % phonetic feature
          0 0.44 0.74; ... % relative pitch
          0.8 0.2 0.2; ... % absolute pitch
          0.92 0.69 0.12];    % both
desel.cols = [0 0 0; colors];
% split up peak rate and phonetic features again for MNI plotting
desel.labels = [{'none'} labels(1) {'peak rate', 'phonetic features'} ...
    labels(3:end)];
for s=unique(encodings.SID)'
    SID = s{1};
    idx = strcmp(encodings.SID, SID);
    desel.(SID).elid = encodings.el(idx);
    desel.(SID).condition = encodings.primary(idx);
    disp(unique(encodings.primary(idx)))
end

% plot histogram of count across y-axis (anterior to posterior)
for ls = 1:3
    
    lsid = find(encodings.ls == ls);
    mni_rh = [];
    if ls<7
        [mni_lh] = plotMNIElec(unique(encodings.SID(lsid)), desel, 'lh', 0);
        [mni_rh] = plotMNIElec(unique(encodings.SID(lsid)), desel, 'rh', 0);
%     else
%          % no RH HS subjects
%         [mni_lh] = plotNativeElec({'HS11'}, desel);
    end
    
    figure;
    
    label = {'onsets', 'phnfeat', 'pitch (F0 + rel F0)'};
    
    if ls < 3
        % normalize before splitting by conditions
        y = normalize(mni_lh.y);       
    else % no warped elecs for mandarin subjects
        % normalize by subject
        y = cell2mat(cellfun(@(x) normalize(mni_lh.y(strcmp(mni_lh.SID, x))), ...
            unique(mni_lh.SID), 'UniformOutput', false));
    end
    stg = strcmp(mni_lh.anatomy, 'superiortemporal');
    mni_conds = mni_lh.cond;
    
    if ~isempty(mni_rh)
        % reflect y in rh so matches A-P direction    
        y_rh = normalize(mni_rh.y);
        %y_rh = -1*mni_rh.y;
        stg_rh = strcmp(mni_rh.anatomy, 'superiortemporal');
        y = [y; y_rh];
        stg = [stg; stg_rh];
        mni_conds = [mni_lh.cond; mni_rh.cond];
    end

    ctr = 1;
    for cond = 1:5   
        [f,xi] = ksdensity(y(mni_conds == cond & stg), ...
            'Function','pdf','Bandwidth',0.3); 

%         plot(xi, f, 'k', 'LineWidth', 1.5);  hold on;
%         box off;
        if cond == 4 % include all pitch
            idx =  ismember(mni_conds, [4 5 6]) & stg;
        else
            idx =  mni_conds == cond & stg;
        end
        
        scatter(y(idx), ones(sum(idx), 1)*ctr - 0.1 + ...
            rand(sum(idx), 1)*0.1, 45, colors(cond, :), ...
            'filled', 'MarkerFaceAlpha', 0.6); hold on;
        line([median(y(idx)) median(y(idx))], [ctr - 0.1 ctr + 0.1], ...
            'Color', colors(cond, :), 'LineWidth', 2)
        ctr = ctr+0.5;    
    end
    ylim([0.8 3.2])
    xlim([prctile(y, 1.5) prctile(y, 98.5)]);
    
    yticks(1:0.5:3);
    yticklabels({'onsets', 'peakrate', 'pitch', 'consonant', 'vowel'})
    xticks([prctile(y, 2.5) prctile(y, 97.5)]);
    xticklabels({'posterior', 'anterior'});
    set(gca, 'FontSize', 13, 'TickLabelInterpreter', 'latex');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* encodings labels;


%% ------------- Neural Encoding: ERPs ------------------------------------
%%  stressed formants ploitting

lims = [159 1100; 500 3000];
stress = {'Unstressed', 'Stressed'};
ctr = 1;
figure;
for formant = 1:2    
    for s = 0:1  
        subplot(2, 2, formant + floor(formant/2)+s)
        histogram(Dvow.formantVals(formant, Dvow.stress==s ...
            & Dvow.meanf0>170), 'EdgeColor', 'none'); hold on;

        histogram(Dvow.formantVals(formant, Dvow.stress==s & ...
            Dvow.meanf0<=170), 'EdgeColor', 'none');
        title(stress{s+1});
        box off;
        set(gca, 'FontSize', 15);
        xlabel(['F' num2str(formant) ' (Hz)']);
        xlim(lims(formant, :));
        
        yticks(0:200:800)
        ylim([0 560])
    end
end
%% Look at /i/ vs. /I/ distinction

[sSIDs, eSIDs, ~, mSIDs] = getSIDinfo();
% [sSIDs, eSIDs, {'HS11', 'HS9', 'HS10'}];

% find all vowel discriminating electrodes
allidx = getElecs(TDvow, [sSIDs eSIDs mSIDs], [], 'timit', 'ftest', struct());

% merge both sets of electrodes into a table
%vpair = {'iy', 'ih'};
vpair = {'ey', 'ih'};
vpair_idx = ismember(TDvow.vowel, vpair)';
varnames = {'SID', 'el', 'ls', 'fval', ['hg_' vpair{1}], ['hg_' vpair{2}], 'maxfstat_tp'};
pair_encoding =  array2table(zeros(0,7), 'VariableNames', varnames);

ctr = 1;
for s = [sSIDs eSIDs mSIDs]
    disp(['loading subject ' num2str(ctr) '/' ...
        num2str(length([sSIDs eSIDs mSIDs]))]);
    SID = s{1}; 
    % spanish - 1, english - 2, mandarin - 3
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs}));
    
    try 
        els = allidx.(SID);

        % extract response from data structure
        resp = TDvow.(SID).resp;

        % remove all nans and all other vowels 
        nanidx = squeeze(any(isnan(resp(els, :, :)), [1 2])) | ...
            isnan(TDvow.vowelType)';

       [Fstat, ~, ~, df1, df2] = Fstat_TIMIT(resp(els, :, vpair_idx & ~nanidx), ...
            TDvow.vowelType(vpair_idx & ~nanidx),  ...
            TDvow.vowelType(vpair_idx & ~nanidx));

        % get significant F threshold
        fthresh = finv(1-0.0001, df1, df2);

        [fvals, tp] = max(Fstat, [], 2);
        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        
        hg = nan(length(els), 2);
        hg(:, 1) = arrayfun(@(x) mean(squeeze(resp(els(x), tp(x), ...
            ismember(TDvow.vowel, vpair{1}))), 'omitnan'),1:length(els)); 
        hg(:, 2) = arrayfun(@(x) mean(squeeze(resp(els(x), tp(x), ...
            ismember(TDvow.vowel, vpair{2}))), 'omitnan'),1:length(els)); 

        tmp = table(sids, els, lss, fvals, hg(:, 1), hg(:, 2), (tp-bef)/100, ...
                'VariableNames', varnames);
        pair_encoding = [pair_encoding; tmp];
    catch
        warning(['No electrode information for subject: ' SID]);
    end
    ctr = ctr + 1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh;

%% comparing fstats across 
pair = {'ih', 'iy'};

% test whether there is /i/ vs. /I/ distinction
figure;
subplot(2, 2, 1);
colors = getColorsCrossComp(1);
for ls = 1:3
    idx = pair_encoding.ls==ls;
    scatter(ls-0.1+ones(1, sum(idx)).*rand(1, sum(idx))*0.2, ...
        pair_encoding.fval(idx), 45, colors(ls, :), ...
        'filled', 'MarkerFaceAlpha', 0.7); hold on;
    line([ls-0.15 ls+0.15], [mean(pair_encoding.fval(idx)) ...
        mean(pair_encoding.fval(idx))], 'LineWidth', 1.5, 'Color', 'k')
end
set(gca, 'yscale', 'log', 'FontSize', 13);
ylabel('/iy/ vs. /ih/ Fstat');
yline(fthresh, '-k');
xticks(1:3);
xlim([0.5 3.5])
xticklabels({'Spanish', 'English', 'Mandarin'});

subplot(2, 2, 2);
idx = pair_encoding.fval>3;
scatter(pair_encoding.hg_iy(idx), pair_encoding.hg_ih(idx), 45, ...
    pair_encoding.ls(idx), 'filled', 'MarkerFaceAlpha', 0.6);
colormap(colors); 
set(gca, 'FontSize', 13);
ylabel('/ih/ HGA peak (z)');
xlabel('/iy/ HGA peak (z)');
xticks([0 2]);
yticks([0 2]);
h=refline(1, 0);
h.Color = 'k';
h.LineWidth = 1.5;
xlim([0 2.2]);
ylim([0 2.2]);

subplot(2, 2, 3);
for ls = 1:3
    idx = pair_encoding.ls==ls;
%     h = histogram(pair_encoding.maxfstat_tp(idx), 15); hold on;
%     h.FaceAlpha = 0.5;
%     h.FaceColor = colors(ls, :);
    scatter(pair_encoding.maxfstat_tp(idx), pair_encoding.fval(idx), 45, ...
        colors(ls, :), 'filled', 'MarkerFaceAlpha', 0.6); hold on;
end
ylabel('Fstat');
xlabel('Time of maximal Fstat (s)');
xline(0, '-k');
yline(fthresh, '-k');
ylim([0 30]);
set(gca, 'FontSize', 13);

figure;
for ls = 1:3
    idx = pair_encoding.ls==ls;
    h = histogram(pair_encoding.(['hg_' pair{1}])(idx) - ...
        pair_encoding.(['hg_' pair{2}])(idx), 15); hold on;
    h.FaceAlpha = 0.5;
    h.FaceColor = colors(ls, :);
end
ylabel('Electrode count');
xlabel(['/' pair{1} '/ HGA - /' pair{2} '/ HGA']);
set(gca, 'FontSize', 13);

% look at ERPs
%plotVowelErp(TDvow, 'EC172',[58 74], 'timit', [], [1 0 0; 0 0 0], {'ih', 'iy'});

%% Look at /v/ vs. /f/ and /z/ vs /s/ distinction

% make Dcons (consonant structures)
%Dcons = loadDDcons('timit', 20, 50, SIDs);
[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
SIDs = [sSIDs, eSIDs, mSIDs];

% find all electrodes with phonetic feature encoding (uv rsq>0.01)
modelnames = {'onset_phnfeatonset_maxDtL', 'onset_maxDtL'};
for s = SIDs
    SID = s{1};
    corpusStrf=loadMultModelStrf(SID, modelnames, 'timit', datapath, 1);
    
    % models without pitch
    uv_phnfeat = corpusStrf{1}.meanTestR.^2 - corpusStrf{2}.meanTestR.^2; 
    allidx.(SID) = find(uv_phnfeat>0.025)';
end

% merge both sets of electrodes into a table
% test /v/ vs /f/ distinction
% test /z/ vs /s/ distinction
pair = {'z', 's'};
pair_idx = ismember(Dcons.phn, pair)';
varnames = {'SID', 'el', 'ls', 'fval', ['hg_' pair{1}], ['hg_' pair{2}], ...
    'fstat_maxtp'};
pair_encoding =  array2table(zeros(0,7), 'VariableNames', varnames);

ctr = 1;
for s = [sSIDs eSIDs mSIDs]
    disp(['loading subject ' num2str(ctr) '/' ...
        num2str(length([sSIDs eSIDs mSIDs]))]);
    SID = s{1}; 
    % spanish - 1, english - 2, mandarin - 3
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs}));
    
    try 
        els = allidx.(SID);

        % extract response from data structure
        resp = Dcons.(SID).resp;

        % remove all nans and all other vowels 
        nanidx = squeeze(any(isnan(resp(els, :, :)), [1 2])) | ...
            isnan(Dcons.phnType)';

       [Fstat, ~, ~, df1, df2] = Fstat_TIMIT(resp(els, :, pair_idx & ~nanidx), ...
            Dcons.phnType(pair_idx & ~nanidx),  ...
            Dcons.phnType(pair_idx & ~nanidx));

        % get significant F threshold
        fthresh = finv(1-0.0001, df1, df2);

        [fvals, tp] = max(Fstat, [], 2);
        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        
        hg = nan(length(els), 2);
        hg(:, 1) = arrayfun(@(x) mean(squeeze(resp(els(x), tp(x), ...
            ismember(Dcons.phn, pair{1}))), 'omitnan'),1:length(els)); 
        hg(:, 2) = arrayfun(@(x) mean(squeeze(resp(els(x), tp(x), ...
            ismember(Dcons.phn, pair{2}))), 'omitnan'),1:length(els)); 

        tmp = table(sids, els, lss, fvals, hg(:, 1), hg(:, 2), (tp-bef)/100, ...
                'VariableNames', varnames);
        pair_encoding = [pair_encoding; tmp];
        clear tmp hg fvals tp sids lss Fstat nanidx resp
    catch
        warning(['No electrode information for subject: ' SID]);
    end
    ctr = ctr + 1;
end

%% visualize whether there is /i/ vs. /I/ distinction

%pair = {'z', 's'};
%pair = {'v', 'f'};
pair = {'iy', 'ih'};
pair_idx = ismember(Dcons.phn, pair)';
figure;
subplot(2, 2, 1);
colors = getColorsCrossComp(1);
for ls = 1:3
    idx = pair_encoding.ls==ls;
    scatter(ls-0.1+ones(1, sum(idx)).*rand(1, sum(idx))*0.2, ...
        pair_encoding.fval(idx), 45, colors(ls, :), ...
        'filled', 'MarkerFaceAlpha', 0.6); hold on;
    line([ls-0.15 ls+0.15], [mean(pair_encoding.fval(idx)) ...
        mean(pair_encoding.fval(idx))], 'LineWidth', 1.5, 'Color', 'k')
end
set(gca, 'yscale', 'log', 'FontSize', 13);
ylabel(['/' pair{1} '/ vs. /' pair{2} '/ Fstat']);
yline(fthresh, '-k');
xticks(1:3);
xlim([0.5 3.5])
xticklabels({'Spanish', 'English', 'Mandarin'});

subplot(2, 2, 2);
scatter(pair_encoding.(['hg_' pair{1}]), pair_encoding.(['hg_' pair{2}]), 45, ...
    pair_encoding.ls, 'filled', 'MarkerFaceAlpha', 0.6);
colormap(colors); 
set(gca, 'FontSize', 13);
xlabel(['/' pair{1} '/ HGA peak (z)']);
ylabel(['/' pair{2} '/ HGA peak (z)']);
xticks([0 2]);
yticks([0 2]);
xlim([0 2]);
ylim([0 2]);
h=refline(1, 0);
h.Color = 'k';
h.LineWidth = 1.5;

subplot(2, 2, 3);
for ls = 1:3
    idx = pair_encoding.ls==ls;
    h = histogram(pair_encoding.fstat_maxtp(idx), 15); hold on;
    h.FaceAlpha = 0.5;
    h.FaceColor = colors(ls, :);
    h.EdgeColor = 'none';
%     scatter(pair_encoding.fstat_maxtp(idx), pair_encoding.fval(idx), 45, ...
%         colors(ls, :), 'filled', 'MarkerFaceAlpha', 0.3); hold on;
end
ylabel('Count');
xlabel('Time of maximal Fstat (s)');
xline(0, '-k');
yticks(0:20:50);
yline(fthresh, '-k');
set(gca, 'FontSize', 13);
box off;

% [h, p] = ttest2(pair_encoding.fval(pair_encoding.ls==1), ...
%     pair_encoding.fval(pair_encoding.ls==2));

% subplot(1, 3, 3);
% figure;
% for ls = 1:3
%     idx = pair_encoding.ls==ls;
%     h = histogram(pair_encoding.(['hg_' pair{1}])(idx) - ...
%         pair_encoding.(['hg_' pair{2}])(idx), 15); hold on;
%     h.FaceAlpha = 0.5;
%     h.FaceColor = colors(ls, :);
% end
% ylabel('Electrode count');
% xlabel(['/' pair{1} '/ HGA - /' pair{2} '/ HGA']);
% set(gca, 'FontSize', 13);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons;
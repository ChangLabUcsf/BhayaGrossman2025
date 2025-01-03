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
SIDs = [sSIDs, eSIDs];

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

% loading in word-level subject data
TDwrd = loadDwrdMulti('timit', bef, aft, {}, timit_details);
Dwrd = loadDwrdMulti('dimex',  bef, aft, {}, dimex_details);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% ----------------- EXPLORATION OF WORD UNIQUE VARIANCE ------------------
%% Load in unique variance  of word features over spectrogram

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

%% Load in unique variance of phoneme surprisal + word boundary over feature model

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

modelnames_dimex={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_spSurpNoOnsBin', ... %remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_spSurpNoOnsBin', ... % remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin'}; % full

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

SIDs = [sSIDs eSIDs bSIDs mSIDs];
[wordsurp_encoding] = loadUniqueVarTbl(modelnames_timit, modelnames_dimex, SIDs);
wordsurp_details.featureOrd ...
    = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
imgall = load_allimgdata;
wordsurp_encoding.hemi = cellfun(@(x) imgall.(x).hemi, ...
    wordsurp_encoding.SID, 'UniformOutput', false);
wordsurp_details.models_dimex = modelnames_dimex;
wordsurp_details.models_timit = modelnames_timit;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% Load in and print anatomy percentages for surprisal + word electrodes

[img] = load_allimgdata;
anatomy = cell(height(wordsurp_encoding), 1);
hemisphere = cell(height(wordsurp_encoding), 1);
for i = 1:height(wordsurp_encoding)

    el = wordsurp_encoding.el(i);
    SID = wordsurp_encoding.SID(i);
    anatomy{i} = img.(SID{1}).imgNative.anatomy(el, 4);
    hemisphere{i} = img.(SID{1}).hemi;
end
wordsurp_encoding.anatomy = anatomy;
wordsurp_encoding.hemi = hemisphere;

feature = 'word+surp';
featureIdx = find(strcmp(wordsurp_details.featureOrd, feature));

for f = {'eng_uv_all', 'sp_uv_all'}
    field = f{1};

    if strcmp(field, 'eng_uv_all')
        disp('------------------ English anatomy ------------------')
    else
        disp('------------------ Spanish anatomy ------------------')
    end

    disp('Left hemisphere');
    tabulate(sort([wordsurp_encoding.anatomy{wordsurp_encoding.(field)(:, 5)>0.001 ...
        & strcmp(wordsurp_encoding.hemi,'lh')}]))

    disp('Right hemisphere');
    tabulate(sort([wordsurp_encoding.anatomy{wordsurp_encoding.(field)(:, 5)>0.001 ...
        & strcmp(wordsurp_encoding.hemi,'rh')}]))
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;


%% Visualize feature unique variance across TIMIT and DIMEx

colors = getColorsCrossComp(1);

% Specify the labels for analysis
labels = {'consonant', 'word+surp'};

% Create a figure for the scatter plots
figure('Position', [200, 100, 900, 900]);

% Set a UV threshold for filtering
uv_thresh = 0.001;

% Initialize subplot counter and axis array
ctr = 1;
ax = nan(3, 1);

% Iterate over the labels
for label = labels
    
    % Find the index of the current label in the featureOrd array
    index = find(ismember(wordsurp_details.featureOrd, label));
    
    % Iterate over the ls values
    for ls = [1 2 4]
        % Create a subplot for the current label and ls value
        ax(ctr) = subplot(length(labels), 3, ctr);
        
        % Get the corresponding unique variance values for English and Spanish
        x = wordsurp_encoding.eng_uv_all(wordsurp_encoding.ls==ls, index);
        y = wordsurp_encoding.sp_uv_all(wordsurp_encoding.ls==ls, index);
        
        % Get the subject IDs for the current ls value
        sid = cellfun(@(x) str2double(x(3:end)), ...
            wordsurp_encoding.SID(wordsurp_encoding.ls==ls));
        
        % Remove data points that do not meet the UV threshold or contain NaN values
        neg = all([x,y]<uv_thresh, 2) | isnan(x) | isnan(y);
        x(neg) = [];
        y(neg) = [];
        sid(neg) = [];

        % Perform permutation testing to compute correlation coefficient and p-value
        maxlim = prctile([x; y], 100);
        minlim = prctile([x; y], 3);  
        [r, p] = corr(x, y, 'Rows', 'complete');

        % Plot the scatter plot
        scatter3(x, y, sid, 15, colors(ls, :), 'filled');
        view(2);
        
        % Set the x and y axis limits and add labels
        xlim([minlim maxlim]);
        ylim([minlim maxlim]);
        xlabel(['English ' label]);
        ylabel(['Spanish ' label]);
        title({['r= ' num2str(r) ','], ['p=' num2str(p, 4)]});
        
        % Increment the subplot counter
        ctr = ctr + 1;
    end
end

% Link axes for phonetic features and word onset/surp
linkaxes(ax(1:3));
linkaxes(ax(4:6));

% Add reference lines and lines at 0
for i = 1:length(labels)*3
    subplot(length(labels), 3, i);
    xline(0, 'Color', 'k', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
end

% Initialize subplot counter and iterate over the labels again
ctr = 1;
for label = labels
    % Find the index of the current label in the featureOrd array
    index = find(ismember(featureOrd, label));
    
    % Iterate over the ls values
    for ls = [1 2 4]
        % Get the corresponding data for English and Spanish
        x = eng_uv_all(wordsurp_encoding.ls==ls, index);
        y = sp_uv_all(wordsurp_encoding.ls==ls, index);

        % Set a UV threshold for filtering
        uv_thresh = 0.001;
        
        % Remove data points that do not meet the UV threshold or contain NaN values
        neg = all([x,y]<uv_thresh, 2) | isnan(x) | isnan(y);
        x(neg) = [];
        y(neg) = [];

        % Get the position of the current subplot
        pos = ax(ctr).Position;
        
        % Create a new axes for the histogram
        h = axes('Parent', gcf, 'Position',[pos(1)+0.1, pos(2)+0.1, 0.1, 0.1]);
        
        % Plot the histogram of the difference between English and Spanish values
        histogram(h, max(0, x)-max(0, y), 10, 'EdgeColor', 'none', 'FaceColor', [0.7 0.7 0.7]);
        hold on;
        
        % Set the x-axis limits based on the index value
        if index < 5
            xlim([-0.04 0.04])
        else
            xlim([-0.02 0.02]);
        end
        xline(0);
        
        % Increment the subplot counter
        ctr = ctr + 1;
    end
end

% TO FIX: Draw a contour cloud over scatter points
% issue: scatter is sparse and it's sparse as you get higher UVs
%         bins = 35;
%         XEDGES = linspace(prctile(x, 1), prctile(x, 100), bins); % -0.03:0.001:0.1
%         YEDGES = linspace(prctile(y, 1), prctile(y, 100), bins);
% 
%         XCENTERS = XEDGES(1:end-1)+mean(diff(XEDGES))/2;
%         YCENTERS = YEDGES(1:end-1)+mean(diff(YEDGES))/2;
% 
%         rng(1)
%         new_x = randsample(minlim:0.0001:maxlim, 440);
%         new_y = interp1(x, y, new_x);
%         DATA1 = histcounts2([x; new_x'], [y; new_y'], XEDGES, YEDGES); % 
% 
%         DATA2 = smooth2a(DATA1, 3, 3);
%         %DATA2 = imgaussfilt(interp2(DATA1, 'cubic'),0.3); 

% XCENTINT = sort([XCENTERS, XCENTERS(1:end-1)+mean(diff(XEDGES))/2]);
%         YCENTINT = sort([YCENTERS, YCENTERS(1:end-1)+mean(diff(YEDGES))/2]);
%         h2=contourf(XCENTERS,YCENTERS, DATA2, ...
%             bins, 'edgecolor','none'); hold on; % [-5 linspace(0, 30, 34)]
%         %imagesc(XCENTERS,YCENTERS, DATA2'); hold on;
% 
%         tmp = [linspace(1, colors(ls, 1), 50); ...
%             linspace(1, colors(ls, 2), 50); linspace(1, colors(ls, 3), 50)]';
%         colormap(ax(ctr), [tmp((50-bins):end, :)]);
%         set(gca, 'YDir', 'normal', 'XDir', 'normal')

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;


%% Visualize feature unique variance across native vs. nonnative as scatter
% Specify the labels for analysis
labels = {'word+surp'};

% Create a figure for the scatter plots
figure('Position', [200, 100, 900, 900]);

% Set a UV threshold for filtering
uv_thresh = 0.001;

% Initialize subplot counter and axis array
ctr = 1;
ax = nan(3, 1);

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
        
        % Get the corresponding unique variance values for English and Spanish
        if ls == 1
            x = wordsurp_encoding.eng_uv_all(wordsurp_encoding.ls==ls, index);
            y = wordsurp_encoding.sp_uv_all(wordsurp_encoding.ls==ls, index);
        else
            x = wordsurp_encoding.sp_uv_all(wordsurp_encoding.ls==ls, index);
            y = wordsurp_encoding.eng_uv_all(wordsurp_encoding.ls==ls, index);
        end
        
        % Get the subject IDs for the current ls value
        sid = cellfun(@(x) str2double(x(3:end)), ...
            wordsurp_encoding.SID(wordsurp_encoding.ls==ls));
        
        % Remove data points that do not meet the UV threshold or contain NaN values
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
        
        scatter(x, y, 25, colors, 'filled', ...
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
        caxis([-0.05 0.05]);
        view(2);
        
        % Set the x and y axis limits and add labels
        xlim([minlim maxlim]);
        ylim([minlim maxlim]);
        xlabel(['Unfamiliar ' label]);
        ylabel(['Native ' label]);
    end

    [r, p] = corr(x_all, y_all, 'Rows', 'complete');
    title({['r= ' num2str(r) ','], ['p=' num2str(p, 4)]});
     % Increment the subplot counter
    ctr = ctr + 1;
end

% Link axes for phonetic features and word onset/surp
linkaxes(ax(1:length(labels)));h = refline(1, 0);

% Add reference lines and lines at 0
for i = 1:length(labels)
    subplot(length(labels), 1, i);
    xline(0, 'Color', 'k', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
    xlim([-0.01 0.1]);
    ylim([-0.01 0.1])
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
end

%% Visualize feature unique variance across Spanish and English speaking subjects

% Set a UV threshold for filtering
uv_thresh = 0.001;

colors = getColorsCrossComp(1);

% full model R^2
% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 

% Set a UV threshold for filtering
label = 'surp';
% label = 'formant';
% label = 'formant';
% label = 'consonant';
% label = 'word boundary';
% label = 'phonetic feature';

showScatter = false;

fields = {'eng_uv_all', 'sp_uv_all'}; 
index = find(ismember(wordsurp_details.featureOrd, label));

figure;
ctr = 1;

% Scatter plots with reference line
for ls = [1 2 4]
    ax(ctr) = subplot(1, 3, ctr);

    % Get data for English and Spanish
    x = wordsurp_encoding.(fields{1})(wordsurp_encoding.ls==ls, index);
    y = wordsurp_encoding.(fields{2})(wordsurp_encoding.ls==ls, index);

    scatter(x, y, 15, colors(ls, :), 'filled');

    maxlim = prctile([x; y], 99);
    xlim([-0.01 maxlim+0.005]);
    ylim([-0.01 maxlim+0.005]);
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    xline(0, 'Color', 'k');
    yline(0, 'Color', 'k');

    xlabel(['English ' label]);
    ylabel(['Spanish ' label]);
    title(['r: ' num2str(nancorr(x, y)) ]);
    ctr = ctr + 1;
end

titles = {'English', 'Spanish'};

% Box charts and pie charts
for f = 1:2
    field = fields{f};
    ctr = 1;
    figure;
    subplot(3, 2, [1 3 5])
    
    % Create box charts
    for ls = [1 2 4]
        lsidx = find(wordsurp_encoding.ls==ls ...
            & wordsurp_encoding.(field)(:, index)>uv_thresh);
        y = wordsurp_encoding.(field)(:, index);

        % Scatter the individual point
        if showScatter
            nsidx = find(wordsurp_encoding.ls==ls & ...
                wordsurp_encoding.(field)(:, index)<uv_thresh);
            scatter(ctr-0.2 + rand(length(nsidx), 1)*0.05, y(nsidx), ...
                5, [0.5 0.5 0.5], 'filled', 'LineWidth', 0.1, 'MarkerFaceAlpha', 0.6); hold on;
            scatter(ctr-0.2 + rand(length(lsidx), 1)*0.05, y(lsidx), 5, colors(ls, :), ...
                'filled', 'LineWidth', 0.1, 'MarkerFaceAlpha', 0.6); hold on;
        end

        % Create the boxcharts
        h = boxchart(ones(length(y(lsidx)), 1)*ctr, y(lsidx), ...
            'BoxFaceColor', colors(ls, :), 'Notch', 'on'); hold on;
        h.JitterOutliers = 'on';
        h.MarkerStyle = '.';
        h.MarkerColor = colors(ls, :);
        ctr = ctr + 1;
    end
    
    % Format the box charts
    legend('off')
    set(gca, 'FontSize', 15);
    xlabel('Subject group');
    ylabel([label ' \DeltaR^2 ']);
    xlim([0.25 3.75]);
    ylim([-0.0 0.045]);
    
    set(gca, 'XTick', 1:3, 'XTickLabels', {'Spanish', 'English', 'Bilingual'}, ...
        'Ytick', -0.01:0.01:0.04,'FontSize', 15, 'Yscale', 'log');
    title(titles{f});
    box off;

    % Show reference lines
    h = yline(0, '--k');
    h.LineWidth = 0.8;
   

    % Add linear mixed effect model
    disp(['----------------- ' titles{f} ' ------------------------'])
    labels = {'MONO', 'ENG to BIL', 'SP to BIL'};
    combs = {[1, 2], [2, 4], [1, 4]};
    pos = {[1, 2], [2, 3], [1, 3]};
    
    % Perform linear mixed effect model and display significance
    for c = 1:3
        word_idx = ismember(wordsurp_encoding.ls, combs{c}) & ...
            wordsurp_encoding.(field)(:, index)>uv_thresh;
        tbl = wordsurp_encoding(word_idx, :);
        tbl.uv = tbl.(field)(:, index);
        lme2 = fitlme(tbl, 'uv~ls+(1|SID)+(1|el:SID)');
        P = lme2.Coefficients{2, 6};
        disp([labels{c} ': ' num2str(P)])

        [str, sig] = getSigStr(P, 2); 
        if ~strcmp(str, 'n.s.')
            line(pos{c}, [0.03 0.03]+c*0.0025, 'Color', 'k', 'LineWidth', 1.5);         
            text(mean(pos{c}), 0.031+c*0.0025, str, 'FontSize', 13);
        end
    end

    ctr = 1;
    
    % Create pie charts
    for ls = [1 2 4]
        lsidx = find(wordsurp_encoding.ls==ls & wordsurp_encoding.(field)>uv_thresh);

        subplot(3, 2, ctr*2);
        p = pie(sum([wordsurp_encoding.(field)(wordsurp_encoding.ls==ls, index)>uv_thresh, ...
            wordsurp_encoding.(field)(wordsurp_encoding.ls==ls, index)<uv_thresh]), ...
            [1, 1]); hold on;
        p(3).FaceColor = [0.7 0.7 0.7];
        p(3).EdgeColor = 'none';
        p(1).FaceColor = colors(ls, :);
        p(1).EdgeColor = 'none';

        ctr = ctr + 1;
    end    
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *elecs;

%% Visualize feature unique variance across native and unfamiliar language conditions

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
fieldname = {'English', 'Spanish'};
conditionLabels = {'native', 'unfamiliar'};
thresh = 0.001;

% label = 'word boundary';
% label = 'phonetic feature';

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% peakrate, phonetic features, word+surp
features = {'onset', 'peakrate', 'formant', 'consonant', 'word', 'surp'};
titles = {'onset', 'peakrate', 'formant', 'consonant', 'word', 'surp'};
% titles = {'peakrate', 'formant', 'consonant', 'word+surp'};
% titles = {'peakrate', 'phonetic features', 'word+surp', 'word', 'surp'};
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
subplts = {1, 2, 3, 4, 5, 6};
% subplts = {1, 2, 3, 4};
ynat_pie = nan(length(subplts), 2);
ynon_pie = nan(length(subplts), 2);

% split by sid
ynat_sid = cell(length(subplts), 2);
ynon_sid = cell(length(subplts), 2);

figure('Renderer', 'Painters');
for s = 1:length(subplts)
    
    % Combining native and unfamiliar subjects
    subplot(1, length(subplts), s)
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

    h=boxchart(ones(sum(y_nat>thresh), 1), y_nat(y_nat>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on'); hold on;
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

    h=boxchart(ones(sum(y_non>thresh), 1)*2, y_non(y_non>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on');
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([y_nat(y_nat>thresh); y_non(y_non>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    xticks([1 2]);
    xticklabels(conditionLabels);
    xlim([0.5 2.5])
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    idx = y_nat>thresh & y_non>thresh;
    tbl.rsq = [y_nat(idx); y_non(idx)];
    
    elecs = repmat(el{1}, rep, 1);
    sids = repmat(sid{1}, rep, 1);
    hemis = repmat(hemi{1}, rep, 1);

    tbl.sid = repmat(sids(idx), 2, 1);
    tbl.elec = repmat(elecs(idx), 2, 1);
    tbl.hemi = repmat(hemis(idx), 2, 1);
    tbl.native = [ones(sum(idx), 1); ones(sum(idx), 1)*2];

    clear elecs sids

    lme2 = fitlme(tbl,'rsq~native+(1|hemi)+(1|sid)+(1|elec:sid)');
    disp(titles{s})
    disp(lme2)
    p = lme2.Coefficients.pValue(2);

    disp(['LME native language p-value = ' num2str(p)])
    line([1, 2], [maxy(2)-0.1 maxy(2)-0.1], 'Color', 'k');
    if isempty(getSigStr(p, 2))
        text(1.4, maxy(2)-0.05, 'n.s.', 'FontSize', 15);
    else
        text(1.4, maxy(2)-0.05, getSigStr(p, 2), 'FontSize', 15);
    end
    title(titles{s});
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2')
end

figure; 
% pie plots showing number that meet threshold
for s = 1:length(subplts)
    y_pie = {ynat_pie(s, :), ynon_pie(s, :)};
    for native = 1:2 % 2 is the unfamiliar case
        % combining native and non-native subjects
        subplot(2, length(subplts), s + length(subplts)*(native-1))
        p = pie(y_pie{native}, [1, 1]); hold on;
        p(3).FaceColor = [0.3 0.3 0.3];
        p(3).EdgeColor = 'none';
        p(1).FaceColor = [0.8 0.8 0.8];
        p(1).EdgeColor = 'none';
        p(2).FontSize = 13;
        p(4).Color = 'w';
        p(4).FontSize = 13;
        title(conditionLabels{native});
    end
    clear y_pie
end 
lgd = legend({'No', 'Yes'});

% figure in which single subjects show percentage of electrodes meeting
% threshold
figure;
boxc = nan(2, 1);
numsubjs = length(unique(sid{1}));
for type = 1:length(titles)
    subplot(1, length(titles), type);
    jitter = rand(numsubjs, 1);
    y = ones(numsubjs, 1) - 0.2 + jitter*0.1;
    cols = getColorsCrossComp(1);

    % Find language background of subjects included in analysis
    lss = cellfun(@(x) wordsurp_encoding.ls(find(strcmp(wordsurp_encoding.SID, x), ...
        1, 'first')), unique(sid{type}, 'stable'));

    % Size indicates how many electrodes included in unique variance
    % analysis per subject
    sz = ynat_sid{type}(:, 2)./max(ynat_sid{type}(:, 2))*45;
    scatter(y, ynat_sid{type}(:, 1), sz, lss, 'filled', ...
        'MarkerFaceAlpha', 0.5, 'HandleVisibility', 'off'); hold on;
    
    b = boxchart(ones(numsubjs, 1)-0.15, ynat_sid{type}(:, 1)', ...
        'BoxFaceColor', 'none', 'MarkerStyle','none', 'BoxLineColor', 'k', ...
        'BoxWidth', 0.2);
    b.HandleVisibility = 'off';

    y = ones(numsubjs, 1) + 0.2 + jitter*0.1;
    sz = ynon_sid{type}(:, 2)./max(ynon_sid{type}(:, 2))*45;

    % Splitting so legend has english and spanish subjects
    for i = 1:2
        scatter(y(lss==i), ynon_sid{type}(lss==i, 1), sz(lss==i), ...
            lss(lss==i), 'filled', 'MarkerFaceAlpha', 0.5);
    end
    line([y-0.4 y]', [ynat_sid{type}(:, 1) ynon_sid{type}(:, 1)]', 'Color', ...
        [0.7 0.7 0.7], 'HandleVisibility', 'off');

    b = boxchart(ones(numsubjs, 1)+0.25, ynon_sid{type}(:, 1)', ...
        'BoxFaceColor', 'none', 'MarkerStyle','none', 'BoxLineColor', 'k', ...
        'BoxWidth', 0.2);
    b.HandleVisibility = 'off';

    % formatting
    ylim([0 1])
    yticks([0 0.5 1])
    xlim([0.6 1.5])
    xticks([1-0.25 1+0.25])
    xticklabels(conditionLabels)

    % TODO: Make this into an LME (?)
    [~, p]=ttest(ynat_sid{type}(:, 1), ynon_sid{type}(:, 1));
    disp(p)

    colormap([0.5 0.5 0.5; 0.8 0.8 0.8]); cols(1:2, :)
    legend({'Spanish mono', 'English mono'})

    ylabel(' \Delta R^2 > 0 (% of speech elecs)')
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* *elecs;

%% Same as above but monolingual vs. bilingual

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
fieldname = {'English', 'Spanish'};
conditionLabels = {'English', 'Spanish'};
thresh = 0.001 ;

% label = 'word boundary';
% label = 'phonetic feature';

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% peakrate, phonetic features, word+surp
features = {'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
titles = {'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
% titles = {'peakrate', 'formant', 'consonant', 'word+surp'};
% titles = {'peakrate', 'phonetic features', 'word+surp', 'word', 'surp'};
fields = {'eng_uv_all', 'sp_uv_all'}; 
ls=4;

ctr = 1;
engl = cell(length(features), 1);
sid_engl = cell(length(features), 1);
span = cell(length(features), 1);
sid_span = cell(length(features), 1);
sid = cell(length(features), 1);
el = cell(length(features), 1);
hemi = cell(length(features), 1);

for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);
    
    sid{ctr} = [sid{ctr}; wordsurp_encoding.SID(wordsurp_encoding.ls==ls)];
    el{ctr} = [el{ctr}; wordsurp_encoding.el(wordsurp_encoding.ls==ls)];
    hemi{ctr} = [hemi{ctr}; wordsurp_encoding.hemi(wordsurp_encoding.ls==ls)];

    for f = 1:2 % Aggregate over English and Spanish fields
        field = fields{f};    
        y = wordsurp_encoding.(field)(:, index);   
        sid_tmp = wordsurp_encoding.SID; 
        
        % Removing bilinguals so the comparison is more balanced
        if f == 1 
            engl{ctr} = [engl{ctr}; y(wordsurp_encoding.ls==ls)];  
            sid_engl{ctr} = [sid_engl{ctr}; sid_tmp(wordsurp_encoding.ls==ls)];
        else
            span{ctr} = [span{ctr}; y(wordsurp_encoding.ls==ls)];
            sid_span{ctr} = [sid_span{ctr}; sid_tmp(wordsurp_encoding.ls==ls)];
        end
    end
    ctr = ctr + 1;
end

% Combine the consonant and vowel features into one
subplts = {1, 2, 3, 4, 5, 6};
% subplts = {1, 2, 3, 4};
ynat_pie = nan(length(subplts), 2);
ynon_pie = nan(length(subplts), 2);

% split by sid
ynat_sid = cell(length(subplts), 2);
ynon_sid = cell(length(subplts), 2);

figure('Renderer', 'Painters');
for s = 1:length(subplts)
    
    % Combining native and unfamiliar subjects
    subplot(1, length(subplts), s)
    rep = length(subplts{s});
    
    % Native boxplot
    y_nat = arrayfun(@(x) engl{x}, [subplts{s}], 'UniformOutput', false);
    y_nat = cat(1, y_nat{:});

    % Per subject, how many electrodes meet threshold?
    y_sid = arrayfun(@(x) sid_engl{x}, [subplts{s}], 'UniformOutput', false);
    y_sid = repmat(y_sid{1}, rep, 1);

    prct = cellfun(@(x) sum(y_nat>thresh & strcmp(y_sid, x)) ...
        ./sum(strcmp(y_sid, x)), unique(y_sid));
    cnt = cellfun(@(x) sum(strcmp(y_sid, x)), unique(y_sid));
    ynat_sid{s, 1} = [prct, cnt];
    clear y_sid prct cnt

    h=boxchart(ones(sum(y_nat>thresh), 1), y_nat(y_nat>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on'); hold on;
                    h.JitterOutliers = 'on';
                    h.MarkerStyle = '.';
                    h.MarkerColor = 'k';

    % Unfamiliar boxplot
    y_non = arrayfun(@(x) span{x}, [subplts{s}], 'UniformOutput', false);
    y_non = cat(1, y_non{:});

    % Per sid, how many electrodes meet threshold?
    y_sid = arrayfun(@(x) sid_span{x}, [subplts{s}], 'UniformOutput', false);
    y_sid = repmat(y_sid{1}, rep, 1);

    prct = cellfun(@(x) sum(y_non>thresh & strcmp(y_sid, x)) ...
        ./sum(strcmp(y_sid, x)), unique(y_sid));
    cnt = cellfun(@(x) sum(strcmp(y_sid, x)), unique(y_sid));
    ynon_sid{s} = [prct, cnt];
    clear y_sid prct cnt

    h=boxchart(ones(sum(y_non>thresh), 1)*2, y_non(y_non>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on');
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([y_nat(y_nat>thresh); y_non(y_non>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    xticks([1 2]);
    xticklabels(conditionLabels);
    xlim([0.5 2.5])
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    idx = y_nat>thresh & y_non>thresh;
    tbl.rsq = [y_nat(idx); y_non(idx)];
    
    elecs = repmat(el{1}, rep, 1);
    sids = repmat(sid{1}, rep, 1);
    hemis = repmat(hemi{1}, rep, 1);

    tbl.sid = repmat(sids(idx), 2, 1);
    tbl.elec = repmat(elecs(idx), 2, 1);
    tbl.hemi = repmat(hemis(idx), 2, 1);
    tbl.native = [ones(sum(idx), 1); ones(sum(idx), 1)*2];

    clear elecs sids

    lme2 = fitlme(tbl,'rsq~native+(1|hemi)+(1|sid)+(1|elec:sid)');
    disp(titles{s})
    disp(lme2)
    p = lme2.Coefficients.pValue(2);

    disp(['LME L1 language p-value = ' num2str(p)])
    line([1, 2], [maxy(2)-0.1 maxy(2)-0.1], 'Color', 'k');
    if isempty(getSigStr(p, 2))
        text(1.4, maxy(2)-0.05, 'n.s.', 'FontSize', 15);
    else
        text(1.4, maxy(2)-0.05, getSigStr(p, 2), 'FontSize', 15);
    end
    title(titles{s});
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2')
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* *elecs;

%% Same as above but really monolingual bilingua

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
fieldname = {'monolingual', 'bilingual'};
conditionLabels = {'monolingual', 'bilingual'};
thresh = 0.001 ;

% label = 'word boundary';
% label = 'phonetic feature';

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% peakrate, phonetic features, word+surp
features = {'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
titles = {'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
% titles = {'peakrate', 'formant', 'consonant', 'word+surp'};
% titles = {'peakrate', 'phonetic features', 'word+surp', 'word', 'surp'};
fields = {'eng_uv_all', 'sp_uv_all'}; 

ctr = 1;
mono = cell(length(features), 1);
sid_mono = cell(length(features), 1);
bil = cell(length(features), 1);
sid_bil = cell(length(features), 1);
sid = cell(length(features), 1);
el = cell(length(features), 1);
hemi = cell(length(features), 1);

[proficiency] = getBilingProf();
profield = {'eng_prof', 'span_prof'};

for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);

    for f = 1:2 % Aggregate over English and Spanish fields
        field = fields{f};    
        y = wordsurp_encoding.(field)(:, index);   
        sid_tmp = wordsurp_encoding.SID; 
        
        % Removing bilinguals so the comparison is more balanced
        if f == 1 
            mono{ctr} = [mono{ctr}; y(ismember(wordsurp_encoding.ls, 2))];  
            sid{ctr} = [sid{ctr}; sid_tmp(ismember(wordsurp_encoding.ls, 2))];
            el{ctr} = [el{ctr}; wordsurp_encoding.el(wordsurp_encoding.ls==2)];
            hemi{ctr} = [hemi{ctr}; wordsurp_encoding.hemi(wordsurp_encoding.ls==2)];
        else
            mono{ctr} = [mono{ctr}; y(ismember(wordsurp_encoding.ls, 1))];  
            sid{ctr} = [sid{ctr}; sid_tmp(ismember(wordsurp_encoding.ls, 1))];
            el{ctr} = [el{ctr}; wordsurp_encoding.el(wordsurp_encoding.ls==1)];
            hemi{ctr} = [hemi{ctr}; wordsurp_encoding.hemi(wordsurp_encoding.ls==1)];
        end

        % subset to bilinguals with full proficiency    
        profSIDs = [proficiency.SID(proficiency.(profield{f})==5)];
        fullprof_biling = wordsurp_encoding.ls==4 & ...
            ismember(wordsurp_encoding.SID, profSIDs);
        bil{ctr} = [bil{ctr}; y(fullprof_biling)];
        sid{ctr} = [sid{ctr}; sid_tmp(fullprof_biling)];
        el{ctr} = [el{ctr}; wordsurp_encoding.el(fullprof_biling)];
        hemi{ctr} = [hemi{ctr}; wordsurp_encoding.hemi(fullprof_biling)];

    end
    ctr = ctr + 1;
end

% Combine the consonant and vowel features into one
subplts = {1, 2, 3, 4, 5, 6};

% split by sid
ynat_sid = cell(length(subplts), 2);
ynon_sid = cell(length(subplts), 2);

figure('Renderer', 'Painters');
for s = 1:length(subplts)
    
    % Combining native and unfamiliar subjects
    subplot(1, length(subplts), s)
    rep = length(subplts{s});
    
    % Native boxplot
    y_nat = arrayfun(@(x)  mono{x}, [subplts{s}], 'UniformOutput', false);
    y_nat = cat(1, y_nat{:});

    h=boxchart(ones(sum(y_nat>thresh), 1), y_nat(y_nat>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on'); hold on;
                    h.JitterOutliers = 'on';
                    h.MarkerStyle = '.';
                    h.MarkerColor = 'k';

    % Unfamiliar boxplot
    y_non = arrayfun(@(x) bil{x}, [subplts{s}], 'UniformOutput', false);
    y_non = cat(1, y_non{:});

    h=boxchart(ones(sum(y_non>thresh), 1)*2, y_non(y_non>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on');
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([y_nat(y_nat>thresh); y_non(y_non>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    xticks([1 2]);
    xticklabels(conditionLabels);
    xlim([0.5 2.5])
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    tbl.rsq = [y_nat; y_non];
    
    elecs = el{1};
    sids = sid{1};
    hemis = hemi{1};

    tbl.sid = sids;
    tbl.elec = elecs;
    tbl.hemi = hemis;
    tbl.native = [ones(length(y_nat), 1); ones(length(y_non), 1)*2];

    clear elecs sids

    lme2 = fitlme(tbl(tbl.rsq>0, :),'rsq~native+(1|hemi)+(1|sid)+(1|elec:sid)');

    disp(titles{s})
    disp(lme2)
    p = lme2.Coefficients.pValue(2);

    disp(['LME L1 language p-value = ' num2str(p)])
%     line([1, 2], [maxy(2)-0.1 maxy(2)-0.1], 'Color', 'k');
    if isempty(getSigStr(p, 2))
        text(1.4, maxy(2)-0.05, 'n.s.', 'FontSize', 15);
    else
        text(1.4, maxy(2)-0.05, getSigStr(p, 2), 'FontSize', 15);
    end
    title(titles{s});
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2')
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* *elecs;

%% Visualize MNI location of electrodes that meet the unique variance threshold

% initialize design electrode structure
fields = {'sp_uv_all', 'eng_uv_all'}; 

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% uv feature order
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant',

figure();
for lang = 1:2
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
        desel.sz = [2; 25*desel.conds']; 
        
        % split up peak rate and phonetic features again for MNI plotting
        desel.labels = [];
        desel.yval = arrayfun(@(x) wordsurp_encoding.(fields{lang})(x, index), ...
            1:height(wordsurp_encoding));
        
        % discretize values
        binedges = [0:0.005:0.01 0.015:0.015:0.045];
        % binedges = [-1, binedges];    
        for s=unique(wordsurp_encoding.SID)'
            SID = s{1};
            idx = strcmp(wordsurp_encoding.SID, SID);
            desel.(SID).elid = wordsurp_encoding.el(idx);
            desel.(SID).condition = discretize(desel.yval(idx), ...
                binedges);
        end
        
        desel.cols = [1 1 1; [linspace(1, featcol(1), length(binedges)); ...
            linspace(1, featcol(2), length(binedges)); ...
            linspace(1, featcol(3), length(binedges))]'];
        if strcmp(feat, 'word+surp') ||  strcmp(feat, 'surp')
            gns = flipud(fpurple(length(binedges)-2));
            desel.cols = [1 1 1 ; gns];
            desel.cols = [1 1 1 ; repmat(gns(3, :),length(binedges), 1) ];
    
            % colorbar
            figure;
            colormap(desel.cols);
            colorbar;            
        end
        
        lsid = find(ismember(wordsurp_encoding.ls, ls));
        mni_lh = plotMNIElec(unique(wordsurp_encoding.SID(lsid)), desel, 'lh', 0);
        sgtitle(lang);
        l = light;
        view(270, 0);   
        set(l,'Style', 'infinite', 'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.6;
        % add a pie
        axes('Position',[.7 .1 .3 .3])
        p = pie([sum(mni_lh.cond>1), sum(mni_lh.cond==1)], [1 1]); 
        p(1).FaceColor = [desel.cols(2, :)];
        p(1).EdgeColor = 'none';
        p(3).FaceColor = [0.6 0.6 0.6];
        p(3).EdgeColor = 'none';
        p(2).FontWeight = 'bold';
        p(2).FontSize = 13;
        p(4).FontWeight = 'bold';
        p(4).FontSize = 13;

        mni_rh = plotMNIElec(unique(wordsurp_encoding.SID(lsid)), desel, 'rh', 0);
        sgtitle(lang);
        l = light;
        view(90, 0);
        set(l,'Style', 'infinite', 'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.6;
        % add a pie
        axes('Position',[.7 .1 .3 .3])
        p = pie([sum(mni_rh.cond>1), sum(mni_rh.cond==1)], [1 1]); 
        p(1).FaceColor = [desel.cols(2, :)];
        p(1).EdgeColor = 'none';
        p(3).FaceColor = [0.6 0.6 0.6];
        p(3).EdgeColor = 'none';
        p(2).FontWeight = 'bold';
        p(2).FontSize = 13;
        p(4).FontWeight = 'bold';
        p(4).FontSize = 13;
     end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd *elecs;

%% Visualize MNI location of electrodes for NATIVE vs. UNFAMILIAR
% initialize design electrode structure
fields = {'sp_uv_all', 'eng_uv_all'}; 

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
        desel.sz = [2; 10; 25*desel.conds']; 
        
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

%           % English natives
            eidx = find(wordsurp_encoding.ls==2);
            desel.yval(eidx) = arrayfun(@(x) wordsurp_encoding.(fields{1})(x, index), eidx);
        end
        
        % discretize values
        % 0.005:0.005:0.01
        binedges = [0.0005:0.005:0.01 0.015:0.015:0.1];
        binedges = [-1, binedges];    
        for s=unique(wordsurp_encoding.SID)'
            SID = s{1};
            idx = strcmp(wordsurp_encoding.SID, SID);
            desel.(SID).elid = wordsurp_encoding.el(idx);
            desel.(SID).condition = discretize(desel.yval(idx), ...
                binedges);
        end
        
        desel.cols = [1 1 1; [linspace(1, featcol(1), length(binedges)); ...
            linspace(1, featcol(2), length(binedges)); ...
            linspace(1, featcol(3), length(binedges))]'];
%         if strcmp(feat, 'word+surp') ||  strcmp(feat, 'surp')
% %             gns = flipud(fpurple(length(binedges)-2));
%             gns = flipud(fpurple(4));
%             desel.cols = [1 1 1 ; gns];
%             desel.cols = [1 1 1 ; repmat(gns(3, :),length(binedges), 1) ];
%         end
        if native
            cls = flipud(blues(8));
        else
            cls = flipud(reds(8));
        end
        desel.cols = [1 1 1 ;cls(3:end, :)];

        figure;
        colormap(desel.cols);
        colorbar;
        
        mni_lh = plotMNIElec(unique(wordsurp_encoding.SID), desel, 'lh', 0);
        sgtitle(native);
        l = light;
        view(270, 0);   
        set(l,'Style', 'infinite', 'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.6;
        % add a pie
        axes('Position',[.6 .15 .3 .3])
        p = pie([sum(mni_lh.cond>1), sum(mni_lh.cond==1)], [1 1]); 
        p(1).FaceColor = [desel.cols(5, :)];
        p(1).EdgeColor = 'none';
        p(3).FaceColor = [0.6 0.6 0.6];
        p(3).EdgeColor = 'none';
        p(2).Color = 'w';
        p(2).FontWeight = 'bold';
        p(2).FontSize = 13;
        p(4).FontWeight = 'bold';
        p(4).Color = 'w';
        p(4).FontSize = 13;

        mni_rh = plotMNIElec(unique(wordsurp_encoding.SID), desel, 'rh', 0);
        sgtitle(native);
        l = light;
        view(90, 0);
        set(l,'Style', 'infinite', 'Position',[1 0 0],'Color',[0.8 0.8 0.8]);
        alpha 0.6;
        % add a pie
        axes('Position',[.6 .15 .3 .3])
        p = pie([sum(mni_rh.cond>1), sum(mni_rh.cond==1)], [1 1]); 
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
    betaInfo* *encoding* *wrd *elecs;

%% Visualize native brain location of electrodes ''

% initialize design electrode structure
fieldnames = {'Spanish', 'English'};
fields = {'sp_uv_all', 'eng_uv_all', '', 'sp_uv_all'}; 

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% uv feature order
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant',

fig = figure();
for lang = 1:2
    for f = 1:length(feats)
        feat = feats{f};
        index = find(ismember(wordsurp_details.featureOrd, feat));
    
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
    
        % Make the desel structure
        desel=struct();
        desel.conds = 1:7;
        
        % Determine size and color
        desel.sz = [2; 35*ones(length(desel.conds), 1)];
        desel.sz = [5; 25*desel.conds']; 
        
        % Split peak rate and phonetic features for MNI plotting wordsurp_encoding.ls(x)
        desel.labels = [];
        desel.yval = arrayfun(@(x) wordsurp_encoding.(fields{lang})(x, index), ...
            1:height(wordsurp_encoding));
        
        % Determine bin-edges linearly
        % yvals = sort(desel.yval(desel.yval>0 & ismember(wordsurp_encoding.ls,ls)'));
        % binedges = yvals(1:ceil(length(yvals)/8):length(yvals));
        % [~, binedges] = discretize(desel.yval(desel.yval>0 & ismember(wordsurp_encoding.ls,ls)'), ...
        %     length(desel.conds)-1);

        % Construct manual, non-linear edges
%         binedges = [-1, 0.00:0.005:0.01 0.015:0.015:0.045];
        binedges = [0.0005:0.005:0.01 0.015:0.015:0.1];
        binedges = [-1, binedges]; 

        for s=unique(wordsurp_encoding.SID)'
            SID = s{1};
            idx = strcmp(wordsurp_encoding.SID, SID);
            desel.(SID).elid = wordsurp_encoding.el(idx);
            desel.(SID).condition = discretize(desel.yval(idx), ...
                binedges);
        end 
        
        desel.cols = [1 1 1; [linspace(1, featcol(1), length(binedges)); ...
            linspace(1, featcol(2), length(binedges)); ...
            linspace(1, featcol(3), length(binedges))]'];
        if strcmp(feat, 'word+surp')
                
            gns = flipud(fpurple(length(binedges)-2));
            desel.cols = [1 1 1 ; gns];
            desel.cols = [0 0 0; repmat(gns(3, :),length(binedges), 1) ];
    
            % colorbar
            figure;
            colormap(desel.cols);
            colorbar;            
        end

        SIDs = {'EC183'};
        desel.(SIDs{1}).selid= 71; %135 
        % EC163 - el55
        ls = 2; % can only do 1-2
        %desel.(SIDs{1}).selid=21; % for EC100

        if lang == ls   
            cls = flipud(blues(8));
        else
            cls = flipud(reds(8));
        end
        desel.cols = [0 0 0;cls(3:end, :)];
         
%         lsid = find(ismember(wordsurp_encoding.ls, ls));
%         SIDs = unique(wordsurp_encoding.SID(lsid));
%         SIDs = {'EC183'};
%         desel.(SIDs{1}).selid=135;

%         SIDs = {'EC183'};
%         desel.(SIDs{1}).selid=71; %135

%           SIDs = {'EC172'};
            
% 
%         SIDs = {'EC260'};
%         desel.(SIDs{1}).selid=71; %135

        % just for native brain and coverage
%         desel.(SIDs{1}).elid=[];
%         desel.(SIDs{1}).condition=[];
        
    
        % desel.cols = [1 1 1; 212/256, 228/256 188/256; ...
        %     54/256, 85/256, 183/256; 64/256 55/256 110/256]; 
        nh = plotNativeElec(SIDs, desel, 1);

        % only works if its on one subject
        l = light; 
        if strcmp(imgall.(SIDs{1}).hemi, 'lh')  
            view(270, 0);   
            set(l,'Style', 'infinite', 'Position',[-1 0 0],'Color',[0.8 0.8 0.8]);
        else
            view(90, 0);
            set(l,'Style', 'infinite', 'Position',[1 0 1],'Color',[0.8 0.8 0.8]);
        end

        %figure('Renderer', 'Painters');
        
        alpha 0.8;
        % add a pie
%         axes('Position',[.7 .1 .3 .3])
%         p = pie([sum(nh.cond>1), sum(nh.cond==1)], [1 1]); 
%         p(1).FaceColor = [desel.cols(2, :)];
%         p(1).EdgeColor = 'none';
%         p(2).FontSize = 15;
%         p(3).FaceColor = [0.6 0.6 0.6];
%         p(3).EdgeColor = 'none';
%         p(4).FontSize = 15;
%         title(fieldnames{lang})

        axes('Position',[.6 .15 .3 .3])
        p = pie([sum(nh.cond>1), sum(nh.cond==1)], [1 1]); 
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
    betaInfo* *encoding* *wrd *elecs;

%% Visualize the correlation of word + surprise unique variance to other features

% Extract unique variance of word+surp from the wordsurp table

% Top row of figure is English (TIMIT), bottom row is Spanish (DIMEx)
fields = {'eng_uv_all', 'sp_uv_all'};

labels = {{'word+surp', 'onset'}, {'word+surp', 'peakrate'}, ...
    {'word+surp', 'formant'}, {'word+surp', 'consonant'}, {'word', 'surp'}};

figure('Position', [200, 100, 900, 900]);
ctr=1;
for field = fields

    % Select only native speakers of the language presented
    if strcmp(field, 'eng_uv_all')
        lsidx = ismember(wordsurp_encoding.ls, [2, 4]);
    else
        lsidx = ismember(wordsurp_encoding.ls, [1, 4]);
    end

    for l = labels
        label = l{1};
        index = [];
        index(1) = find(ismember(wordsurp_details.featureOrd, label{1}));   
        index(2) = find(ismember(wordsurp_details.featureOrd, label(2)));   
        
        ax(ctr)=subplot(2, length(labels), ctr);
       
        x = wordsurp_encoding.(field{1})(lsidx, index(2));
        y = wordsurp_encoding.(field{1})(lsidx, index(1));
%         sid = cellfun(@(x) str2double(x(3:end)), ...
%             wordsurp_encoding.SID(wordsurp_encoding.ls==ls));

        % Ensure positive unique variance for both features being compared
        pos = x>0 & y>0;
        uv_thresh = 0.001;

        [r, p] = corr(x(pos), y(pos), 'Type', 'Spearman'); 
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
        ylabel([label{1} ' \Delta R^2']);
        title({['r= ' num2str(r, 3) ',']},{['p=' num2str(p, 3)]});

        set(gca, 'FontSize', 13);
        ctr=ctr+1;
    end
end
% linkaxes(ax(1:3));   % phonetic feat
% linkaxes(ax(4:6));   % word onset/surp 

for i=1:length(labels)*2
    subplot(2, length(labels), i);
%     h=refline(1, 0);
    h = lsline;
    h.LineWidth = 2;
    h.Color = 'k';
    xline(0, 'Color', 'k', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% plot in alternate form (above)

% Extract unique variance of word+surp from the wordsurp table
eng_uv_all = wordsurp_encoding.eng_uv_all;
sp_uv_all = wordsurp_encoding.sp_uv_all;

% Top row of figure is English (TIMIT), bottom row is Spanish (DIMEx)
fields = {'sp_uv_all', 'eng_uv_all'};
featureOrd = {'onset', 'peakrate', 'formant', 'consonant', ...
            'word+surp', 'word', 'surp'};

labels = {{'word+surp', 'onset'}, {'word+surp', 'peakrate'}, ...
    {'word+surp', 'formant'}, {'word+surp', 'consonant'}}; % , {'word', 'surp'}

r = nan(2, 4);
p = nan(2, 4);
for l = labels

    corp=1;
    label = l{1};
    index = [];
    index(1) = find(ismember(wordsurp_details.featureOrd, label{1}));   
    index(2) = find(ismember(wordsurp_details.featureOrd, label{2}));   

    for field = fields
        % Select only native speakers of the language presented
        if strcmp(field, 'eng_uv_all')
            lsidx = ismember(wordsurp_encoding.ls, [2, 4]);
        else
            lsidx = ismember(wordsurp_encoding.ls, [1, 4]);
        end
       
        x = wordsurp_encoding.(field{1})(lsidx, index(2));
        y = wordsurp_encoding.(field{1})(lsidx, index(1));

        % Ensure positive unique variance for both features being compared
        pos = x>0 & y>0;
        uv_thresh = 0.001;

        [r(corp, index(2)), p(corp, index(2))] = corr(x(pos), y(pos), ...
            'Type', 'Spearman'); 
        corp=corp+1;
    end
end

figure;
% compare TIMIT and DIMEx
bar(r', 'grouped', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
for corp = 1:2
    for i = 1:4
        text(i+0.2*(corp-1), r(corp, i)+0.1, ...
            getSigStr(p(corp, i), 2), 'FontSize', 13);
    end
end

ylim([-0.1 0.6]);
box off;
set(gca, 'FontSize', 15, 'Xtick', 1:4, 'Xticklabel', ...
    featureOrd(1:4), 'Ytick', [0 0.5]);
ylabel({'Spearman corr with', 'word+surp \Delta R^2'})
legend({'Spanish speech', 'English speech'})


%% Compare prediction performance in cross-language tested TRF models

% Repeated sentences in DIMEx and TIMIT
repsentName = {'s00104','s00804', 's01904', 's03004', 's05004', ...
    's06104', 's06804', 's07404', 's09004', ...
    'fcaj0_si1479', 'fcaj0_si1804', 'fdfb0_si1948', 'fdxw0_si2141', ...
    'fisb0_si2209', 'mbbr0_si2315', 'mdlc2_si2244', 'mdls0_si998', ...
    'mjdh0_si1984', 'mjmm0_si625'};

modelname='onset_phnfeatConsOnset_maxDtL_formantMedOnset'; 
SIDs = [sSIDs, eSIDs bSIDs];
% SIDs(ismember(SIDs, {'EC152', 'EC225', 'EC252', 'EC161', 'EC200', ...
%     'EC212', 'EC221'})) = [];

% Subjects with no repeated sentences for DIMEx so unable to calculate
% cross-language predictions
SIDs(ismember(SIDs, {'EC252', 'EC152'})) = [];
samelang = 'timit';

if strcmp(samelang, 'dimex')
    predcols = {'b', 'r'}; 
    crosslang = 'timit';
    sentdet = dimex_details.sentdet(ismember({dimex_details.sentdet.name}, ...
        repsentName));
elseif strcmp(samelang, 'timit')
    predcols = {'r', 'b'};  
    crosslang = 'dimex';
    sentdet = timit_details.sentdet(ismember({timit_details.sentdet.name}, ...
        repsentName));
end

% Initialize all the correlation structures
samelang_corr = cell(length(SIDs), 1);
crosslang_corr = cell(length(SIDs), 1);
pred_corr = cell(length(SIDs), 1);

testR_same = nan(height(wordsurp_encoding), 1);
testR_cross = nan(height(wordsurp_encoding), 1);

ctr = 1;
for s = SIDs
    SID = s{1};

    idx = strcmp(wordsurp_encoding.SID, SID);

    % Calculating same language predictions
    [out_same, testR] = out_addStrfpred(SID, samelang, modelname, 1, ...
        sentdet, samelang);
    testR_same(idx) = testR(wordsurp_encoding.el(idx));
   
    % Calculating cross language predictions
    [out_cross, testR] = out_addStrfpred(SID, samelang, modelname, 1, ...
        sentdet, crosslang);
    testR_cross(idx) = testR(wordsurp_encoding.el(idx));

%     % Correlate predictions from same and cross language tested STRF models
%     meanrepresp = cellfun(@(x) mean(x, 3), {out_same.resp}, 'UniformOutput', false);
%     resp = mean(cat(3, [meanrepresp{:}]), 3);
%     same_pred = cat(3, [out_same.predResp]);
%     cross_pred = cat(3, [out_cross.predResp]);
% 
%     % Prediction to response and prediction to prediction correlations
%     samelang_corr{ctr} = diag(corr(same_pred', resp'));
%     crosslang_corr{ctr} = diag(corr(cross_pred', resp'));
%     pred_corr{ctr} = diag(corr(same_pred', cross_pred'));

    % Plot several example sentence predictions
    if strcmp(SID, 'EC100')
        % [~, maxels] = max(mean(cat(3, testR_cross(1:minel, :), 
        % testR_same(1:minel, :)), [2, 3],'omitnan'));

        % Use the same electrodes from above for cross-prediction example
        maxels = [22 150];
        
        % Plot each electrode example sentences as a new figure
        for maxel = maxels
            figure('renderer', 'painters');

            for sent = 1:3          
                subplot(1, 3, sent);                            

                % Smooth neural response 
                data = smoothdata(squeeze(out_same(sent).resp(maxel, :, :)), ...
                    'gaussian', 'SmoothingFactor', 0.1);
                x = -0.49:0.01:0.01*size(out_same(sent).resp, 2)-0.5;        
                
                % Plot true neural response
                if size(data, 1) == 1
                    plot(x, data, 'LineWidth', 1.75, 'Color', predcols{1});
                    hold on;
                else
                    data = data - mean(data(1:55, :), [1, 2]);
                    addpath(genpath('../../../ecog_scripts'));
                    shadedErrorBar(x, data', {@mean,@nansem}, ...
                        {'color', predcols{1}, 'linewidth', 1.75, ...
                        'DisplayName', 'response'}, 0.3);
                    hold on;
                end
                
                % Plot predicted responses
                yyaxis right
                plot(x, out_same(sent).predResp(maxel, :, :), 'LineStyle', '--', ...
                    'LineWidth', 2, 'Color', predcols{1}, 'DisplayName', 'same-prediction');
                plot(x, out_cross(sent).predResp(maxel, :, :), 'LineStyle', '--', ...
                    'LineWidth', 2, 'Color', predcols{2}, 'DisplayName', 'cross-prediction');
                xlim([0.2 1.5]);

                yticks([]);
                sentidx = strcmp({sentdet(:).name},  out_cross(sent).name);     
                title(join(sentdet(sentidx).wordList, ' '));
            end                                                 
        end
    end
    ctr = ctr + 1;
end
wordsurp_encoding.([samelang '_trained']) = [testR_same, testR_cross]; 

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Scatter plot of mean same-language and cross-langauge predictions on test set

ctr=1;
figure;
for ls = [1, 2, 4]
    subplot(1, 3, ctr);

    % Aggregate all same language and cross-language predictions
    lsidx = find(wordsurp_encoding.ls==ls);
    testR_Same = mean([wordsurp_encoding.dimex_trained(:, 1), ...
        wordsurp_encoding.timit_trained(:, 1)], 2, 'omitnan'); 
    testR_Cross = mean([wordsurp_encoding.dimex_trained(:, 2), ...
        wordsurp_encoding.timit_trained(:, 2)], 2, 'omitnan'); 

    % Scatter point of R^2 values
    scatter(testR_Same(lsidx).^2, testR_Cross(lsidx).^2, 15, ...
        [0.6 0.6 0.6], 'filled');
    hold on;

    % Calculate correlation coefficient
    [rho, ~] = corr(testR_Same(lsidx).^2, testR_Cross(lsidx).^2, ...
        'rows', 'pairwise');
    text(0.15, 0.7, ['r=' num2str(rho, 2)], 'FontSize', 15);

    % Formatting
    ylim([0.1 0.75]);
    xlim([0.1 0.75]);   
    xlabel({'Mean same-language','prediction R^2'});
    h = refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';
    set(gca, 'FontSize', 15);

    if ctr==1
        ylabel({'Mean cross-language','prediction R^2'});
    end

    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% --------------------------- all subjects -------------------------------
%% ---------------- OLD: STRF surprisal unique variance -------------------
modelnames_timit={'onset_phnfeatonset_maxDtL', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpNoOnsBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_engAvgBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_phnPos_engNormSurpBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_phnPos', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_engTriSurpBin', ...
    'onset_phnfeatonset', ...
    'onset_maxDtL', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_engGPTPrevBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpBin_engCrossTriSurpBin', ...
    };  
% engSurpBin_engEntBin, engTriSurpBin, engTriSurpBin_engTriEntBin

modelnames_dimex={'onset_phnfeatonset_maxDtL', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpNoOnsBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_spAvgBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_phnPos_spNormSurpBin', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_phnPos', ...
    'onset_phnfeatonset_maxDtL_wordOns_wordL_spTriSurpBin', ...
    'onset_phnfeatonset', 'onset_maxDtL'}; 

modelnames = {modelnames_timit, modelnames_dimex};

labels = {'word', 'surprisal', 'surpword' 'surpaverage', 'phnpos', 'normsurp', ...
    'biphn', 'triphn', 'GPT', 'peakrate', 'phnfeat'};
lang = {'eng', 'sp'};
corpus = {'timit', 'dimex'};

timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");


% determine unique variance per feature and primary encoding
varnames = {'SID', 'el', 'ls', 'eng_base_rsq', ...
    ['eng_rsq_' labels{1}], ...
    ['eng_rsq_' labels{2}], ...
    ['eng_uv_' labels{1}], ...
    ['eng_uv_' labels{2}], ...
    ['engcross_uv_' labels{2}], ...
    ['eng_uv_' labels{3}], ...
    ['eng_uv_' labels{4}], ...
    ['eng_uv_' labels{5}], ...
    ['eng_uv_' labels{7}], ...
    ['eng_uv_' labels{8}], ...
    ['eng_uv_' labels{9}], ...
    ['eng_uv_' labels{10}], ...
    ['eng_uv_' labels{11}], ...
    'sp_base_rsq', ...
    ['sp_rsq_' labels{1}], ...
    ['sp_rsq_' labels{2}], ...
    ['sp_uv_' labels{1}], ...
    ['sp_uv_' labels{2}], ...
    ['sp_uv_' labels{3}], ...
    ['sp_uv_' labels{4}], ...
    ['sp_uv_' labels{5}], ...
    ['sp_uv_' labels{7}], ...
    ['sp_uv_' labels{8}], ...
    ['sp_uv_' labels{10}], ...
    ['sp_uv_' labels{11}], ...
    };
surp_encoding =  array2table(zeros(0, length(varnames)), 'VariableNames', varnames);

[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
for s = [sSIDs eSIDs bSIDs mSIDs]
    SID = s{1}; 
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs})); %  mSIDs
    for l = 1:2
        corpusStrf{l} = loadMultModelStrf(SID, modelnames{l}, corpus{l}, ...
            datapath, 1, 'v5');
    end

    if ~any(cellfun(@(x) isempty(x), [corpusStrf{1} corpusStrf{2}]))

        % find which electrodes meet threshold rsq for the full bin model
        minel = min(cellfun(@(x) length(x{3}.meanTestR), corpusStrf));
%         els = find(corpusStrf{1}{3}.meanTestR(1:minel).^2' > 0.1 | ...
%             corpusStrf{2}{3}.meanTestR(1:minel).^2' > 0.1);
        if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            els = intersect(timit_elecs.allidx.(SID), dimex_elecs.allidx.(SID));
        elseif isfield(dimex_elecs.allidx, SID)
            els = dimex_elecs.allidx.(SID);
        elseif isfield(timit_elecs.allidx, SID)
            els = timit_elecs.allidx.(SID);
        end

        base = cell(2, 1);
        surp = cell(2, 1);
        surpBin = cell(2, 1);      
        avgsurpBin = cell(2, 1);
        phnpos = cell(2, 1);
        normsurpBin = cell(2, 1);
        biphnBin = cell(2, 1);
        triphnBin = cell(2, 1);

        % timit only
        surpCrossBin = cell(1, 1);
        surpGPT = cell(1, 1);

        uvsurp = cell(2, 1);
        uvpeakrate = cell(2, 1);
        uvsurpword = cell(2, 1);
        uvsurpbin = cell(2, 1);
        uvsurpbincross = cell(2, 1);
        uvavgsurpbin = cell(2, 1);
        uvphnpos = cell(2, 1);
        uvnormsurp = cell(2, 1);
        uvtriphn = cell(2, 1);

        for l = 1:2
            % models without pitch
            base{l} = (corpusStrf{l}{1}.meanTestR.^2)';
            surp{l} = (corpusStrf{l}{2}.meanTestR.^2)';
            surpBin{l} = (corpusStrf{l}{3}.meanTestR.^2)';
            avgsurpBin{l} = (corpusStrf{l}{4}.meanTestR.^2)'; 
            normsurpBin{l} = (corpusStrf{l}{5}.meanTestR.^2)'; 
            phnpos{l} = (corpusStrf{l}{6}.meanTestR.^2)'; 
            triphnBin{l} = (corpusStrf{l}{7}.meanTestR.^2)'; 
            remPr{l} = (corpusStrf{l}{8}.meanTestR.^2)'; 
            remPhn{l} = (corpusStrf{l}{9}.meanTestR.^2)'; 
    
            % unique variances
            uvsurp{l} = surp{l} - base{l};
            uvsurpword{l} = surpBin{l} - base{l};
            uvsurpbin{l} = surpBin{l} - surp{l};
            uvavgsurpbin{l} = avgsurpBin{l} - surp{l};
            uvphnpos{l} = phnpos{l} - surp{l};
            uvnormsurp{l} = normsurpBin{l} - surp{l};
            uvtriphn{l} = triphnBin{l} - surp{l};
            uvpeakrate{l} = base{l} - remPr{l};
            uvphnfeat{l} = base{l} - remPhn{l};

            % cross language case
            if strcmp(corpus{l}, 'timit')
                surpGPT{l} = (corpusStrf{l}{9}.meanTestR.^2)';     
                uvGPT{l}= surpGPT{l} - surp{l};

                surpCrossBin{l} = (corpusStrf{l}{10}.meanTestR.^2)';     
                uvsurpbincross{l}= surpCrossBin{l} - surp{l};
            end
        end

        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        tmp = table(sids, els, lss, ...
            base{1}(els), surp{1}(els), surpBin{1}(els),  ...
            uvsurp{1}(els), uvsurpbin{1}(els), uvsurpbincross{1}(els), ...
            uvsurpword{1}(els), uvavgsurpbin{1}(els), ...
            uvphnpos{1}(els), uvnormsurp{1}(els), ...
            uvtriphn{1}(els), uvGPT{1}(els), uvpeakrate{1}(els), uvphnfeat{1}(els), ...
            base{2}(els), surp{2}(els), surpBin{2}(els), ...
            uvsurp{2}(els), uvsurpbin{2}(els), uvsurpword{2}(els), ...
            uvavgsurpbin{2}(els), ...
            uvphnpos{2}(els), uvnormsurp{2}(els), ...
            uvtriphn{2}(els), uvpeakrate{2}(els), uvphnfeat{2}(els), ...
            'VariableNames', varnames);

%         tmp = table(sids, els, lss, ...
%             base{1}(els), surp{1}(els), surpBin{1}(els), ...
%             uvsurp{1}(els), uvsurpbin{1}(els), ...
%             uvsurpbincross{1}(els), uvavgsurpbin{1}(els), ...
%             nan(length(els), 1), nan(length(els), 1), nan(length(els), 1), ...
%             nan(length(els), 1), nan(length(els), 1), nan(length(els), 1), ...
%             'VariableNames', varnames);
        surp_encoding = [surp_encoding; tmp];
    else
        warning(['Missing subject ' SID]);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *wrd;

%% vis: surprisal vs. word unique variance

cols = getColorsCrossComp(1);

figure;
% maxls = 5;
% ax = subplot(2, 3, 1);
ctr = 1;
ax = [];
for ls = [1, 2, 4]
    

    label = 'surprisal';
    fields_surp = {'eng_uv_surprisal', 'sp_uv_surprisal'};
    fields_word = {'eng_uv_word', 'sp_uv_word'};

    for lang = 1:2 % english, spanish
         % english
        uv_surp = surp_encoding.(fields_surp{lang});
        uv_word = surp_encoding.(fields_word{lang});
    
        ax = [ax subplot(2, 3, ctr+(lang-1)*3)];
        idx = surp_encoding.ls==ls & (uv_surp> 0 | uv_word> 0);
        x = uv_surp(idx);
        y = uv_word(idx);
    
        scatter(x(x>0 & y>0), y(x>0 & y>0), 15, ...
            cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.8); hold on;
        scatter(x(x<0 | y<0), y(x<0 | y<0), 10, ...
            [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.6); hold on;
    
        ylim([-0.001 0.028])
        xlim([-0.001 0.01])

        % correlation
        [r, p] = corr(x(x>0 & y>0), y(x>0 & y>0), 'Type', 'Spearman');
        title(['r=' num2str(r, 2) ', ' getSigStr(p, 2)], 'FontWeight', 'normal')

        % counts in each quadrant
        cnts = [sum(x<0 & y>0), sum(x>0 & y>0); sum(x<0 & y<0), sum(x>0 & y<0)];
        axes('Position',[.2+(ctr-1)*0.3 .8 .1 .1])
        imagesc(cnts)
        colormap(gray(8))
        caxis([0 50]);
    end
    ctr = ctr+1;
end
% linkaxes()


for ctr=1:6
    subplot(2, 3, ctr);

    set(ax(ctr), 'FontSize', 13)

    xline(0); 
    yline(0);

    % identity line
    h=refline(1, 0);
    h.LineWidth = 1.5;    
    h.Color = 'k';

    if ctr == 1
        xlabel(['English ' label '  \DeltaR^2']);
        ylabel('English word \DeltaR^2');
    elseif ctr == 4
        xlabel(['Spanish ' label ' \DeltaR^2']);
        ylabel('Spanish word \DeltaR^2');
    end

    chi=get(gca, 'Children');


    %Reverse the stacking order so that the patch overlays the line
    set(gca, 'Children',flipud(chi));
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;

%% vis: surprisal/word in English vs. Spanish
cols = getColorsCrossComp(1);

figure;
eng_surp = surp_encoding.eng_uv_surprisal;
sp_surp  = surp_encoding.sp_uv_surprisal;
thresh = 0;

ctr = 1;
for ls = [1, 2, 4]
    subplot(2, 3, ctr)
    idx = ismember(surp_encoding.ls, ls)  & ...
        (eng_surp > thresh | sp_surp > thresh);
    x = sp_surp(idx);
    y = eng_surp(idx);
    scatter3(x, y, find(idx), 25, ...
        cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.6);
    cnts = [sum(x<0 & y>0), sum(x>0 & y>0); sum(x<0 & y<0), sum(x>0 & y<0)];
    
    % formatting
    xline(0); 
    yline(0);
    view(2);
    grid off;

    if ctr==3
        xlabel('Spanish surprisal \DeltaR^2');
        ylabel('English surprisal \DeltaR^2');
    else
        xticks([]);
        yticks([]);
    end
    
    set(gca, 'FontSize', 13);
    ylim([-0.02 0.02]);
    xlim([-0.02 0.02]);
    h=refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';

    subplot(2, 3, ctr+3);
    imagesc(cnts);
    colormap(gray(5));
    caxis([0 60])
    for x=1:2
        for y=1:2
            if cnts(x, y)>30
                text(y, x, num2str(cnts(x, y)), 'FontSize', 13);         
            else
                text(y, x, num2str(cnts(x, y)), 'FontSize', 13, 'Color', 'w');         
            end
        end
    end
    yticks([]);
    xticks([]);

    ctr = ctr+1;
end

figure;
eng_word = surp_encoding.eng_uv_word;
sp_word  = surp_encoding.sp_uv_word;

ctr = 1;
for ls = [1, 2, 4]
    subplot(2, 3, ctr)
    idx = ismember(surp_encoding.ls, ls) & ...
        (eng_word > thresh | sp_word > thresh);
    x = sp_word(idx);
    y = eng_word(idx);
    cnts = [sum(x<0 & y>0), sum(x>0 & y>0); sum(x<0 & y<0), sum(x>0 & y<0)];

    scatter3(x, y, find(idx), 25, ...
        cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.6);
    % formatting
    xline(0); 
    yline(0);
    grid off;
    view(2)

    if ctr==3
        xlabel('Spanish word \DeltaR^2');
        ylabel('English word \DeltaR^2');
    else
        xticks([]);
        yticks([]);
    end
%     colormap([cols([1:2 4], :)]); % cols(4, :)
    set(gca, 'FontSize', 13);
    ylim([-0.01 0.03]);
    xlim([-0.01 0.038]);
    h=refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';
    ylim([-0.01 0.03]);   

    subplot(2, 3, ctr+3);
    imagesc(cnts);
    colormap(gray(5));
    caxis([0 160])
    for x=1:2
        for y=1:2
            if cnts(x, y)>60
                text(y, x, num2str(cnts(x, y)), 'FontSize', 13);         
            else
                text(y, x, num2str(cnts(x, y)), 'FontSize', 13, 'Color', 'w');         
            end
        end
    end
    yticks([]);
    xticks([]);

    ctr = ctr+1;

end

figure('Renderer','painters');
for ls = [1, 2]
    idx = ismember(surp_encoding.ls, ls) & ...
        (eng_word > thresh | sp_word > thresh);
    scatter3(sp_word(idx), eng_word(idx), find(idx), 45, ...
        cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.6); hold on;
end

% formatting
xline(0); 
yline(0);
view(2)
xlabel('Spanish word \DeltaR^2');
ylabel('English word \DeltaR^2');

colormap([cols([1:2 4], :)]); % cols(4, :)
set(gca, 'FontSize', 13);
ylim([-0.01 0.03]);
xlim([-0.01 0.038]);
h=refline(1, 0);
h.LineWidth = 1.5;
h.Color = 'k';
h.HandleVisibility = 'off';
ylim([-0.01 0.03]);
legend({'Spanish', 'English'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;

%% vis: example surprisal beta weights

elecs = containers.Map;
% english
elecs('EC163') = 8; % , 105, 106, 185
% elecs('EC105') = [166, 168, 106, 169];
elecs('EC183') = [55, 71, 72];
% spanish
elecs('EC172') = [105, 122];
elecs('EC100') = 118;

% elecs('EC172') = 122;
% elecs('EC214') = [69];
% elecs('EC100') = [5 21 69 118];
% 
% % english
% elecs('EC183') = [56 72]; % E1 and E2, ERP electrode
% % elecs('EC222') = [117 85]; % E1 and E2, ERP electrode
% 
% % bilingual
% % elecs('EC139') = [147 178 ]; % E1 and E2, ERP electrode
% 
% % for lab meeting: 
% % elecs('EC100') = 69; elecs('EC163') = 8; % spanish
% elecs('EC222') = 117; spanish

modelnames = {{'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpNoOnsBin'}, ...
     {'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpNoOnsBin'}};
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};

for s = elecs.keys()
    % selected electrodes
    SID = s{1};

    for l = 1:2
        corpusStrf{l} = loadMultModelStrf(SID, modelnames{l}, corpus{l}, ...
            datapath, 1, 'v5');
    end

    for el = elecs(SID)
%         figure;
%         subplot(1, 2, 1);
%         imagesc(corpusStrf{1}{1}.meanStrf(19:26, :, el))
% 
%         subplot(1, 2, 2);
%         imagesc(corpusStrf{2}{1}.meanStrf(20:27, :, el));
%         colormap(flipud(brewermap(50, 'Spectral')));       
        for ls = [1, 2]
            feats = details{ls}.features.names;
            surpfeat = 5+length(feats);
            strf = smoothdata(corpusStrf{ls}{1}.meanStrf, 2);
            [~, tp] = max(mean(abs(strf(surpfeat:surpfeat+7, :, el))));
            y = mean(strf(surpfeat:surpfeat+7, max(1, tp-5):min(tp+5, 61), el), 2);

            figure;
            subplot(1, 3, 1);
            imagesc(strf(2:2+length(feats), :, el));
            yticks(1:length(feats)+1);
            yticklabels([feats; {'peakRate'}]);
            b = balanced(8);
            colormap([b(1:4, :); [1 1 1; 1 1 1; 1 1 1]; b(5:8, :)]);
            caxis([-0.5 0.5]);

            cols = flipud(brewermap(8, 'Spectral'));

            subplot(1, 3, 2);
            [r, p] = corr((1:8)', y, 'type', 'Spearman');           
            scatter(1:8, y, 35, cols, 'filled');hold on;
            l=lsline();
            l.LineWidth = 1;
            text(1, 0.1, {['r=' num2str(r) ],[' p=' num2str(p, 2)]});
            xlim([0.5 8.5])

        
            subplot(1, 3, 3);   
            for i = 1:8
                x = -0.6:0.01:0;
                plot(x, flipud(smoothdata(strf(surpfeat+(i-1), :, el)')), ...
                    'Color', cols(i, :), 'LineWidth', 2); hold on;   
                xline(x(61-tp));
            end
            ylabel('beta weight');
            xlabel('time (s)');
            ylim([-0.25 0.25]);
            %ylim([-5 5]);
            set(gca, 'FontSize', 15);
            box off;


            sgtitle([SID ', el - ' num2str(el) ': ' corpus{ls}]);
        end
    end
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;


%% vis: summary of direction of surprisal beta weights

modelnames = {{'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpNoOnsBin'}, ...
     {'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpNoOnsBin'}};
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};

rho = nan(2, height(surp_encoding));
p = nan(2, height(surp_encoding));
for l = 1:2
    
    feats = details{l}.features.names;  
    surpfeat = 5+length(feats);            
    
    for s = unique(surp_encoding.SID)'
        SID = s{1};
        sididx = strcmp(surp_encoding.SID, SID);
    
        corpusStrf = loadMultModelStrf(SID, modelnames{l}, corpus{l}, ...
                datapath, 1, 'v5');
        strf = smoothdata(corpusStrf{1}.meanStrf, 2);
    
    
        for i = find(sididx)'
            el = surp_encoding.el(i);
            [~, tp] = max(mean(abs(strf(surpfeat:surpfeat+7, :, el))));
            y = mean(strf(surpfeat:surpfeat+7, max(1, tp-5):min(tp+5, 61), el), 2);
            if length(max(1, tp-5):min(tp+5, 61))>8
                [rho(l, i), p(l, i)] = corr((1:8)', y, 'type', 'Spearman');
    
                if p(l, i)<0.0001
                    figure;
                    subplot(1, 3, 1);
                    imagesc(zscore(strf(2:2+length(feats), :, el)));
                    yticks(1:length(feats));
                    yticklabels(feats);
                    colormap(balanced);
                    caxis([-10 10])
%                     colormap([ones(20, 1), linspace(1, 0, 20)', linspace(1, 0, 20)']);
    
                    cols = flipud(brewermap(8, 'Spectral'));
        
                    subplot(1, 3, 2);
                    scatter(1:8, y, 35, cols, 'filled');
                    lsline
                
                    subplot(1, 3, 3);   
                    for j = 1:7
                        x = -0.6:0.01:0;
                        plot(x, flipud(smoothdata(strf(surpfeat+(j-1), :, el)')), ...
                            'Color', cols(j, :), 'LineWidth', 2); hold on;   
                        xline(x(61-tp));
                    end
                    ylabel('beta weight');
                    xlabel('time (s)');
                    ylim([-0.25 0.25]);
                    set(gca, 'FontSize', 15);
                    box off;
    
                    sgtitle([SID ', el - ' num2str(el) ': ' corpus{l}]);
                end
            end
        end
    end
end

fields = {'eng_uv_surprisal', 'sp_uv_surprisal'};
cols = [1, 0, 0; 0, 0, 1];
figure; 
for l = 1:2
    subplot(1, 2, l);
    kshist(rho(l, surp_encoding.(fields{l})>0), 12, cols(l, :));
    xline(0, 'LineWidth', 1.75, 'Color', 'k', 'LineStyle', '--');
    xlim([-1 1]);
    yyaxis left
    ylim([0 0.21])
    

    set(gca, 'FontSize', 15);
end


%% vis: maximal phonetic feature encoding per surprisal electrode

modelnames = {{'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpNoOnsBin'}, ...
     {'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpNoOnsBin'}};
corpus = {'timit', 'dimex'};
fields = {'eng_uv_surpword', 'sp_uv_surpword'};
featval=nan(height(surp_encoding), 2);
featid=nan(height(surp_encoding), 2);
pr= nan(height(surp_encoding), 2);

for l = 1:2
    
    feats = details{l}.features.names;  
%     surpfeat = 5+length(feats);            
    
    for s = unique(surp_encoding.SID)'
        SID = s{1};
        sididx = strcmp(surp_encoding.SID, SID);

    
        corpusStrf = loadMultModelStrf(SID, modelnames{l}, corpus{l}, ...
                datapath, 1, 'v5');
        strf = smoothdata(corpusStrf{1}.meanStrf, 2);
    
        for i = find(sididx)'
            el = surp_encoding.el(i);
            [featval(i, :), featid(i, :)] = ...
                maxk(sum(abs(strf(2:1+length(feats), :, el)), 2), 2);
            pr(i, l) =  mean(strf(2+length(feats), :, el), 2);
        end
    end

    figure; 
    idx = abs(diff(featval'))>prctile(abs(diff(featval')), 50) ...
        & surp_encoding.(fields{l})'>0; 
    histogram(featid(idx, 1), 'facecolor', 'k', 'facealpha', ...
        0.5, 'edgecolor', 'none')
    xticks(1:length(feats))
    xticklabels(feats)
end


%% vis: MNI Brain for word/surprisal encoding 

% initialize design electrode structure
desel=struct();
desel.conds = 1:6;
ls = [1, 2, 4];
% size indicates the unique variance explained by matching predictor
desel.sz = [5; 35*ones(length(desel.conds), 1)]; %1:20:200; %ones(1, 10)*0.00000001; %1

% split up peak rate and phonetic features again for MNI plotting
desel.labels = [];
field = 'sp_uv_surpword';
posidx = surp_encoding.(field)>-1;

% [~, binedges] = discretize(surp_encoding.(field)(posidx), length(desel.conds));
[~, binedges] = discretize(surp_encoding.(field)(surp_encoding.(field)>0 & surp_encoding.ls==ls), ...
    length(desel.conds)-1);
binedges = [-1, binedges];

for s=unique(surp_encoding.SID)'
    SID = s{1};
    idx = strcmp(surp_encoding.SID, SID) & posidx;
    desel.(SID).elid = surp_encoding.el(idx);
    desel.(SID).condition = discretize(surp_encoding.(field)(idx), ...
        binedges);
end

% legend for size
x = [1 3 ]; %5 7 10
scatter(x, ones(length(x), 1), desel.sz(x), 'MarkerEdgeColor','k');
text([0, 9], [0.85 0.85], {'min', 'max'}, 'FontSize',14)
% text(x-1, ones(length(x), 1)*0.8, ...
%     split(num2str(binedges(x)+diff(binedges(x:x+1))/2)), 'FontSize', 13);
box off;
axis off;
xlim([-5 15]);
text(4, 0.7, 'Word \Delta R^2', 'FontSize',14);

% plot histogram of count across y-axis (anterior to posterior)
colors = getColors(1);
desel.cols = [linspace(1, colors(ls, 1), 10); ...
    linspace(1, colors(ls, 2), 10); linspace(1, colors(ls, 3), 10)]';
desel.cols = [1 1 1; desel.cols(3:end, :)];
% desel.cols = [0 0 0 ; flipud(greens(10))];
tmp = flipud(internet(5));
desel.cols = [1 1 1 ; tmp];

% colorbar
figure;
colormap(desel.cols);
colorbar

% desel.cols = [1 1 1; 212/256, 228/256 188/256; ...
%     54/256, 85/256, 183/256; 64/256 55/256 110/256];
% lsid = find(ismember(surp_encoding.ls, [1, 2, 4]));
lsid = find(ismember(surp_encoding.ls, ls));
[mni_lh] = plotMNIElec(unique(surp_encoding.SID(lsid)), desel, 'lh', 0);
[mni_rh] = plotMNIElec(unique(surp_encoding.SID(lsid)), desel, 'rh', 0);

sids = {'EC260', 'EC266', 'EC100', 'EC105', 'EC172'}; %names(startsWith(names, 'EC'));
desel.EC260.selid = [221 236 205 229 228];
desel.EC266.selid = [178 150];

[native_plot] = plotNativeElec(sids, desel);

% [mni_lh] = plotNativeElec(unique(surp_encoding.SID(lsid)), desel);

% for ls = [1, 2, 4]
%     desel.cols = [repmat([cols(ls, :)], 10, 1)];
%     lsid = find(surp_encoding.ls == ls);
%     mni_rh = [];
%     if ismember(ls, [1, 2, 4]) 
%         [mni_lh] = plotMNIElec(unique(surp_encoding.SID(lsid)), desel, 'lh', 0);
%         [mni_rh] = plotMNIElec(unique(surp_encoding.SID(lsid)), desel, 'rh', 0);
%     else
%          % no RH HS subjects
%         [mni_lh] = plotNativeElec({'HS11'}, desel);
%     end
% end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: full model fit, subject group comparison

figure; 
% full model R^2
% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
% fieldname = {'English', 'Spanish'};
thresh = 0.0001;

% label = {'Spanish', 'English', '', 'Bilingual'};
% ctr = 1;
% for ls = [1:2 4]
%     subplot(1, 3, ctr)
%     rsq = cell2mat(arrayfun(@(x) ...
%         surp_encoding.(field{x})(surp_encoding.ls==ls), ...
%         1:2, 'UniformOutput',false));
% 
%     scatter(rsq(:, 2), rsq(:, 1), 35, [0.5 0.5 0.5]);
%     set(gca, 'FontSize', 15)
%     box off;
%     ylabel('full English model R^2'); 
%     xlabel('full Spanish model R^2');
%     title([label{ls} ' Subjects']);
%     h = refline(1, 0);
%     h.Color = 'k';
%     h.LineWidth = 0.8;
%     ctr = ctr+1;
% end

colors = getColors(1);
% fig = figure();
% p(1) = figure();
% p(2) = figure();
% 
% label = 'word boundary';
% fields = {'eng_uv_surpword', 'sp_uv_surpword'};

label = 'peakrate ';
fields = {'eng_uv_peakrate', 'sp_uv_peakrate'}; 

% label = 'phonetic feature ';
% fields = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 

%fields = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 

% fields = {'eng_uv_word', 'sp_uv_word'};
% label = 'word';

figure;
ctr=1;
for ls = [1 2 4]
    ax(ctr)=subplot(1, 3, ctr);
    scatter(surp_encoding.(fields{1})(surp_encoding.ls==ls), ...
        surp_encoding.(fields{2})(surp_encoding.ls==ls), 15, ...
        colors(ls, :), 'filled');
    r = nancorr(surp_encoding.(fields{1})(surp_encoding.ls==ls), ...
        surp_encoding.(fields{2})(surp_encoding.ls==ls));

    maxlim = prctile([surp_encoding.(fields{1})(surp_encoding.ls==ls); ...
        surp_encoding.(fields{2})(surp_encoding.ls==ls)], 99);
    xlim([-0.01 maxlim+0.005]);
    ylim([-0.01 maxlim+0.005]);
    h=refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    xline(0, 'Color', 'k');
    yline(0, 'Color', 'k');

    xlabel(['English ' label]);
    ylabel(['Spanish ' label]);
    title(['r: ' num2str(r) ]);
    ctr=ctr+1;
end

titles = {'English', 'Spanish'};

for f = 1:2
    field = fields{f};    
    ctr = 1;
    figure;
    subplot(3, 2, [1 3 5])
    for ls = [1 2 4]
        lsidx = find(surp_encoding.ls==ls & surp_encoding.(field)>thresh);
    %         z = cellfun(@(x) str2double(x(3:end)), surp_encoding.SID);
    %         z = surp_encoding.el;
        
        scatter(ctr-0.2 + rand(length(lsidx), 1)*0.05, surp_encoding.(field)(lsidx), ...
            10, colors(ls, :), 'filled', 'LineWidth', 0.1, 'MarkerFaceAlpha', 0.6); hold on;
        
        nsidx = find(surp_encoding.ls==ls & surp_encoding.(field)<thresh);
        scatter(ctr-0.2 + rand(length(nsidx), 1)*0.05, surp_encoding.(field)(nsidx), ...
            10, [0.5 0.5 0.5], 'filled', 'LineWidth', 0.1, 'MarkerFaceAlpha', 0.6); hold on;
        boxplot(surp_encoding.(field)(lsidx), 'Colors',colors(ls, :), 'Positions', ...
            ctr, 'widths', 0.45, 'Symbol', '', 'BoxStyle', 'filled'); hold on;

        % clip the violin outliers
        %         lsidx = find(surp_encoding.ls==ls & abs(surp_encoding.(field))<0.015);
        %         violin(surp_encoding.(field)(lsidx), 'facecolor',colors(ls, :), 'x', ...
        %             ctr, 'medc', '', 'bw', 0.001); hold on;
        %         v = violin(surp_encoding.(field)(lsidx), 'facecolor', colors(ls, :), ...
        %             'medc', [], 'bw', 0.0015);
        %         v.XData = v.XData+ctr-1;
        ctr = ctr + 1;
    end
    
    % formatting
    legend('off')
    set(gca, 'FontSize', 15); % , , 'Yscale', 'log'
    xlabel('Subject group');
    ylabel([label ' \DeltaR^2 ']);
    xlim([0.25 3.75]);
    ylim([-0.015 0.045]);
    xticks(1:3);

    xticklabels({'Spanish', 'English', 'Bilingual'});
    box off;
    
    h=yline(0, '--k');
    h.LineWidth = 0.8;
    yticks(-0.01:0.01:0.04)
    title(titles{f});

    idx = surp_encoding.(field)>-1;
    inc = 0.005;
    
    disp(['----------------- ' titles{f} ' ------------------------'])
%     [~, P, ~] = kstest2(surp_encoding.(field)(surp_encoding.ls==1 & idx), ...
%         surp_encoding.(field)(surp_encoding.ls==2 & idx), 'Tail', 'unequal');
    [P, ~, ~] = ranksum(surp_encoding.(field)(surp_encoding.ls==1 & idx), ...
        surp_encoding.(field)(surp_encoding.ls==2 & idx)); % ), 'Tail', 'larger'
    disp(['rank sum p-val [MONO]: ' num2str(P)]);
    [str, sig] = getSigStr(P, 2); 
    if sig, line([1 2], [0.03 0.03], 'Color', 'k', 'LineWidth', 1.5); end
    text( 1.5, 0.03+0.005, str, 'FontSize', 13);
    
%     [~, P, ~] = kstest2(surp_encoding.(field)(surp_encoding.ls==2 & idx), ...
%         surp_encoding.(field)(surp_encoding.ls==4 & idx), 'Tail', 'unequal');
    [P, ~, ~] = ranksum(surp_encoding.(field)(surp_encoding.ls==2 & idx), ...
        surp_encoding.(field)(surp_encoding.ls==4 & idx));
    disp(['rank sum p-val [ENG to BIL]: ' num2str(P)]);
    [str, sig] = getSigStr(P, 2); 
    if sig, line([2 3], [0.03+inc 0.03+inc], 'Color', 'k', 'LineWidth', 1.5); end
    text( 2.5, 0.03+inc+0.005, str, 'FontSize', 13);
    
%     [~, P, ~] = kstest2(surp_encoding.(field)(surp_encoding.ls==1 & idx), ...
%         surp_encoding.(field)(surp_encoding.ls==4 & idx), 'Tail', 'unequal');
    [P, ~] = ranksum(surp_encoding.(field)(surp_encoding.ls==1 & idx), ...
        surp_encoding.(field)(surp_encoding.ls==4 & idx));
    disp(['rank sum p-val [SP to BIL]: ' num2str(P)]);
    [str, sig] = getSigStr(P, 2); 
    if sig, line([1 3], [0.03+inc*2 0.03+inc*2], 'Color', 'k', 'LineWidth', 1.5); end
    text( 2, 0.03+inc*2+0.005, str, 'FontSize', 13);

    ctr=1;
    for ls = [1 2 4]
        lsidx = find(surp_encoding.ls==ls & surp_encoding.(field)>thresh);

        subplot(3, 2, ctr*2);
        p = pie(sum([surp_encoding.(field)(surp_encoding.ls==ls)>0, ...
            surp_encoding.(field)(surp_encoding.ls==ls)<0]), ...
            [1, 1]); hold on;
        p(3).FaceColor = [0.7 0.7 0.7];
        p(3).EdgeColor = 'none';
        p(1).FaceColor = colors(ls, :);
        p(1).EdgeColor = 'none';

        ctr = ctr + 1;
    end    
end
% add linear mixed effect model

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;

%% vis: word onset compared to aud unique variance
% look into EC222 ec 123
cols = getColorsCrossComp(1);
corpus = {'Listening to English', 'Listening to Spanish'};
figure;
for c = 1:2
    uvwrd = word_encoding.uv_word(:, c);
    base = word_encoding.base_rsq(:, c);

    subplot(1, 2, c)
    colors = getColorsCrossComp(1);
    for ls = 1:3
        idx = word_encoding.ls==ls;
        scatter(ls-0.1+ones(1, sum(idx)).*rand(1, sum(idx))*0.2, ...
            uvwrd(idx), 45, colors(ls, :),'filled', 'MarkerFaceAlpha', ...
            0.6); hold on;
        boxplot(uvwrd(idx), 'Position', ls, 'Width', 0.5, 'Color', 'k');
        line([ls-0.15 ls+0.15], [median(uvwrd(idx), 'omitnan') ...
            median(uvwrd(idx), 'omitnan')], ...
            'LineWidth', 1.5, 'Color', 'k')
    end
    set(gca,  'FontSize', 13); % 'yscale', 'log',
    xticks(1:3);
    xlim([0.5 2.5]);
    ylim([-0.0001 0.0016]);
    xticklabels({'Spanish', 'English', 'Mandarin'});
    ylabel('Word \DeltaR^2')
    title(corpus{c});
    yline(0);
end

figure;
% sid = cellfun(@(x) str2double(x(3:end)), ...
%     word_encoding.SID);
idx = word_encoding.ls<3;
x = word_encoding.uv_word(:, 2); % spanish onset UV
y = word_encoding.uv_word(:, 1); % english onset UV
scatter3(x(idx), y(idx), word_encoding.el(idx), 42, ...
    word_encoding.ls(idx), 'filled', 'MarkerFaceAlpha', 0.6);

%formatting
xlabel('Spanish word \DeltaR^2');
ylabel('English word \DeltaR^2');
colormap(cols(1:2, :));
ylim([-0.0006 0.0015]);
xlim([-0.0006 0.0015]);
yline(0);
xline(0);
refline(1, 0);
view(2);
set(gca, 'FontSize', 13);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;

%% vis: same as above, size as peakRate variance
% look into EC222 ec 123
cols = getColorsCrossComp(1);
corpus = {'Listening to English', 'Listening to Spanish'};
figure;
for c = 1:2
    uvwrd = word_encoding.uv_word(:, c);
    base = word_encoding.base_rsq(:, c);

%     subplot(1, 3, 3);
    subplot(1, 2, c);
    for ls = 1:3
        idx = word_encoding.ls==ls;
        scatter(ls-0.1+ones(1, sum(idx)).*rand(1, sum(idx))*0.2, ...
            uvwrd(idx), 45, [0.5 0.5 0.5],'filled', 'MarkerFaceAlpha', ...
            0.3); hold on;
        boxplot(uvwrd(idx), 'Colors','k', 'Positions', ls, 'widths', 0.45, ...
            'OutlierSize', 2);
    end

    ylabel('Word \DeltaR^2');
    colormap(cols(1:2, :));
    xlim([0.5 2.5]);
    ylim([-0.0006 0.002]);
    yline(0);
    xline(0);
    set(gca, 'FontSize', 13); % 'yscale', 'log'

    xticks(1:3);
    xticklabels({'Spanish', 'English', 'Mandarin'});
    title(corpus{c});
    yline(0);
    box off;
end

figure;
z = cellfun(@(x) str2double(x(3:end)), ...
    word_encoding.SID);
y = word_encoding.uv_word(:, 1);
x = word_encoding.uv_word(:, 2);
sz = abs(word_encoding.uv_peakRate(:, 1)*20^5);
idx = word_encoding.ls<3;
scatter3(x(idx), y(idx), word_encoding.el(idx), sz(idx), ...
    word_encoding.ls(idx), 'filled', 'MarkerFaceAlpha', 0.6); hold on;

%formatting
colormap(cols(1:2, :));
view(2);
yline(0);
xline(0);
h = refline(1, 0);
h.Color = 'k';

ylim([-0.0001 0.003]);
xlim([-0.0001 0.003]);
xlabel('Spanish Word \DeltaR^2');
ylabel('English Word \DeltaR^2');
title(corpus{c});

set(gca, 'FontSize', 13); % , 'yscale', 'log'
box off;

% selected electrodes
elecs = containers.Map;
elecs('EC100') = 21; 
elecs('EC172') = 122; 
elecs('EC222') = [101, 90];
elecs('EC235') = [249, 250];
ctr = 1;
for k = elecs.keys()
    % selected electrodes
    key = k{1};
    for e = elecs(key)
        idx = strcmp(word_encoding.SID, key) & ...
            word_encoding.el == e;
        text(x(idx), y(idx), ['E' num2str(ctr)], 'FontSize', 13);
        disp(['E' num2str(ctr) ': ' key ', el = ' num2str(e)]);
        ctr = ctr+1;
    end
end

clearvars -except *all subj *vow* *cons* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% vis: TRF weights for surprisal and word UV

row = 23; % 11, 43
modelnames={'onset_phnfeatonset_maxDtL', 'onset_aud'};   % 'onset_phnfeatonset_maxDtL'{'onset_maxDtL_aud'};   % 
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};


el = word_encoding.el(row);
SID = word_encoding.SID{row};

figure
for l = 1:2
    corpusStrf = loadMultModelStrf(SID, modelnames, corpus{l}, ...
        datapath, 1);
    featnames = details{l}.features.names;
    
    subplot(2, 2, l);
    imagesc(squeeze(corpusStrf{1}.meanStrf(2:end-1, :, el)));
    if l == 1
        yticks(1:2+length(featnames));
        yticklabels(featnames') % {'onset'}, , {'peakRate'}
    else
        yticks([])
    end

    subplot(2, 2, l+2);
    imagesc(squeeze(corpusStrf{2}.meanStrf(:, :, el)));
    if l == 2
        yticks([])
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *cons*;

%% vis: showing frequency selectivity across word encoding types

modelnames={'onset_phnfeatonset_maxDtL', 'onset_aud'};   % 'onset_phnfeatonset_maxDtL'{'onset_maxDtL_aud'};   % 
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};

aud_weights = nan(height(word_encoding), 80); % max tp frequency selectivity
feat_weights = nan(height(word_encoding), 14);

word_rsq = [];
rsq_idx = [];
for s = unique(word_encoding.SID)'
    SID = s{1};
    corpusStrf = loadMultModelStrf(SID, modelnames, corpus{1}, ...
        datapath, 1);

    if ~isempty(corpusStrf{2})

        els = strcmp(word_encoding.SID, SID);
        word_rsq = [word_rsq corpusStrf{1}.meanTestR(word_encoding.el(els)).^2];
        rsq_idx = [rsq_idx find(els)'];

        for row = find(els)'
            el = word_encoding.el(row);
            featnames = details{1}.features.names;       
            [~, maxtp] = max(mean(corpusStrf{2}.meanStrf(2:end, :, el), [1, 3]));
            wind = 3;
            range = max(1, maxtp-wind):min(60, maxtp+wind);
            % smoothed mean frequency selectivity at maximum tp
            y = smoothdata(mean(corpusStrf{2}.meanStrf(2:end, ...
                range+1, el), [2, 3]), 'movmean', 'SmoothingFactor',0.08);
            y = rescale(y,'InputMin', min(y), 'InputMax', max(y));  

            [~, maxtp] = max(mean(corpusStrf{1}.meanStrf(2:end-1, :, el), [1, 3]));
            range_feat = max(1, maxtp-wind):min(60, maxtp+wind);
            y_feat = smoothdata(squeeze(corpusStrf{1}.meanStrf(2:end-1, :, el)), ...
                2, 'movmean', 'SmoothingFactor',0.08);
            [~, max_feat(row)] = max(mean(y_feat(:, range_feat), 2));

            debug = 0;
            if debug
                figure;
                subplot(2, 2, 1);
                imagesc(squeeze(corpusStrf{2}.meanStrf(:, :, el))); hold on;
                rectangle('Position', [range(1), 1, ...
                    length(range), 80], 'LineWidth', 2);

                subplot(2, 2, 2);
                plot(y);
                title(corpusStrf{2}.meanTestR(el)^2)

                subplot(2, 2, 3)
                imagesc(y_feat); hold on;
                rectangle('Position', [range_feat(1), 0.25, ...
                    length(range_feat), 14.5], 'LineWidth', 2);

                subplot(2, 2, 4)
                plot(y_feat'); hold on;
                rectangle('Position', [range_feat(1), min(y_feat(:))-0.05, ...
                    length(range_feat), abs(min(y_feat(:)))+max(y_feat(:))], 'LineWidth', 2);
            end

            if corpusStrf{2}.meanTestR(el)^2>0.0001
                aud_weights(row, :)  = y;
            end
        end
    end
end

[g, ~, ~, labels] = crosstab(max_feat);

x = cellfun(@(x) str2double(x), labels);
b = bar(x, g);
b(1).FaceColor = [0.5 0.5 0.5];

xticks(1:14);
xlim([0.4 14.8]);
xticklabels(featnames);

ylabel('Electrode count');
xlabel('Primary Phonetic Feature');
set(gca, 'FontSize', 15);
box off;

% r-squared comparison
figure;
scatter(word_encoding.base_rsq(rsq_idx, 1), word_rsq, 55, 'k', ...
     'filled', 'MarkerFaceAlpha', 0.5); 

ylim([0 0.55]);
xlim([0 0.55]);

h=refline(1, 0);
h.LineWidth = 1.5;
h.LineStyle = '-';
h.Color = 'k';
ylabel('phonetic feature R^2');
xlabel('mel spectrogram R^2');
set(gca, 'FontSize', 15);

% beta weight comparison (what does the average frequency selectivity look like?)
figure;   
[~, tpidx] = max(aud_weights, [], 2);
[~, sortidx] = sort(tpidx);
imagesc(aud_weights(sortidx, :)); %.*rhos(v_minus(sortidx)));
ylabel('Electrode count');
xlabel('Frequency Bin');
set(gca, 'FontSize', 15);
box off;

cm = brewermap(50, 'Greys');
colormap(cm(1:end-8, :));

cbh = colorbar();
ylabel(cbh, 'STRF weight (% max)');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *cons*;


%% vis: full model fit, subject group comparison for CROSS CASE

figure;
fields = {'engcross_uv_surprisal'};
label = 'phoneme surprisal';
% fields = {'eng_uv_word', 'sp_uv_word'};
% label = 'word';

titles = {'Listening to: English', 'Listening to: Spanish'};

field = fields{1};
ctr = 1;
for ls = [1:2 4]
    lsidx = find(surp_encoding.ls==ls & surp_encoding.(field)>0);
%         z = cellfun(@(x) str2double(x(3:end)), surp_encoding.SID);
%         z = surp_encoding.el;

    scatter(ctr-0.1 + rand(length(lsidx), 1)*0.2, surp_encoding.(field)(lsidx), ...
        45, [0.7 0.7 0.7], 'LineWidth', 1.75); hold on;
    boxplot(surp_encoding.(field)(lsidx), 'Colors','k', 'Positions', ...
        ctr, 'widths', 0.45, 'Symbol','');
    ctr = ctr + 1;
end
view(2)
xlim([0.25 3.75]);
ylim([-0.002 0.015]);
xticks(1:3);
xticklabels({'Spanish', 'English', 'Bilingual'});
set(gca, 'FontSize', 15); % , , 'Yscale', 'log'
xlabel('Subject group');

box off;
ylabel([label ' \DeltaR^2 ']);
h=yline(0, '--k');
h.LineWidth = 0.8;
yticks(0:0.01:0.04)
title(titles{1})

idx = surp_encoding.(field)>0;

disp(['----------------- ' titles{1} ' ------------------------'])
[P,~, ~] = ranksum(surp_encoding.(field)(surp_encoding.ls==1 & idx), ...
    surp_encoding.(field)(surp_encoding.ls==2 & idx));
disp(['rank sum p-val [MONO]: ' num2str(P)]);
if P<0.05 
    line([1 2], [0.03 0.03], 'Color', 'k', 'LineWidth', 1.5);
    text( 1.5, 0.0305, getSigStr(P), 'FontSize', 13);
end

[P,~, ~] = ranksum(surp_encoding.(field)(surp_encoding.ls==2 & idx), ...
    surp_encoding.(field)(surp_encoding.ls==4 & idx));
disp(['rank sum p-val [ENG to BIL]: ' num2str(P)]);
if P<0.05 
    line([2 3], [0.032 0.032], 'Color', 'k', 'LineWidth', 1.5);
    text( 2.5, 0.0325, getSigStr(P), 'FontSize', 13);
end

[P,~, ~] = ranksum(surp_encoding.(field)(surp_encoding.ls==1 & idx), ...
    surp_encoding.(field)(surp_encoding.ls==4 & idx));
disp(['rank sum p-val [SP to BIL]: ' num2str(P)]);
if P<0.05 
    line([1 3], [0.034 0.034], 'Color', 'k', 'LineWidth', 1.5);
    text( 2, 0.0345, getSigStr(P), 'FontSize', 13);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd;


%% vis: uv surprisal versus surprisal average

cols = getColorsCrossComp(1);

figure;
maxls = 5;
ax = subplot(2, 3, 1);
ctr = 1;
for ls = [1, 2, 4]
    subplot(2, 3, ctr)

    uv_surp = surp_encoding.eng_uv_surprisal;
    uv_word = surp_encoding.eng_uv_surpaverage;

    idx = surp_encoding.ls==ls & (uv_surp > 0 | uv_word> 0);
    scatter(uv_surp(idx), uv_word(idx), 25, cols(ls, :), ...
        'filled', 'MarkerFaceAlpha', 0.6);

    % formatting
    xline(0); 
    yline(0);
    if ctr == 1
        xlabel('English surprisal \DeltaR^2');
        ylabel('English surprisal average \DeltaR^2');
    end
    colormap(cols(:, :));
    set(gca, 'FontSize', 13);
    xlim([-0.02 0.045]);
    ylim([-0.01 0.045]);
    h=refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';
    ylim([-0.01 0.045]);

    ax=subplot(2, 3, ctr+3);
    uv_surp = surp_encoding.sp_uv_surprisal;
    uv_word = surp_encoding.sp_uv_surpaverage;
    idx = surp_encoding.ls==ls & (uv_surp > 0 | uv_word > 0);
    scatter(uv_surp(idx), uv_word(idx), 45, ...
        cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.6);
    
    % formatting
    xline(0); 
    yline(0);
    if ctr == 1
        xlabel('Spanish surprisal \DeltaR^2');
        ylabel('Spanish surprisal average \DeltaR^2');
    end
    colormap(cols(:, :));
    set(gca, 'FontSize', 13);
    xlim([-0.02 0.045]);
    ylim([-0.01 0.045]);
    h=refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';
    ylim([-0.01 0.045]);

    ctr = ctr+1;

end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;


%% --------------------------- SCRATCH ------------------------------------
% scratch --- as struct
corpus_details = timit_details;
befaft = [1 1]; % information around boundary

% add surprise to the out struct
[corpus_details.sentdet]=makeSurprisal(corpus_details.sentdet, 8);

% find sentences where surprise is defined
% concatenate all sentences where surprise is defined
surpSent = zeros(1, length(corpus_details.sentdet));
allsent = struct();
allsent.aud=[];
allsent.surprise = [];
allsent.env = [];
allsent.wrdOns = [];
allsent.wrdLn = [];
allsent.precwrdLn = [];
allsent.precwrdSyl = [];
allsent.sentOns = [];
allsent.phn = [];
allsent.stress = [];
allsent.sentid = {};

for sent = 1:length(corpus_details.sentdet)

    % find all syllable onsets

    sent_info = corpus_details.sentdet(sent);
    dataf = sent_info.dataf;
    
    % find where surprisal is defined
    if  ~isempty(sent_info.befaft) && ~strcmp(sent_info.name, 's08204')
        sentp = (sent_info.befaft(1)*dataf):...
                size(sent_info.onsOff, 2)-sent_info.befaft(2)*dataf;
        wrdons = find(sent_info.wordOns(:, sentp)>0);
        % remove silence before and after sentence
    
    %     if any(sent_info.engSurp~=0)
    %         surpSent(sent) = 1;        
    %         allsent.surprise = cat(2, allsent.surprise, sent_info.engSurp(:, sentp));         
    %     end
    
        allsent.aud = cat(2, allsent.aud, sent_info.aud(:, sentp));
        allsent.env = cat(2, allsent.env, sent_info.loudness(:, sentp));
        allsent.sentOns = cat(2, allsent.sentOns, sent_info.onsOff(1, sentp));
        allsent.wrdOns = cat(2, allsent.wrdOns, sent_info.wordOns(:, sentp));

         % adding feature with current word length
        wrdln = zeros(1, length(sentp));
        wrdln(1, wrdons) = cellfun(@(x) length(x), sent_info.wordList);
        allsent.wrdLn = cat(2, allsent.wrdLn, wrdln);

        % adding feature with preceding word length
        pwrdln = zeros(1, length(sentp)); 
        pwrdln(1, wrdons(2:end)) = cellfun(@(x) length(x), ...
            sent_info.wordList(1:end-1));
        allsent.precwrdLn = cat(2, allsent.precwrdLn, pwrdln);

        tmp = arrayfun(@(x) find(sent_info.phnmatonset(:, x)), sentp, ...
            'UniformOutput', false);
        phns = -1*ones(1, length(sentp));
        % check that exactly one phoneme is occurring at time point ...
        % (ignoring cases where two phonemes present)
        idx = cellfun(@(x) length(x)==1, tmp);
        phns(idx) = [tmp{idx}];
        allsent.phn = cat(2, allsent.phn, phns);
        clear tmp phns
    
        % syllOns
        allsent.stress = cat(2, allsent.stress, ...
            sent_info.syltype(sentp));
        allsent.sentid = [allsent.sentid sent_info.name];

        debug = 0;
        if debug
            imagesc([sent_info.wordOns(:, sentp); ...
                sent_info.syltype(sentp)>0; ...
                wrdln; pwrdln])
        end
    end
end
surpSent = find(surpSent);

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons*;

% exploring what each predictor looks like at word onsets

field = 'env';
dataf = 100;
idx = allsent.wrdOns & ~allsent.sentOns & allsent.precwrdLn>0;
fieldlen = length(allsent.(field));
envWordOns = arrayfun(@(x) allsent.(field)(:, max(x-befaft(1)*dataf, 1):...
    min(x+befaft(2)*dataf, fieldlen)), find(idx), 'UniformOutput',false);
% find all word onsets segments with the same length
envWordOns = envWordOns(cellfun(@(x) size(x, 2) == sum(befaft*dataf)+1,...
    envWordOns));
envWordOns = cell2mat(envWordOns');
[~, sortidx] = sort(allsent.precwrdLn(idx));

% at syllOns
envSyllOns = arrayfun(@(x) allsent.(field)(:, max(x-befaft(1)*dataf, 1):...
    min(x+befaft(2)*dataf, fieldlen)), ...
    find(allsent.syllOns & ~allsent.wrdOns & ~allsent.sentOns), 'UniformOutput',false);
% find all word onsets segments with the same length
envSyllOns = envSyllOns(cellfun(@(x) size(x, 2) == sum(befaft*dataf)+1,...
    envSyllOns));
envSyllOns = cell2mat(envSyllOns');

figure;
subplot(2, 4, [1 2]);
y_wrd = normalize(envWordOns, 2, 'zscore')';
y_wrd = y_wrd(:, sortidx(sortidx<length(y_wrd)));
y_syll = normalize(envSyllOns, 2, 'zscore')';
leg = {'word', 'onset'};
if ~isempty(y_syll)
    shadedErrorBar(-befaft(1):0.01:befaft(2), mean(y_syll, 2), ...
        nansem(y_syll, 2), 'lineprops', {'r'});
    leg = {'syllable', 'word', 'onset'};
end
shadedErrorBar(-befaft(1):0.01:befaft(2), mean(y_wrd, 2), ...
    nansem(y_wrd, 2), 'lineprops', {'b'});
xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
legend(leg);
ylabel('norm env (z-score)'); xlabel('time from ons (s)')
set(gca, 'FontSize', 13);
ylim([-0.6 0.4]);

subplot(2, 4, 3);
if ~isempty(y_syll)
    imagesc('Xdata', -befaft(1):0.01:befaft(2), 'Cdata', y_syll');
    ylim([0 size(y_syll, 2)]);
    xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
    title('Syllable Onset')
end
subplot(2, 4, 4);
imagesc('Xdata', -befaft(1):0.01:befaft(2), 'Cdata', y_wrd');
ylim([0 size(y_wrd, 2)]);
xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
title('Word Onset')
xlabel('Time from onset (s)');
colormap(flipud(brewermap(45, 'RdBu')));
ylabel('trials'); xlabel('time from ons (s)')
set(gca, 'FontSize', 13);
% sliding window binomial logistic regression for word onset

subplot(2, 4, [5 6 7]);
C = categorical(allsent.phn(allsent.wrdOns==1 & allsent.phn~=-1), ...
    1:length(corpus_details.phnnames), corpus_details.phnnames);
histogram(C);
%xlim([categorical({'G'}), categorical({'x'})]);
ylabel('count'); xlabel('phoneme at onset')
set(gca, 'FontSize', 13);

subplot(2, 4, 8);
matidx = [find(strcmp(corpus_details.features.names, {'sonorant'})) ...
    find(strcmp(corpus_details.features.names, {'obstruent'}))];
obson = corpus_details.features.mat(matidx, allsent.phn(allsent.wrdOns==1 ...
    & allsent.phn~=-1));
C = categorical(obson(1, :) + 2.*obson (2, :), 0:2, {'na', 'sonorant', 'obstruent'});
histogram(C, 'BarWidth', 0.4, 'FaceColor','k');
xlim([categorical({'sonorant'}), categorical({'obstruent'})]);
ylabel('count'); xlabel('sound type at onset')
set(gca, 'FontSize', 13);

%% functions


function kshist(x, bins, col)

    if nargin<3, col = 'k'; end
    histogram(x, bins, 'Normalization', 'Probability', ...
        'EdgeColor','none', 'FaceAlpha', 0.3, 'FaceColor',col); hold on;
    [f, xi]=ksdensity(x, 'bandwidth', 0.28);
    yyaxis right
    plot(xi, f, 'Color',col, 'LineWidth',2.5);
    yticks([]);

end

% function [colors] = getColors(type)
%     switch type
%         case 1 % language type
%             % spanish, english, bilingual mandarin
%             colors = brewermap(4, 'Dark2');
%             colors(4, :) = [35 170 225]./256;
%         case 2 % primary feature encoding map
%             colors = brewermap(3, 'Dark2');
%         case 3 % word syllable distinction
%             colors = [0.6 0 0.6; 0.1 0.7 0.2];
%     end
% end

function [nancol, nanrow, X] = getMaximalResp(X, factor)
        
        % remove columns and rows where data is missing
        % trials for which over electrodes do not have data
        % remanining electrodes with missing trials
        % A is electrode/time x trials
        nancol = sum(isnan(X))>factor*size(X, 1);
        X(:, nancol) = [];
        nanrow = sum(isnan(X), 2)>0;
        X(nanrow, :) = [];
end

function subjs = getSubjectLabels(ls)

    lsorder = {'Spanish', 'English', 'Mandarin', 'Bilingual'};
    subjs = cell(length(ls), 1);

    ctr = 1;
    for l = ls
        subjs(ctr) = lsorder(l);
        ctr = ctr+1;
    end
end


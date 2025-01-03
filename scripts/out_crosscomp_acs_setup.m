%% Set up
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath brewer
addpath(genpath('util'))
zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;
addpath(genpath(datapath))

bef=50;
aft=50;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs] = getSIDinfo();
% SIDs=sSIDs;

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
asccd_details = load('out_sentence_details_asccd_all_loudness.mat');
tps = 50:55;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% ----------------------- Corpus phoneme counts ----------------------- %%
sd = {timit_details.sentdet, dimex_details.sentdet};
load('util/arpabet2ipa.mat', 'arpabet', 'ipa');
timit_details.mapping = [arpabet; ipa]';
clear ipa arpabet

% load in dimex phoneme mapping
load('util/mexbet2ipa.mat', 'mexbet', 'ipa');
dimex_details.mapping = [mexbet; ipa]';
clear ipa mexbet

allphn = {{}, {}};
for i = 1:length(sd{1})    
    allphn{1} = [allphn{1} convertToIPA(timit_details.mapping, ...
        sd{1}(i).phnnames)];
    % for dimex, getting phn names is more complicated
    % count total number of phones in the corpus
    totalphn = sum(sd{2}(i).phnmatonset, 2);
    phnnames = {};
    for j = 1:length(totalphn)
        phnnames = [phnnames repmat(dimex_details.phnnames(j), ...
            1, totalphn(j))];
    end
    
    allphn{2} = [allphn{2} convertToIPA(dimex_details.mapping, ...
        phnnames)];
    clear phnnames j totalphn
end

% remove all closure annotations
allphn{1}(endsWith(allphn{1}, 'cl')) = [];
allphn{2}(endsWith(allphn{2}, '_c')) = [];

[tbl_timit,~, ~,labels_timit] = crosstab(allphn{1});
[tbl_dimex,~, ~,labels_dimex] = crosstab(allphn{2});

% Calculate the relative counts and labels for overlapping labels
overlapLabels = intersect(labels_timit, labels_dimex);
overlapCounts_timit = tbl_timit(ismember(labels_timit, overlapLabels));
overlapCounts_dimex = tbl_dimex(ismember(labels_dimex, overlapLabels));
overlapRelCounts_timit = overlapCounts_timit ./ sum(overlapCounts_timit);
overlapRelCounts_dimex = overlapCounts_dimex ./ sum(overlapCounts_dimex);

% Calculate the relative counts and labels for TIMIT-only labels
timitOnlyLabels = setdiff(labels_timit, labels_dimex);
timitOnlyCounts = tbl_timit(ismember(labels_timit, timitOnlyLabels));
timitOnlyRelCounts = timitOnlyCounts ./ sum(timitOnlyCounts);

% Calculate the relative counts and labels for DIMEX-only labels
dimexOnlyLabels = setdiff(labels_dimex, labels_timit);
dimexOnlyCounts = tbl_dimex(ismember(labels_dimex, dimexOnlyLabels));
dimexOnlyRelCounts = dimexOnlyCounts ./ sum(dimexOnlyCounts);

% Plotting
figure;

% Subplot for overlapping labels
subplot(1, 3, 1);
h = bar([overlapRelCounts_timit overlapRelCounts_dimex]);
h(1).FaceColor = 'red';
h(2).FaceColor = 'blue';
legend('TIMIT', 'DIMEX');
xlabel('Labels')
ylabel('Relative Counts')
title('Overlapping Labels');
box off;
set(gca, 'XTick', 1:length(overlapLabels), 'XTickLabel', overlapLabels);

% Subplot for TIMIT-only labels
subplot(1, 3, 2);
bar(timitOnlyRelCounts, 'FaceColor','red');
xlabel('Labels')
ylabel('Relative Counts')
title('TIMIT-Only Labels');
box off;
set(gca, 'XTick', 1:length(timitOnlyLabels), 'XTickLabel', timitOnlyLabels)

% Subplot for DIMEX-only labels
subplot(1, 3, 3);
bar(dimexOnlyRelCounts, 'FaceColor','blue');
xlabel('Labels')
ylabel('Relative Counts')
title('DIMEx-Only Labels')
box off;
set(gca, 'XTick', 1:length(dimexOnlyLabels), 'XTickLabel', dimexOnlyLabels)

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% ----------------------- Sentence starting characteristics ----------- %%

figure;

% Subplot 1: Plotting loudness for TIMIT sentences
subplot(2, 2, 1);
for i = 1:400
    plot(timit_details.sentdet(i).loudness, 'Color', [0.5 0.5 0.5]); % Plot loudness
    hold on;
end
xline(50, 'LineWidth', 2); % Add a vertical line at x = 50
title('TIMIT');
xlim([30 80]); % Set the x-axis limits

% Subplot 2: Plotting loudness for DIMEx sentences
subplot(2, 2, 2);
for i = 1:400
    % Plot loudness
    plot(dimex_details.sentdet(i).loudness, 'Color', [0.5 0.5 0.5]); 
    hold on;
end
xline(50, 'LineWidth', 2); % Add a vertical line at x = 50
title('DIMEx');
xlim([30 80]); % Set the x-axis limits

% Find the first peak
timit_locs = nan(1, 400);
dimex_locs = nan(1, 400);
for i = 1:400
    % Find peaks in the difference of loudness
    [~, locs] = findpeaks(diff(timit_details.sentdet(i).loudness)); 
    timit_locs(i) = locs(2); % Store the second peak location

    % Find peaks in the difference of loudness
    [~, locs] = findpeaks(diff(dimex_details.sentdet(i).loudness)); 
    dimex_locs(i) = locs(2); % Store the second peak location
end

% Subplot 3: Histogram of TIMIT peak locations
% Subplot 4: Histogram of DIMEx peak locations
subplot(2, 2, [3 4]);
% Plot histogram of TIMIT peak locations
histogram(timit_locs, 'FaceColor', [0.6 0.01 0.01], 'EdgeColor', 'none'); 
hold on;

% Plot histogram of DIMEx peak locations
histogram(dimex_locs, 'FaceColor', [0.01 0.01 0.9], 'EdgeColor', 'none'); 
legend('TIMIT', 'DIMEx')
title('Corpus peak locations')

% Display the average difference between TIMIT and DIMEx onsets
disp(['Average difference between TIMIT/DIMEx onsets: ' ...
    num2str(mean(timit_locs)-mean(dimex_locs))]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% ----------------------- Modulation power spectra -------------------- %%

corpora = {timit_details, dimex_details, asccd_details};
corpuslabels = {'timit', 'dimex',  'asccd'};

figure;
for c = 1:length(corpora)

    % Load the corpus details
    corpus_details = corpora{c};
    corpus = corpuslabels{c};

    % Create two NaN matrices of size 450x4096
    p = nan(450, 4096);
    d = nan(450, 4096);

    % Create a NaN matrix of size 100x257x512
    mtf = nan(100, 257, 512);
    ctr = 1;

    % Iterate over the first 200 sentences
    for i = 1:58 % minimum number of sentences
        fs = 16000;
        % fs = 44100;

        % Calculate the power spectrum of the sound in the i-th sentence
        [p(i, :),~] = pspectrum(corpus_details.sentdet(i).sound, fs);

        % Define start positions for segmentation
        strts = 50:60:size(corpus_details.sentdet(i).aud, 2)-120;

        % Add F0 (fundamental frequency) for the sentence from corpus_details
        f0 = addF0(corpus_details.sentdet(i), corpus);

        % Remove zero values from f0
        f0(f0==0) = [];

        % Check if the mean of f0 is greater than 170
        if mean(f0) > 170
            % Iterate over the start positions
            for j = 1:length(strts)

                % Extract a segment of audio data and compute the MTF
                segaud = corpus_details.sentdet(i).aud(:, strts(j):strts(j)+60);
                [mtf(ctr, :, :), xax, yax, btm, bsm] = get_MTF_bsm_btm(segaud, ...
                    corpus_details.sentdet(i).dataf, 8);

                % Increment the counter
                ctr = ctr + 1;
            end
        end
    end

    % Display the average MTF as an image
    subplot(1, 3, c);
    imagesc(xax, fliplr(yax), squeeze(mean(mtf, 'omitnan')));
    title(corpus);
    hold on;
    set(gca, 'Ydir', 'normal');
    xlim([-8 8]);
    ylim([0, 0.75]);
    colormap(inferno);
    brighten(0.6)

    % Add contour lines to the image
    % contour(xax, fliplr(yax), squeeze(mean(mtf)), 'Color', [0 0 0], 'LineWidth', 1.75);

    % Set labels and font size
    ylabel('Spectral modulation (cycles/oct)');
    xlabel('Temporal modulation (Hz)');
    set(gca, 'FontSize', 15);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% ----------------------- Example sentence with envelope / syllable --- %%

%sents = {[42, 27, 186], [333, 400 450]};
%sents = {[42, 27, 186, 60, 40, 150], [33, 57, 98, 110, 130, 144]};
x = find(strcmp({timit_details.sentdet.name}, 'fcaj0_si1804'));

% find a sentence you know the words for
wordlists = cellfun(@(x) join(x, ' '), {timit_details.sentdet.wordList});
x = find(strcmp(wordlists, 'the wagons were burning fiercely'));

sents = {42, x, 1};

% Call makeSurprisal function for TIMIT and DIMEx details
[timit_details.sentdet] = makeSurprisal(timit_details.sentdet, 8, 'timit');
[dimex_details.sentdet] = makeSurprisal(dimex_details.sentdet, 8, 'dimex');

% Assign corpus details to a cell array
corpus_details = {dimex_details, timit_details, asccd_details};
corpuslabels = {'dimex', 'timit', 'asccd'};

% Load mexbet2ipa mapping
load('mexbet2ipa.mat');

% Create a new figure
figure;

% Set the x-axis limits for the plots
xl = [0, 1.8];

% Iterate over the two corpora (DIMEx and TIMIT)
numcorp = 3;
for c = 1:numcorp
    details = corpus_details{c};

    disp(details.sentdet(sents{c}).name)
    
    % Iterate over the selected sentences for the current corpus
    for s = 1:length(sents{c})
        % Setup for current sentence
        sent = sents{c}(s);
        sent_info = details.sentdet(sent);

        % save the audio of the file to the current folder
        % save only the first 2 seconds
        audiowrite(['example_' sent_info.name '_' corpuslabels{c} '.wav'], ...
            details.sentdet(sent).sound(1:2.5*details.sentdet(sent).soundf), ...
            details.sentdet(sent).soundf);

        % Plot amplitude (waveform) of the sound
        subplot(5, length(corpus_details), c);        
        x = ((1:length(sent_info.sound))/sent_info.soundf)-0.5;
        if startsWith(sent_info.name, 'f')
            factor = 1;
        else
            factor = 3;
        end
        plot(downsample(x, 3), downsample(sent_info.sound, 3)*factor, ...
            'Color', [0.8 0.8 0.8]);
        hold on;

        % Resample x
        x = ((downsample(1:length(sent_info.sound), ...
            sent_info.soundf./sent_info.dataf))/sent_info.soundf)-0.5;
        offby = length(x)-length(sent_info.loudness);
        x = x(1+ceil(offby/2):end-floor(offby/2));

        % Plot amplitude envelope
        plot(x, sent_info.loudness, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
        
        % Formatting for the current subplot
        xticks(0:2);
        xlim(xl)
        ylim([0 1])
        box off;
        xlabel('Time (s)');
        set(gca, 'FontSize', 15);

        % Plot the spectrogram
        ax = subplot(5, numcorp, numcorp*1+c);
        imagesc(x, 1:81, sent_info.aud);
        set(gca, 'YDir', 'normal');
        colormap(ax, inferno);
        xlim(xl)
        xticks(0:2);

        % Plot the phonetic features
        ax = subplot(5,numcorp, numcorp*2+c);

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
        subplot(5, numcorp, numcorp*3+c);
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

        % Plot word boundaries and envelope
        subplot(5, numcorp, numcorp*4+c);
        plot(x, sent_info.loudness, 'Color', [0.1 0.1 0.1], 'LineWidth', 1.5);
        xline(x(sent_info.syltype>0), 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1.5);
        xline(x(sent_info.wordOns>0), 'LineWidth', 2);
        str = join(sent_info.wordList, ' ');
        title([num2str(sent) ': ' str{:}]);   
        xlim(xl);
        box off;
    end   
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% Example sentence in English

sents = [42];



%% Example sentence in Mandarin

sents = [1];

ctr=0;
xmax = 2.1;
for sent = sents
    sent_info = asccd_details.sentdet(sent);
    voweltimes = sent_info.vowelTimes(1, sent_info.vowelTimes(1, :)<(xmax*sent_info.dataf));
    voweltimes = voweltimes/sent_info.dataf - 0.5;

    % Plot amplitude (waveform) of the sound
    subplot(3, length(sents), 1+ctr);        
    x = ((1:length(sent_info.sound))/sent_info.soundf)-0.5;
    plot(downsample(x, 3), downsample(sent_info.sound, 3)*2.5, ...
        'Color', [0.8 0.8 0.8]);
    hold on;

    % Resample x
    x = ((downsample(1:length(sent_info.sound), ...
        sent_info.soundf./sent_info.dataf))/sent_info.soundf)-0.5;
    offby = length(x)-length(sent_info.loudness);
    x = x(1+ceil(offby/2):end-floor(offby/2));
    
    % Plot amplitude envelope
    plot(x, sent_info.loudness, 'Color', [0.3 0.3 0.3], 'LineWidth', 1.5);
    xline(voweltimes, 'LineWidth', 2, 'Color', 'k');

    % Formatting for the current subplot
    xticks(0:2);
    xlim([0 xmax])
    ylim([0 1])
    box off;
    xlabel('Time (s)');
    set(gca, 'FontSize', 15);

    % Plot the spectrogram
    ax = subplot(3, length(sents), 2+ctr);
    imagesc(x, 1:81, sent_info.aud);
    box off;
    set(gca, 'YDir', 'normal');
    colormap(ax, flipud(gray));
    xlim([0 xmax])
    xticks(0:2);
    xline(voweltimes, 'LineWidth', 2, 'Color', 'k');

    % Plot the phonetic features
    ax = subplot(3, length(sents), 3+ctr);

    numfeats = length(asccd_details.features.names);
    imagesc(x, 1:numfeats, sent_info.phnfeat(1:numfeats, :));
    yticks(1:numfeats);
    yticklabels(asccd_details.features.names);
    set(gca, 'YDir', 'normal');

    colormap(ax, [1 1 1; 0 0 0]);
    xlim([0 xmax]);
    xticks(0:2);
    xline(voweltimes, 'LineWidth', 2, 'Color', 'k');
    box off;

    phnmatons = sent_info.phnmatonset(:, 1:xmax*sent_info.dataf);
    [phnidx, ~] = find(phnmatons);
    % convert to phn namess
    phns = arrayfun(@(x) asccd_details.phnnames{x}, ...
        phnidx, 'UniformOutput', false);

    % remove 'sil'
    phns(strcmp(phns, 'sil')) = [];

    chars = '大家實驗了不同的方法';
    sgtitle(chars);
    pinyin = 'Dàjiā shì yǎn liǎo bùtóng de fāngfǎ';
    title(pinyin, 'FontWeight', 'normal');
    ctr=ctr+1;
end

%% ----------------------- Speaker rates ------------------------------- %%

% approximate silence is 1 second for all sentences (some jitter in start)
sil = 1;

corpus_details = {dimex_details, timit_details};
legends = {'DIMEx', 'TIMIT'};
colors = [0 0 1; 1 0 0];

phonrate = cell(2, 1);
syllrate = cell(2, 1);
wordrate = cell(2, 1);
durations = cell(2, 1);
envpeakrate = cell(2, 1);

% For both corpora
for corp = 1:2

    sentdet = corpus_details{corp}.sentdet;

    % sentence durations
    durations{corp} = arrayfun(@(x) sentdet(x).duration-sil, 1:length(sentdet));
    
    % phoneme rate
    phonrate{corp} = arrayfun(@(x) sum(sentdet(x).phnmatonset>0, "all") ...
        /durations{corp}(x), 1:length(sentdet));

    % syllable rate
    syllrate{corp} = arrayfun(@(x) sum(sentdet(x).syltype>0) ...
        /durations{corp}(x), 1:length(sentdet));

    % env peak rate rate
    envpeakrate{corp} = arrayfun(@(x) length(findpeaks(sentdet(x).loudness)) ...
        /durations{corp}(x), 1:length(sentdet));

    % word rate
    wordrate{corp} = arrayfun(@(x) length(sentdet(x).wordList) ...
        /durations{corp}(x), 1:length(sentdet));

    % plotting rate comparisons
    subplot(1, 4, 1);
    histogram(phonrate{corp}, DisplayName=legends{corp}, ...
         HandleVisibility='on', EdgeColor='none', Binedges=3:1:25, ...
         FaceColor=colors(corp, :)); hold on; 
    title('Phoneme rate (pps)');

    subplot(1, 4, 2);
    histogram(syllrate{corp}, DisplayName=legends{corp}, ...
        HandleVisibility='on', EdgeColor='none', Binedges=2:0.5:8, ...
        FaceColor=colors(corp, :)); hold on; 
    title('Syllable rate (sps)');

    subplot(1, 4, 3);
    histogram(envpeakrate{corp}, DisplayName=legends{corp}, ...
        HandleVisibility='on', EdgeColor='none', Binedges=2:0.5:6, ...
        FaceColor=colors(corp, :)); hold on; 
    title('Envelope peak rate (eps)');

    subplot(1, 4, 4);
    histogram(wordrate{corp}, DisplayName=legends{corp}, ...
        HandleVisibility='on', EdgeColor='none', Binedges=1:0.25:6, ...
        FaceColor=colors(corp, :)); hold on;
    title('Word rate (wps)');
end

for i = 1:4
    subplot(1, 4, i);
    ylabel('sentence count')
    legend;
    box off;
end

% plot sentence duration differences
figure;
histogram(durations{1}, Binedges=0.7:0.25:5, FaceColor=colors(1, :), ...
    EdgeColor='none'); hold on;
histogram(durations{2}, Binedges=0.7:0.25:5, FaceColor=colors(2, :), ...
    EdgeColor='none'); hold on;
ylabel('sentence count');
xlabel('sentence duration (s)');
legend;
box off;
set(gca, 'FontSize', 15);
legend('DIMEx', 'TIMIT');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% ----------------------- Correlating sentence features cross-corpus -- %%

wind = bef + 1 : bef + 200;

aud = cell(1, 2);
env = cell(1, 2);
feat = cell(1, 2);
surp = cell(1, 2);
word = cell(1, 2);

% Compute surprisal values for TIMIT and DIMEx datasets
[timit_details.sentdet] = makeSurprisal(timit_details.sentdet, 8, 'timit');
[dimex_details.sentdet] = makeSurprisal(dimex_details.sentdet, 8, 'dimex');

lang = 1;
feature_overlap = intersect(timit_details.features.names, ...
    dimex_details.features.names);

for details = {timit_details, dimex_details}
    corpus_details = details{1};

    % Initialize NaN matrices for audio, envelope, feature, surprisal, and word onset
    aud{lang} = nan(size(corpus_details.sentdet, 2), 80, length(wind));
    env{lang} = nan(size(corpus_details.sentdet, 2), length(wind));
    feat{lang} = nan(size(corpus_details.sentdet, 2), ...
        length(feature_overlap), length(wind));
    surp{lang} = nan(size(corpus_details.sentdet, 2), length(wind));

    feat_idx = find(ismember(dimex_details.features.names, feature_overlap));
    for i = 1:size(corpus_details.sentdet, 2)

        % Check if the sentence duration is longer than one second
        if size(corpus_details.sentdet(i).aud, 2) > wind(end) + aft
            aud_tmp = corpus_details.sentdet(i).aud(:, wind);
            loudness_tmp = corpus_details.sentdet(i).loudness(wind);
            phnfeat_tmp = corpus_details.sentdet(i).phnfeatonset(feat_idx, wind);

            if lang == 1
                surp_tmp = corpus_details.sentdet(i).engSurp(wind);
            else
                surp_tmp = corpus_details.sentdet(i).spSurp(wind);
            end

            word_tmp = corpus_details.sentdet(i).wordOns(wind);

            % Find the time of phonetic features
            [~, c] = find(phnfeat_tmp);
            time_phns = sort(unique(c));

            % Store the data in the respective matrices
            aud{lang}(i, :, 1:length(time_phns)) = aud_tmp(:, time_phns);
            env{lang}(i, 1:length(time_phns)) = loudness_tmp(:, time_phns);
            feat{lang}(i, :, 1:length(time_phns)) = phnfeat_tmp(:, time_phns);
            surp{lang}(i, 1:length(time_phns)) = surp_tmp(:, time_phns);
            word{lang}(i, 1:length(time_phns)) = word_tmp(:, time_phns);

            clear *_tmp
        end
    end

    % Reshape and remove NaN rows from the matrices
    aud{lang} = reshape(aud{lang}, [], 80 * length(wind));
    aud{lang}(sum(isnan(aud{lang}), 2) == size(aud{lang}, 2), :) = [];

    feat{lang} = reshape(feat{lang}, [], length(feature_overlap) * length(wind));
    feat{lang}(sum(isnan(feat{lang}), 2) == size(feat{lang}, 2), :) = [];

    lang = lang + 1;
end

figure;
ctr = 1;
features = {aud, env, feat, surp};
numfeat = length(features);
for f = features
    level = f{1};

    % Look at normalized histograms of the feature
    if ctr ~= 3
        subplot(2, numfeat, ctr + numfeat);
        l = level{1}(:);
        l(isnan(l) | l == 0) = [];
        histogram(zscore(l)); hold on;

        l = level{2}(:);
        l(isnan(l) | l == 0) = [];
        histogram(zscore(l));
        box off;
    end

    % TIMIT only correlation
    eng_corr = triu(corr(level{1}', level{1}', 'rows', 'pairwise'), 1);
    eng_corr(eng_corr == 0) = NaN;

    % DIMEx only correlation
    sp_corr = triu(corr(level{2}', level{2}', 'rows', 'pairwise'), 1);
    sp_corr(sp_corr == 0) = NaN;

    % Cross corpus correlation
    cross_corr = triu(corr(level{1}', level{2}', 'rows', 'pairwise'), 1);
    cross_corr(cross_corr == 0) = NaN;

    subplot(2, numfeat, ctr);

    % Compute means and standard errors
    means = [mean(eng_corr(:), 'omitnan'), mean(sp_corr(:), 'omitnan'), ...
        mean(cross_corr(:), 'omitnan')];
    errors = [std(eng_corr(:), 'omitnan')/sqrt(sum(~isnan(eng_corr(:)))), ...
              std(sp_corr(:), 'omitnan')/sqrt(sum(~isnan(sp_corr(:)))), ...
              std(cross_corr(:), 'omitnan')/sqrt(sum(~isnan(cross_corr(:))))];
    
    % Create bar plot with error bars
    bar(means, 'FaceColor',[0.6, 0.6, 0.6], 'EdgeColor','none');
    hold on;
    numgroups = 1;
    numbars = length(means);
    groupwidth = min(0.5, numbars/(numbars + 1.5));
    x = (1:numbars) - groupwidth/2 + groupwidth / (2*numgroups);
    for i = 1:numbars
        errorbar(x(i), means(i), errors(i), 'k', 'linewidth', 1.5);
    end
    box off;
    set(gca, 'XTick', 1:numbars, 'XTickLabels', {'TIMIT only', ...
        'DIMEX only', 'Cross corpus'})

    ctr = ctr + 1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;


%% ----------------------- Average sentence properties ----------------- %%

details = {timit_details, dimex_details};
corpora = {'timit', 'dimex'};

% plot envelope overlaid with start times plotTRFPred
figure;
env_aligned = cell(1, 2);
spec_aligned = cell(1, 2);
phn_aligned = cell(1, 2);
formant_aligned = cell(1, 2);
for corpus = 1:2

    % load in the first 150 ms of the sentence from the sentence onset
    env = cell2mat(arrayfun(@(x) details{corpus}.sentdet(x).loudness(1:149), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false)');
    % load in the last 150 ms of the sentence from the sentence offset
    env = cat(2, env, cell2mat(arrayfun(@(x) details{corpus}.sentdet(x).loudness(end-100:end), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false)'));

    % get the start times of the sentences
    dataf = details{corpus}.sentdet(1).dataf;

    % onset and offset assuming no silence at the beginning
    startps = arrayfun(@(x) round(details{corpus}.sentdet(x).soundOns(1).*dataf), 1:length(details{corpus}.sentdet));
    endtps = arrayfun(@(x) round(details{corpus}.sentdet(x).soundOns(2).*dataf), 1:length(details{corpus}.sentdet));

    % startps and endtps for the stitched sentences 
    startps_stitch = startps+50;
    endtps_stitch = arrayfun(@(x) endtps(x)-length(details{corpus}.sentdet(x).loudness)+300, ...
        1:length(details{corpus}.sentdet));

    % align envelope to 50ms before sentence onset
    env_aligned{corpus} = cell2mat(arrayfun(@(x) details{corpus}.sentdet(x).loudness(startps(x):startps(x)+149), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false)');
    env_aligned{corpus} = cat(2, env_aligned{corpus}, cell2mat(arrayfun(@(x) details{corpus}.sentdet(x).loudness(endtps(x):endtps(x)+99), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false)'));

    % align spectrogram to 50ms before sentence onset
    spec_tmp = arrayfun(@(x) details{corpus}.sentdet(x).aud(:, startps(x):startps(x)+149), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false);
    spec_aligned{corpus} = cat(3, spec_tmp{:});
    spec_tmp = arrayfun(@(x) details{corpus}.sentdet(x).aud(:, endtps(x):endtps(x)+99), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false);
    spec_aligned{corpus} = cat(2,  spec_aligned{corpus}, cat(3, spec_tmp{:}));

    % align phonetic features to 50ms before sentence onset
    phn_tmp = arrayfun(@(x) details{corpus}.sentdet(x).phnfeat(:, startps(x):startps(x)+149), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false);
    phn_aligned{corpus} = cat(3, phn_tmp{:});
    phn_tmp = arrayfun(@(x) details{corpus}.sentdet(x).phnfeat(:, endtps(x):endtps(x)+99), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false);
    phn_aligned{corpus} = cat(2,  phn_aligned{corpus}, cat(3, phn_tmp{:}));

    % align formant features to 50ms before sentence onset
    formant_tmp = arrayfun(@(x) details{corpus}.sentdet(x).formants(:, startps(x):startps(x)+149), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false);
    formant_aligned{corpus} = cat(3, formant_tmp{:});
    formant_tmp = arrayfun(@(x) details{corpus}.sentdet(x).formants(:, endtps(x):endtps(x)+99), ...
            1:length(details{corpus}.sentdet), 'UniformOutput', false);
    formant_aligned{corpus} = cat(2,  formant_aligned{corpus}, cat(3, formant_tmp{:}));
             
    % only sentences 1-200
    sents = 1:200; 

    subplot(2, 2, corpus);
    imagesc(env(sents, :));
    hold on;
    scatter(startps_stitch(sents), 1:length(sents), 15, 'r', 'filled');
    scatter(endtps_stitch(sents), 1:length(sents), 15, 'r', 'filled');
    xline([50, 200], 'LineWidth', 2);

    % Formatting
    colormap(flipud(gray));
    xlim([0 250]);
    ylim([0 length(sents)]);
    ylabel('Sentence');
    xlabel('Time (ms)');
    title(corpora{corpus});

    % plot average envelope
    subplot(2, 2, corpus+2);
    plot(1:250, mean(env, 1)); hold on;
    plot(1:250, mean(env_aligned{corpus}, 1));
    xline([50, 200], 'LineWidth', 2);
    xlim([0 250]);
    ylim([0 1]);
    ylabel('Envelope');
    xlabel('Time (ms)');
    title(corpora{corpus});
end

% plot aligned envelope with each language overlaid
figure;

timit_env = normalize(replaceMat(env_aligned{1}, 0, NaN), 2, 'range');
dimex_env = normalize(replaceMat(env_aligned{2}, 0, NaN), 2, 'range');
plot(1:250, mean(timit_env, 1, "omitnan"), 'LineWidth', 2, 'Color', 'r'); hold on;
plot(1:250, mean(dimex_env, 1, "omitnan"), 'LineWidth', 2, 'Color', 'b');
xline([50, 150, 200], 'LineWidth', 2);
xlim([0 250]);
ylim([0 1]);
ylabel('Norm Envelope');
xlabel('Time (ms)');
title('Aligned Envelope');
legend({'TIMIT', 'DIMEx'});

% plot average spectrogram and phonetic features
figure;
for corpus = 1:2
    % plot average spectrogram
    subplot(2, 2, corpus);
    imagesc(mean(spec_aligned{corpus}, 3));
    set(gca,'YDir','normal');
    colormap(gca, internet);

    % Formatting
    xline([50, 150, 200], 'LineWidth', 2);
    xlim([0 250]);
    ylim([0 50]);
    ylabel('Frequency (Hz)');
    xlabel('Time (ms)');
    title(corpora{corpus});
    caxis([1 5]);

    % plot average phonetic features
    subplot(2, 2, corpus+2);
    imagesc(mean(phn_aligned{corpus}, 3));
    set(gca,'YDir','normal');
    colormap(gca, flipud(gray));

    % Formatting
    xline([50, 150, 200], 'LineWidth', 2);
    xlim([0 250]);
    ylabel('Feature');
    xlabel('Time (ms)');
    yticks(1:length(details{corpus}.features.names));
    yticklabels(details{corpus}.features.names);
    ylim([0.5 14.5]);
    colorbar;
    caxis([0 0.7]);
    box off;
end 

% plot average formant
figure;
colors = [1 0 0; 0 0 1];
for formant = 1:4
    % plot average aligned formant
    subplot(1, 4, formant);
    for corpus = 1:2
        %normform = normalize(formant_aligned{corpus}(formant, :, :), 2, 'range');
        normform = formant_aligned{corpus}(formant, :, :);
        normform = smoothdata(normform, 2, "gaussian", 'SmoothingFactor',0.1);
        plot(mean(normform, 3), 'LineWidth', 2, 'Color', colors(corpus, :)); hold on;
        clear normform
    end

    % Formatting
    xline([50, 150, 200], 'LineWidth', 2);
    xlim([50 200]);
    xlabel('Time (ms)');
    title(['F' num2str(formant)]);
    box off;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% ----------------------- Acoustic classifiers load data -------------- %%

% preprocess both cons structures to include only relevant phoneme classes
remphn = ismember(timit_details.phnnames, {'ax-h', 'axr', 'p', 't', 'k', 'g', ...
    'b', 'd', 'epi', 'dx'});
phns = timit_details.phnnames(~remphn);
[TDcons] = loadDDcons('timit', bef, aft, {}, [], [], timit_details.sentdet, phns);

remphn = ismember(dimex_details.phnnames, {'p', 't', 'k', 'g', 'b', 'd'});
% TODO: Fix the dimex_details voicing category?!
dimex_details.features.voiced = [dimex_details.features.voiced, {'a', 'e', 'i', 'o', 'u'}];
phns = dimex_details.phnnames(~remphn);
[Dcons] = loadDDcons('dimex', bef, aft, {}, [], [], dimex_details.sentdet, phns);
clear phns remphns

details = {timit_details, dimex_details};
Scons = {TDcons, Dcons};
corpnames = {'timit', 'dimex'};

%% ----------------------- Classifier for each phonetic feature -------- %%
% decoding details
timing = 20;
tps = 51-round(timing/2):51+round(timing/2); % around  
% 
nreps = 20;
AUC_all=cell(1, 2);
T_all = cell(2, 12);
featnames_all=cell(1, 2);

% pre-loaded data and prediction matrix for cross-language test
trlidx = cell(3, 2);
X = cell(3, 2);
y = cell(3, 2);
phns = cell(3, 2);
trialset = {{'dorsal', 'coronal', 'labial', 'plosive', 'fricative', 'nasal'}, ...
    {'high', 'front', 'low', 'back'}, {'sonorant', 'voiced', ...
    'obstruent'}};

for corp = 1:2

    % load in all phoneme labels
    Scon = Scons{corp};

     % get phn feat matrix 
    phnfeat = arrayfun(@(x) details{corp}.features.mat(:, x), Scon.phnType, ...
        'UniformOutput', false); 
    phnfeat = cat(2, phnfeat{:});

     % remove obstruent (since sonorant is almost exact opposite)
    featrem = {'obstruent', 'approximant', 'syllabic'};
    featnames = details{corp}.features.names( ...
        ~ismember(details{corp}.features.names, featrem));
    phnfeat(ismember(details{corp}.features.names, featrem), :) = [];

    if corp == 2
        % manually add voicing to vowels in dimex
        vidx = find(strcmp(featnames, 'voiced'));
        vowelidx = find(ismember(Scon.phn, {'a', 'e', 'i', 'o', 'u'}));
        phnfeat(vidx, vowelidx) = 1;
    end

    % construct different samples for consonants and vowels
    for trltype = 1:3

        % select 5000 samples randomly from the set of all phns such that
        % at least one feature (on) is represented per phn 
        rng(1);
        tmp = cellfun(@(x) details{corp}.features.(x), ...
           trialset{trltype}, 'UniformOutput', false);
        tidx = find(ismember(Scon.phn, unique(cat(2, tmp{:}))));
        trlidx{trltype, corp} = randsample(tidx, 2000);
            
        X{trltype, corp} = reshape(Scon.aud(:, tps, ...
            trlidx{trltype, corp}), 80*length(tps), []);
        y{trltype, corp} = phnfeat(:, trlidx{trltype, corp});
        phns{trltype, corp} = Scon.phn(trlidx{trltype, corp});
    end
end

for corp = 1:2 
    disp(corpnames{corp})
    
    acc = nan(length(featnames), nreps); % test accuracy
    weights = cell(length(featnames), 1); % aud weights
    for feat = 1:length(featnames)
        % find trialtype
       trialtype = find(cellfun(@(x) ismember(featnames{feat}, x), trialset));
    
        % same language comparison -- TRAINED ON LANG 1, TESTED ON LANG 1
        [~, ~, AUC_all{corp}(1, feat, :), ~, ~, acc(feat, :), ...
            weights{feat}, ~, Mdl, comp, mu] = logistic(X{trialtype, corp}', ...
            y{trialtype, corp}(feat, :)', 1, [], tps, nreps); 
        
        % cross language comparison -- TRAINED ON LANG 1, TESTED ON LANG 2
        crosscorp = mod(corp, 2)+1;
        [AUC_all{crosscorp}(2, feat, :), pcaX, ypred] = crosstestlogistic(X{trialtype, crosscorp}', ...
            y{trialtype, crosscorp}(feat, :), Mdl, comp, mu, nreps);

        disp(['Computing acoustic clarity score for ' ...
            featnames{feat} '... median AUC = ' num2str(median(AUC_all{corp}(1, feat, :))) ])
        disp(['Cross-language... median AUC = ' num2str(median(AUC_all{crosscorp}(2, feat, :))) ]);

        % look at which classes are the most misclassified
        %if strcmp(featnames{feat}, 'dorsal')            
            mislabel = cell2table(tabulate(phns{trialtype, crosscorp}(~ypred)), ...
                'VariableNames', {'Phn','Count','Percent'});
            alllabel = cell2table(tabulate(phns{trialtype, crosscorp}), ...
                'VariableNames', {'Phn','Count','Percent'});
            T = outerjoin(mislabel, alllabel, 'Keys', 'Phn');
            T.mislabeldivall = T.Count_mislabel./T.Count_alllabel;

            % remove all 
            T(T.Count_alllabel<10, :) = [];
            T.label = cellfun(@(x) ismember(x, ...
                details{crosscorp}.features.(featnames{feat})), T.Phn_alllabel);
            T_all{corp, feat} = sortrows(T, 7);
       %end
    end
    featnames_all{corp} = featnames;
end

% plotting
tmp = brewermap(5, 'YlGnBu');
tmp2 = brewermap(5, 'YlOrRd');

% plot the auc values
cols = [{tmp2([4, 2], :)}, {tmp([4, 2], :)}];

% plot the auc values
figure; 
auc = nan(2, nreps);
p = nan(2, length(featnames));
for ls = 1:2
    for feat = 1:length(featnames)
        for predtype = 1:2 % 1 is same language prediction, 2 is cross-language
            auc(predtype, :) = squeeze(AUC_all{ls}(predtype, feat, :));
            boxplot(auc(predtype, :), 'Position', ...
                feat-0.5+0.45*ls+0.15*predtype, 'Color', ...
                cols{ls}(predtype, :), 'BoxStyle', ...
                'filled', 'OutlierSize', 0.01); hold on;
        end
        p(ls, feat) = ranksum(auc(1, :), auc(2, :));

        % bonferroni corrected!
        text(feat+0.5, 1-0.1*ls, getSigStr(p(ls, feat), 3, 26), 'FontSize', 15);    
    end
end
xticks(1:length(featnames))
xticklabels(featnames);


% plotting

featgroups = {{'dorsal', 'coronal', 'labial', 'plosive', 'fricative', 'nasal'}, ...
    {'high', 'front', 'low', 'back'}, {'sonorant', 'voiced'}};
featitles = {'consonant', 'vowel', 'sonority'};
figure;

cols = [{tmp2([4, 2], :)}, {tmp([4, 2], :)}];
for trialtype = 1:3
    feats = featgroups{trialtype};
    fidx = cellfun(@(x) find(strcmp(featnames, x)), feats);

    cdata = [median(AUC_all{1}(:, fidx, :), [2, 3]), ...
        median(AUC_all{2}(:, fidx, :), [2, 3]);]';
    err = [nansem(median(AUC_all{1}(:, fidx, :), 3), 2), ...
        nansem(median(AUC_all{2}(:, fidx, :), 3), 2);]';
    subplot(3, 1, trialtype);
    b = bar(cdata, 'FaceColor', 'none'); hold on;
    
    % Calculate the number of groups and number of bars in each group
    [ngroups,nbars] = size(cdata);

    % Get the x coordinate of the bars
    x = nan(nbars, ngroups);
    for i = 1:nbars
        x(i,:) = b(i).XEndPoints;
    end
    % Plot the errorbars
    errorbar(x',cdata,err,'k','linestyle', ...
        'none', 'LineWidth', 1, 'CapSize',0.5);
    hold off

    hatchfill2(b(2),'cross','HatchAngle',45, ...
        'hatchcolor',[0 0 0],'HatchDensity', 20);
    
    set(gca, 'FontSize', 13, 'Xtick', 1:2, 'Xticklabel', ...
        {'English', 'Spanish'});
    ylabel('AUC')
    title(featitles{trialtype}, 'FontWeight', 'normal');
    ylim([0.5 1]);
    yticks([0.5 1]);
    box off; 
end
xlabel('Trained on');
legend({'Same-language', 'Cross-language'});

% plotting pt. 2

featgroups = {{'dorsal', 'coronal', 'labial', 'plosive', 'fricative', 'nasal'}, ...
    {'high', 'front', 'low', 'back'}, {'sonorant', 'voiced'}};
for trialtype = 1:3
    feats = featgroups{trialtype};
    fidx = cellfun(@(x) find(strcmp(featnames, x)), feats);

    cdata = [diff(flipud(median(AUC_all{1}(:, fidx, :), 3))); ...
        diff(flipud(median(AUC_all{2}(:, fidx, :), 3)));]';
    figure;
    imagesc(cdata); 
    
    for ls = 1:2
        for feat = 1:length(feats)
            if p(ls, fidx(feat))<(0.01/26)
                fontweight = 'bold';
            else
                fontweight = 'normal';
            end
    
            if cdata(feat, ls)<0.09
                text(ls-0.2, feat, num2str(cdata(feat, ls), 2), 'FontSize', 13, ...
                    'FontWeight', fontweight); hold on;
            else
                text(ls-0.2, feat, num2str(cdata(feat, ls), 2), 'FontSize', 13, ...
                   'FontWeight', fontweight, 'Color', 'w'); hold on;
            end
        end         
    end
    
    set(gca, 'FontSize', 15, 'Ytick', 1:length(feats), ...
        'Yticklabel', feats, 'Xtick', 1:2, 'Xticklabel', ...
        {'English speech', 'Spanish speech'});
    xlabel('Trained on')
    box off; 
    colormap(flipud(gray));
    cbh = colorbar;
    cbh.Ticks = [0, 0.3];
    cbh.Location = "northoutside" ;
    caxis([0 0.25]);
    ylabel(cbh, 'Tested on: Same language - Different language');
end

% figure;
% barh(cdata);
% set(gca, 'YDir', 'normal')

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections *cons;

%% ----------------------- Voiced phoneme durations -------------------- %%

histogram(diff(Dcons.phnTimes(:, ...
    ismember(Dcons.phn, dimex_details.features.voiced))/100), 'BinWidth', 0.01, ...
    'Normalization','probability'); 
hold on;
histogram(diff(TDcons.phnTimes(:, ...
    ismember(TDcons.phn, timit_details.features.voiced))/100), 'BinWidth', 0.01, ...
    'Normalization','probability'); 
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% Other ideas on quantifying language similarity...

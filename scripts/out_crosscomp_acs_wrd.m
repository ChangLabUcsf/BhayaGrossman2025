%% Set up

% Ilina Bhaya-Grossman
% 01.08.2022
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('util'))

zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs,bSIDs, mSIDs] = getSIDinfo();
SIDs = [sSIDs, eSIDs, {'HS11', 'HS9', 'HS10'}];

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
% asccd_details = load('stim_info/out_sentence_details_acssd_loudness.mat');

bef=50;
aft=50;

% eSIDs sSIDs
Dwrd = loadDwrdMulti('dimex',  bef, aft, [], dimex_details); % {'EC260', 'EC266'}
TDwrd = loadDwrdMulti('timit', bef, aft, [], timit_details); % {'EC260', 'EC266'}

window = 20;
fieldname = 'aud';
Dwrd.ambiguity = getAAI(window, Dwrd, fieldname);
TDwrd.ambiguity = getAAI(window, TDwrd, fieldname);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% --------------- ACOUSTIC EXPLORATION OF WORD BOUNDARY ------------------
%% Visualize the difference between syllable and word onsets

Swrds = {Dwrd, TDwrd}; % Contains the variables Dwrd and TDwrd

figure;
for s = Swrds

    % For every sentence, aggregate the consecutive differences between syllables and words
    Swrd = s{1};
    syl2syl = [];
    wrd2wrd = [];
    cr = [];

    for si = unique(Swrd.sentid)'
        sid = si{1};
        tmp = Swrd(strcmp(Swrd.sentid, sid), ["tps" "wordOns"]);        

        % Find all consecutive 0-1 or 1-0 transitions (syllable/word to word/syllable)
        cidx = [find(diff(tmp.wordOns)) find(diff(tmp.wordOns))+1]; 
        cr = [cr arrayfun(@(x) diff(tmp.tps(cidx(x, :))), 1:size(cidx, 1))];      

        % Find all consecutive 0-0 (multisyllable) or 1-1 (word-word) transitions
        sidx = find(diff(tmp.wordOns) == 0); % All cases where consecutive values are the same
        % 0-0 transitions (syllable to syllable)
        syllidx = [sidx(tmp.wordOns(sidx) == 0) sidx(tmp.wordOns(sidx) == 0)+1]; 

        % 1-1 transitions (word to word)
        wrdidx = [sidx(tmp.wordOns(sidx) == 1) sidx(tmp.wordOns(sidx) == 1)+1]; 

        % Populate the syllable to syllable and word to word arrays
        syl2syl = [syl2syl arrayfun(@(x) diff(tmp.tps(syllidx(x, :))), 1:size(syllidx, 1))];
        wrd2wrd = [wrd2wrd arrayfun(@(x) diff(tmp.tps(wrdidx(x, :))), 1:size(wrdidx, 1))];

        clear cidx tmp sid
    end

    subplot(1, 3, 1);
    histogram(cr, 20, 'Normalization', 'probability', 'EdgeColor', 'none'); hold on;

    subplot(1, 3, 2);
    histogram(syl2syl, 20, 'Normalization', 'probability', 'EdgeColor', 'none'); hold on;

    subplot(1, 3, 3);
    histogram(wrd2wrd, 20, 'Normalization', 'probability', 'EdgeColor', 'none'); hold on;
end

titles = {'Cross (syll-word, word-syll)', 'Syll to Syll', 'Wrd to Wrd'};
for i = 1:3
    subplot(1, 3, i);
    xlabel('Time between events (10ms)');
    ylabel('Probability');
    set(gca, 'FontSize', 15);
    title(titles{i});
    xlim([0 80]);
    legend({'Spanish', 'English'});
    box off;
    xline(5, 'HandleVisibility', 'off', 'LineWidth', 1.5, 'LineStyle', '--');
end

clear Swrd wrd2wrd 

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% Visualize average envelope for words and syllables

befaft = [0.5 0.5];
x = -befaft(1):0.01:befaft(2);

colors = getColorsCrossComp(3);
addpath(genpath('shadederror/'))

figure;
ctr=1;
for i = {TDwrd, Dwrd}
    subplot(1, 2, ctr)
    ywrd = squeeze(cat(3, i{1}.env{logical(i{1}.wordOns)}))';
    ysyll = squeeze(cat(3, i{1}.env{~i{1}.wordOns}))';
    
    
    % Plot the average envelope for word and syllable
    shadedErrorBar(x, mean(ysyll, 1), nansem(ysyll, 1), ...
        'lineprops', {'color', colors(1, :), 'linewidth', 1.8});
    shadedErrorBar(x, mean(ywrd, 1), nansem(ywrd, 1), ...
        'lineprops', {'color', colors(2, :), 'linewidth', 1.8}); hold on;

    % Calculating the F-stat values over time
    [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
        reshape([ysyll; ywrd]', 1, size([ysyll; ywrd], 2), size([ysyll; ywrd], 1)), ...
        [ones(1, size(ysyll, 1)) 2*ones(1, size(ywrd, 1))], [1, 2]);
    fthresh = finv(1-0.0001, df1, df2);  

    % Visualizing the F-stat values on average envelope plot
    scatter(x(fvals>fthresh), 0.1*ones(1, sum(fvals>fthresh)), 15, ...
        fvals(fvals>fthresh), 'filled');
    cm = colormap("gray");
    colormap(flipud(cm(1:200, :)));

    % Reference lines and formatting
    xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
    ylabel('norm envelope'); 
    xlabel('time from ons (s)');
    legend({'syllable', 'word'});
    set(gca, 'FontSize', 13);
    ylim([-0.01 1]);

    ctr=ctr+1;
end

%
figure;
ctr=1;

befaft = [0.5 0.5];
x = -befaft(1):0.01:befaft(2);

colors = getColorsCrossComp(3);
addpath(genpath('shadederror/'))

for i = {TDwrd, Dwrd}
    subplot(1, 2, ctr)
    ywrd = squeeze(cat(3, i{1}.env{logical(i{1}.wordOns) & ~i{1}.sentOns}))';
    ysyll = squeeze(cat(3, i{1}.env{~i{1}.wordOns & ~i{1}.sentOns}))';
    
    % Plot all envelopes for word and syllable
    h=plot(x, ysyll(1:40, :), 'color', colors(1, :), ...
        'linewidth', 0.1); hold on;
    for j = 1:40, h(j).Color(4) = 0.9; end

    h=plot(x, ywrd(1:40, :), 'color', colors(2, :), ...
        'linewidth', 0.1); hold on;
    for j = 1:40, h(j).Color(4) = 0.4; end

    % Reference lines and formatting
    xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
    ylabel('norm envelope'); 
    xlabel('time from ons (s)');
    legend({'syllable', 'word'});
    set(gca, 'FontSize', 13);
    xlim([-0.2 0.2])
%     ylim([-0.01 1]);

    ctr=ctr+1;
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;


%% Visualize envelope for words and syllables with varying length

befaft = [0.5 0.5];
x = -befaft(1):0.01:befaft(2);
addpath(genpath('shadederror/'))

for i = {TDwrd, Dwrd}

    allsent_tbl = i{1}; 

    ywrd = squeeze(cat(3, i{1}.env{logical(i{1}.wordOns)}))';
    ysyll = squeeze(cat(3, i{1}.env{~i{1}.wordOns}))';
    
    % remove all envelope instances that are naned
    envidx = cellfun(@(x) ~isempty(x), [allsent_tbl.env]);
    
    % words
    widx = envidx & allsent_tbl.wordOns;
    sidx = envidx & ~allsent_tbl.wordOns;
    
    % stressed and unstressed syllables
    % y_str = cell2mat([allsent_tbl.env(envidx & allsent_tbl.stress==2 ...
    %     & allsent_tbl.wordOns==0)]);
    % y_unstr = cell2mat([allsent_tbl.env(envidx & allsent_tbl.stress==1 ...
    %     & allsent_tbl.wordOns==0)]);

    % Visualize effect of preceding word length on envelope
    numbins = 3;
    colors = getColorsCrossComp(2);
    len_mean = nan(3, 1);
    
    figure;
    bins = [0 prctile(allsent_tbl.precWordLen, 33), ...
        prctile(allsent_tbl.precWordLen, 66), 20]+1;
    allplen = discretize(allsent_tbl.precWordLen, bins); % [0 2 3 6]
    for len = 1:numbins
        precwrd = allplen == len & widx;
        if sum(precwrd)>20
            y = cell2mat(allsent_tbl.env(precwrd));
            shadedErrorBar(x, mean(y, 1), nansem(y, 1), ...
                'lineprops', {'Color', colors(len, :), 'linewidth', 2});
            len_mean(len, 1) = round(mean(allsent_tbl.precWordLen(precwrd)));
        end
        clear y precwrd 
    end

    % legend({'1 syll', '2 syll', '3+ syll'});
    legend({['short: ' num2str(len_mean(1))] , ...
        ['medium: ' num2str(len_mean(2))], ['long: ' num2str(len_mean(3))]});
    xline(0, 'LineWidth', 2, 'LineStyle', '--', ...
        'Color','k', 'HandleVisibility', 'off');
    set(gca, 'FontSize', 13);

    % Average envelope for word and syllable overlaid on single instances
    colors = getColorsCrossComp(3);
    figure;
    ax = subplot(2, 1, 1);
    numint = 100;
    cdata = normr(ysyll(1:numint, :));
    imagesc('Xdata', x, 'Cdata', smoothdata(cdata, 1)); hold on;

    set(gca, 'FontSize', 13);
    xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
    ylim([0 numint]);
    xlim([x(1) x(end)]);

    title('Syllable Onset');
    colormap(ax, brighten(flipud(gray), 0.6));
    
    yyaxis right
    plot(x, smoothdata(mean(cdata), 'SmoothingFactor', 0.15), ...
        'Color', [0 0 0], 'linewidth', 3);    
    yticks([]);
    
    ax = subplot(2, 1, 2);
    [syls, sortidx] = sort(allsent_tbl.precWordSyl(widx), 'descend');
    sortidx(isnan(syls)|allsent_tbl.precWordLen(sortidx)==0)=[];
    cdata = normr(ywrd(sortidx(1:numint), :));
    imagesc('Xdata', x, 'Cdata', smoothdata(cdata, 1)); hold on;
    
    % formatting
    ylim([0 numint]);
    xlim([x(1) x(end)]);
    xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
    title('Word Onset')
    xlabel('time from onset (s)');
    ylabel('trials'); 
    set(gca, 'FontSize', 13);
    colormap(ax, brighten(flipud(fpurple), 0.6));

    yyaxis right
    plot(x, smoothdata(mean(cdata)), ...
        'Color', colors(2, :), 'linewidth', 3);   
    yticks([]);
end
    
clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% Visualizing onset and offset phonemes

labels = {'onset', 'offset'};
corpus = {'TIMIT', 'DIMEx'};

% Load conversion tables for ARPAbet to IPA phonemes
timit_conv = load("arpabet2ipa.mat");
timit_conv.bet = timit_conv.arpabet;
dimex_conv = load('mexbet2ipa.mat');
dimex_conv.bet = dimex_conv.mexbet;
convs = {timit_conv, dimex_conv};

% Access corpus details for TIMIT and DIMEx
corpus_details = {timit_details, dimex_details};

ctr = 1;
for i = {TDwrd, Dwrd}
    allsent_tbl = i{1}; % DIMEx
    
    % Extract onset and offset phonemes from the word table
    nonempty = cellfun(@(x) size(x, 2)>0, allsent_tbl.phns);
    onset_phns = cellfun(@(x) x(1), allsent_tbl.phns(allsent_tbl.wordOns>0 & nonempty));
    offset_phns = cellfun(@(x) x(end), allsent_tbl.phns(allsent_tbl.wordOns>0 & nonempty));    
    phntype = {onset_phns, offset_phns};

    figure;
    for t = 1:2
        subplot(2, 4, (t-1)*4+[1 2 3]);

        % Map ARPAbet to IPA phonemes based on the corpus type
        ipapp = {'b', 'p', 'g', 'k', 't', 'd', 'h', '', '', 'ə'};
        if startsWith(allsent_tbl.sentid(1), 's')
            append = {'b_c', 'p_c', 'g_c', 'k_c', 't_c', 'd_c', 'h', 'epi', '', ''};
        else
            append = {'bcl', 'pcl', 'gcl', 'kcl', 'tcl', 'dcl', 'hv', 'epi', 'eng', 'ax-h'};
        end

        arpabet = [convs{ctr}.bet append];
        ipa = [convs{ctr}.ipa ipapp];
        ipas = cellfun(@(x) join(string(ipa(strcmp(arpabet, x))), ''), phntype{t}, ...
            'UniformOutput', false);
        C = categorical(cellstr(ipas));
        histogram(C, 'EdgeColor', 'none', 'FaceColor', 'k', 'FaceAlpha', 0.3);

        % Formatting
        ylabel('count');
        xlabel(['phoneme at ' labels{t}]);
        set(gca, 'FontSize', 13);
        box off;

        subplot(2, 4, t*4);
        
        % Determine sound types (sonorant or obstruent) for the phonemes
        son = ismember(phntype{t}, corpus_details{ctr}.features.sonorant);
        obs = ismember(phntype{t}, corpus_details{ctr}.features.obstruent);
        C = categorical(son + 2.*obs, 0:2, {'na', 'sonorant', 'obstruent'});
        histogram(C, 'BarWidth', 0.4, 'FaceColor', 'k', ...
            'Normalization', 'Probability', 'EdgeColor', 'none', 'FaceAlpha', 0.3);
        
        % Formatting
        ylabel('probability');
        xlabel(['sound type at ' labels{t}]);
        xlim([categorical({'sonorant'}), categorical({'obstruent'})]);
        set(gca, 'FontSize', 13);
        box off;
    end 

    sgtitle(corpus{ctr});
    ctr = ctr + 1;
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;
%% Visualize probability of sonority around word onset

details = {timit_details, dimex_details};
titles = {'TIMIT', 'DIMEx'};
swrds = {TDwrd, Dwrd};
colors = getColorsCrossComp(3);

figure;
for s = 1:2
    Swrd = swrds{s};
    phnids = nan(height(Swrd), 11);
    sonority = nan(height(Swrd), 11);
    for i = 1:height(Swrd)
        precphns = Swrd.precPhns{i};
        if ~isempty(precphns)        
            phnids(i, 7-length(precphns):6) = precphns;
        end

        onsphns = Swrd.onsPhns{i};
        if ~isempty(onsphns)        
            phnids(i, 6:5+length(onsphns)) = onsphns;
        end

        % Convert phoneme ids to phoneme names to sonority values (0 / 1)
        for c = 1:size(phnids(i, :), 2)
            phnid = phnids(i, c);
            if ~isnan(phnid)
                % convert to phn
                charphn = details{s}.phnnames(phnid);

                % convert to sonority
                sonority(i, c) = ismember(charphn, ...
                    details{s}.features.sonorant);
            end
        end        
        clear onsphns precphns
    end    

    subplot(1, 2, s);

    % Plot probability of sonority at each position around word onset
    wrdsonor = sonority(logical(Swrd.wordOns), :);
    wrdprct = sum(wrdsonor, 'omitnan')./size(wrdsonor, 1);
    plot(smoothdata(wrdprct, 'SmoothingFactor', 0), ...
        'LineWidth', 2, 'Color', colors(2, :)); hold on;

    % Plot probability of sonority at each position around syllable onset
    syllsonor = sonority(logical(~Swrd.wordOns), :);
    syllprct = sum(syllsonor, 'omitnan')./size(syllsonor, 1);
    plot(smoothdata(syllprct, 'SmoothingFactor', 0), ...
        'LineWidth', 2, 'Color', colors(1, :)); hold on;

    % Formatting
    ylim([0.28 0.82]);
    xlim([1 11]);
    xlabel('# phns from word onset');
    ylabel('probability of sonorant');
    xline(6, 'LineWidth', 2, 'LineStyle','--');
    title(titles{s})
    box off;
    
    set(gca, 'FontSize', 15, 'Xtick', 1:5:11, ...
        'Xticklabels', split(num2str(-5:5:5)));
    legend({'Word', 'Syllable'});
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;


%% Visualize word length comparison in multisyllabic cases

colors = [0 0 1; 1 0 0];

figure;
subplot(1, 2, 1);
histogram(cellfun(@(x) length(x), TDwrd.phns(~isnan(TDwrd.syll))), 'Normalization', ...
    'probability', 'FaceColor', colors(2, :), 'EdgeColor','none'); hold on;

histogram(cellfun(@(x) length(x), Dwrd.phns(~isnan(Dwrd.syll))), 'Normalization', ...
    'probability', 'FaceColor', colors(1, :), 'EdgeColor','none');
legend({'TIMIT', 'DIMEx'});

% Formatting
ylabel('probability');
xlabel('# phns');
xlim([0 15]);
set(gca, 'FontSize', 15);
box off;

% Plotting result (length and syllable count)
subplot(1, 2, 2);
histogram(TDwrd.syll(~isnan(TDwrd.syll)), 'Normalization', ...
    'probability', 'FaceColor', colors(2, :), 'EdgeColor','none'); 
hold on;

histogram(Dwrd.syll(~isnan(Dwrd.syll)), 'Normalization', ...
    'probability', 'FaceColor', colors(1, :), 'EdgeColor','none'); 

% Formatting
ylabel('probability');
set(gca, 'FontSize', 15);
xlabel('# syllables');
xlim([1.5 6]);
box off

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% Visualize example sentence surprisal values

% Load in surprise values for TIMIT / DIMEx
[timit_details.sentdet]=makeSurprisal(timit_details.sentdet, 8,...
            'timit');
[dimex_details.sentdet]=makeSurprisal(dimex_details.sentdet, 8,...
            'dimex');

% Load the corpus annotation to IPA mappers
timit_conv = load("arpabet2ipa.mat");
dimex_conv = load('mexbet2ipa.mat');

% Select example sentences
sent = [80 42];
details = {timit_details, dimex_details, timit_details};
convs = {timit_conv, dimex_conv, timit_conv};
slabel = {'engSurp', 'spSurp', 'engAvgSurp'};
for l = 1:2 % TIMIT / DIMEx
    figure;
    sentinf = details{l}.sentdet(sent(l));

    % for conversion to ipa
    ipa = convs{l}.ipa;
    surp = sentinf.(slabel{l});
    if l == 1
        arpabet = convs{l}.arpabet;
    else
        arpabet = convs{l}.mexbet;
    end

    wrd = sentinf.wordOns;
    [r, tp] = find(sentinf.phnmatonset);
    phns  = arrayfun(@(x) details{l}.phnnames(x), r);
      
    ipas = cellfun(@(x) ipa(strcmp(arpabet, x)), phns, ...
        'UniformOutput',false);
    
    text(tp, 10*ones(1, length(phns)), ipas, 'FontSize',20); hold on;
    xlim([50 size(sentinf.aud, 2)-50]);
    
    % Plot surprisal values in a stem plot
    stem(find(surp>0), surp(surp>0), 'Color','k', 'LineWidth',2); hold on;
    xline(find(wrd>0)-1, 'LineWidth', 1.5, 'LineStyle', '--');
    xlim([50 size(sentinf.aud, 2)-50]);
    
    ylim([0 11]);
    title(strjoin(sentinf.wordList, ' '))
end


figure;
style = {'-', '--', ':'};
cols = {'r', 'b', 'm'};
allsurp=cell(1, 2);
for l = 1:2
    tmp = cellfun(@(x) x(x>0), {details{l}.sentdet.(slabel{l})}', ...
        'UniformOutput',false);
    spSurp = nan(length(tmp), max(cellfun(@(x) length(x), tmp)));
    for i = 1:length(tmp)
        spSurp(i, 1:length(tmp{i})) = tmp{i};
    end
    
    % Single sentence trial surprise (currently not plotted)
    [~,firstnan]=max(isnan(spSurp),[],2);
    [~, idx] = sort(firstnan);  
    subplot(1, 2, 1);
    imagesc(spSurp(idx, :));
    
    subplot(2, 1, 1);
    stem((1:15)+(l-2)/3, mean(spSurp(:, 1:15), 'omitnan'), 'Color', 'k', ...
        'LineWidth', 2, 'LineStyle', style{l}); hold on;
    
    wordSurp = [];
    for i = 1:length(tmp)
        tmp = details{l}.sentdet(i).(slabel{l});
        wordOns = [find(details{l}.sentdet(i).wordOns>0) ...
            length(tmp)];
        wordSurp = cat(2, wordSurp, arrayfun(@(x) tmp(wordOns(x):wordOns(x+1)-1), ...
            1:length(wordOns)-1, 'UniformOutput', false));
    end
    wordSurp = cellfun(@(x) x(x>0), wordSurp, 'UniformOutput', false);

    % Remove all empty lists
    wordSurp(cellfun(@(x) isempty(x), wordSurp)) = [];
    surp = nan(length(wordSurp), max(cellfun(@(x) length(x), wordSurp)));

    % Aggregate all word level surprisal values
    for i = 1:length(wordSurp)
        surp(i, 1:length(wordSurp{i})) = wordSurp{i};
    end
    
    subplot(2, 1, 2);
    stem((1:8)+(l-2)/3, rescale(mean(surp(:, 1:8), 'omitnan')), 'Color', cols{l}, ...
        'LineWidth', 2, 'LineStyle', style{l}); hold on;

    allsurp{l}=surp;
end


% Formatting
subplot(2, 1, 1);
caxis([0 10]);
xlim([0 16]);
ylabel('surprise');
xlabel('# phone');
set(gca, 'FontSize', 13, 'Ytick', 0:5);
box off;

subplot(2, 1, 2);
caxis([0 10]);
xlim([0 9]);
ylim([-0.05 1.05]);
ylabel('surprise');
xlabel('# phone');
set(gca, 'FontSize', 13, 'Ytick', 0:5);
box off;

labels = {'English', 'Spanish', 'cross-English'};
legend(labels);

% single word trial surprise
figure;
for l = 1:2

    [~,firstnan]=max(isnan(allsurp{l}),[],2);
    [~, idx] = sort(firstnan); 
    % suprrisal feature over all
    subplot(2, 2, l);
    imagesc(allsurp{l}(idx, :));
    caxis([0 10]);
    title(labels(l))

    % normalized surprisal feature
    subplot(2, 2, l+2);
    norms = allsurp{l}(idx, 1:8)./mean(surp(:, 1:8), 'omitnan');
    imagesc(norms);
    caxis([0 5])
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;


%% Visualize the envelope and waveform for several example syllables and words

sents = {[42, 27, 186], [333, 50, 88, 144]}; 
wrd = {{'cosas', 'mundo', 'gusta'}, {'touchdown', 'dollar', 'hopeful', 'every'}};
corpus_details = {dimex_details, timit_details};
for c = 1:2
    details = corpus_details{c};

    figure;
    for s = 1:length(sents{c})
        sent = sents{c}(s);
        subplot(length(sents{c}), 1, s);
        sent_info = details.sentdet(sent);

        % Find where the word starts
        widx = find(strcmp(sent_info.wordList, wrd{c}{s}));
        wtps = find(sent_info.wordOns);
        wtp = wtps(widx);

        % To make the word onset at time zero on the x-axis
        tosub = wtp/100; 

        x = ((1:length(sent_info.sound))/sent_info.soundf)-tosub;
        plot(x, sent_info.sound, 'Color', [0.5 0.5 0.5]); hold on;
        
        x = ((downsample(1:length(sent_info.sound), ...
            sent_info.soundf./sent_info.dataf))/sent_info.soundf)-tosub;
        offby = length(x)-length(sent_info.loudness);
        
        x = x(1+ceil(offby/2):end-floor(offby/2));

        % Plots the envelope over the waveform
        yyaxis right
        env = sent_info.loudness;               

        % Only plots syllables for this word
        sylons = find(sent_info.syltype>0);
        wrdons = find(sent_info.wordOns>0);
        cols = getColorsCrossComp(3);

        plot(x(wtp-15:wtp+10), env(wtp-15:wtp+10), ...
            'Color', cols(2, :), 'LineWidth', 2); hold on;
        plot(x(wtp+10:wtp+35), env(wtp+10:wtp+35), ...
            'Color', cols(1, :), 'LineWidth', 2, 'LineStyle', '-');

        xline(x(sylons(sylons>wtp))', 'Color', 'k', 'LineStyle', ...
            '--', 'LineWidth', 1.5);
        if ~isempty(wrdons(wrdons>wtp))
            xline(x(wrdons(wrdons>wtp))', 'LineWidth', 2);
        end

        % Formatting
        ylim([min(env(wtp-15:wtp+35)) max(env(wtp-15:wtp+35))]);
        yticks([])

        yyaxis left
        title([join(sent_info.wordList, ' ')]);
        xlim([x(wtp)-0.2 x(wtp)+0.5]);
        xline(x(wtp), 'k', 'LineWidth', 2.5) 
        set(gca, 'FontSize', 15, 'Xtick', [], 'Ytick', []);
        box off;
    end
    xlabel('Time (s)');
    set(gca, 'FontSize', 15, 'Xtick', -0.4:0.2:1, 'Ytick', []);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% Visualizing average pitch and phoneme surprise

befaft = [0.5 0.5];
x = -befaft(1):0.01:befaft(2);

allsent_tbl = Dwrd; 
pitchidx = cellfun(@(x) ~isempty(x), [allsent_tbl.pitch]);
wrd_idx = pitchidx & allsent_tbl.wordOns==1 & ~allsent_tbl.sentOns;
pitchsyl = cell2mat([allsent_tbl.pitch(pitchidx & allsent_tbl.wordOns==0)]);
pitchwrd = cell2mat([allsent_tbl.pitch(wrd_idx)]);

subplot(2, 2, 1);
leg = {'word', 'onset'};
if ~isempty(pitchsyl)
    shadedErrorBar(x, mean(pitchsyl, 1, "omitnan"), nansem(pitchsyl, 1), ...
        'lineprops', {'r'});
    leg = {'syl', 'word', 'onset'};
end
shadedErrorBar(x, mean(pitchwrd, 1, 'omitnan'), nansem(pitchwrd, 1), ...
    'lineprops', {'b'});
legend(leg);

% Formatting
xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
ylabel('norm pitch (z-score)'); 
xlabel('time from ons (s)')
set(gca, 'FontSize', 13);

% Visualizing preceding word length
subplot(2, 2, 2);
numbins = 3;
cols = brewermap(numbins, 'Dark2');
len_mean = nan(3, 1);
allplen = discretize(allsent_tbl.precWordLen, [0 3 5 15]); % [0 2 3 6]
for len = 1:numbins
    if sum(allplen == len & wrd_idx)>20
        y_len = cell2mat(allsent_tbl.pitch(allplen == len & wrd_idx));
        shadedErrorBar(x, mean(y_len, 1, 'omitnan'), nansem(y_len, 1), ...
            'lineprops', {'Color', cols(len, :)});
        len_mean(len, 1) = round(mean(allsent_tbl.precWordLen(allplen == ...
            len & wrd_idx)));
    end
end

% Formatting
% legend({'1 syll', '2 syll', '3+ syll'});
legend({['short:' num2str(len_mean(1))] , ...
    ['medium:' num2str(len_mean(2))], ['long:' num2str(len_mean(3))]});
xline(0, 'LineWidth', 2, 'LineStyle', '--', ...
    'Color','k', 'HandleVisibility', 'off');
set(gca, 'FontSize', 13);

subplot(2, 2, 4);
imagesc('Xdata', x, 'Cdata', pitchwrd);

% Formatting
ylim([0 size(pitchwrd, 1)]);
xlim([x(1) x(end)]);
xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
title('Word Onset');

ylabel('trials'); 
xlabel('time from ons (s)');
colormap([1 1 1; flipud(brewermap(8, 'RdBu'))]);
set(gca, 'FontSize', 13);

if ~isempty(pitchsyl)
    subplot(2, 2, 3);
    imagesc('Xdata', x, 'Cdata', pitchsyl);
    ylim([0 size(pitchsyl, 1)]);
    xlim([x(1) x(end)]);
    xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
    title('Syllable Onset')
    xlabel('Time from onset (s)');
    colormap([1 1 1; flipud(brewermap(8, 'RdBu'))]);
    ylabel('trials'); 
    xlabel('time from ons (s)')
    set(gca, 'FontSize', 13);
end

% visualizing result (surprise)
figure;
surpidx = cellfun(@(x) ~isempty(x), [allsent_tbl.surprise]);
surpsylresp = cell2mat([allsent_tbl.surprise(surpidx & allsent_tbl.wordOns==0)]);
shadedErrorBar(x, mean(surpsylresp, 1), nansem(surpsylresp, 1), 'lineprops', {'r'});
surpwrdresp = cell2mat([allsent_tbl.surprise(surpidx & allsent_tbl.wordOns==1)]);
shadedErrorBar(x, mean(surpwrdresp, 1), nansem(surpwrdresp, 1), 'lineprops', {'b'});

% Formatting 
ylabel('Surprise');
xlabel('Time (s)');
set(gca, 'FontSize', 13);

figure;
numbins = 3;
cols = brewermap(numbins, 'Dark2');
len_mean = nan(3, 1);
allplen = discretize(allsent_tbl.precWordLen, [0 3 6 20]); % [0 2 3 6]
for len = 1:numbins
    if sum(allplen == len & surpidx)>20
        y_len = cell2mat(allsent_tbl.surprise(allplen == len & surpidx));
        shadedErrorBar(x, mean(y_len, 1, 'omitnan'), nansem(y_len, 1), ...
            'lineprops', {'Color', cols(len, :)});
        len_mean(len, 1) = round(mean(allsent_tbl.precWordLen(allplen == ...
            len & surpidx)));
    end
end
legend({['short:' num2str(len_mean(1))] , ...
    ['medium:' num2str(len_mean(2))], ['long:' num2str(len_mean(3))]});

% Formatting and reference lines
set(gca, 'FontSize', 13);
ylabel('Surprise');
xlabel('Time (s)');
xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k', ...
    'HandleVisibility', 'off');
 
clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% ------------- FIX: ACOUSTIC CLASSIFICATION OF WORD BOUNDARY -----------------
%% Logistic on aud/envelope, 500ms preceding or surrounding onset

% to set
tps = 51-25:51+25; % 51-50:51
for i = {TDwrd, Dwrd}
    Swrd = i{1}; % dimex

    % remove all envelope instances that are naned
    envidx = cellfun(@(x) ~isempty(x), [Swrd.env]);
    addpath(genpath('shadederror/'))
    
    % using full window for 
    X_env = cell2mat([Swrd.env(envidx)]);
    X_env = X_env(:, tps);
    y = Swrd.wordOns(envidx)>0;
    x = (tps-bef)./100;

    [comp, score, ~, ~, exp] = pca(X_env);
    n_comp = find(diff(cumsum(exp)>80));
    [fp_env, tp_env, AUC_env, pcaX, scores, ~, ~, thresh] = ...
        logistic(X_env, y, 0, comp, tps, 10);

    % pca weights
    figure;
    subplot(1, 4, [1 2]);
    imagesc(comp(:, 1:n_comp)', 'XData', x)
    yticks([])
    ylabel('PCA components');
    xlabel('Time (s)');
    cbh = colorbar;
    cbh.Ticks = -0.2:0.1:0.2;
    colormap(flipud(brewermap(256, 'Spectral')));
    box off; set(gca, 'FontSize', 15);

    % ROC curve
    figure;
    subplot(1, 2, 1);
    plot(fp_env, tp_env, 'LineWidth', 2, 'Color', 'k'); hold on;
    h = refline(1, 0);
    h.LineStyle = '--';
    h.Color = 'k';
    xlabel('False positive rate') 
    ylabel('True positive rate')
    title(['AUC: ' num2str(mean(AUC_env))]);
    box off; set(gca, 'FontSize', 15);

    % FP/TP vs FN/TN examples
    subplot(1, 2, 2);
    tp = scores>=mean(thresh) & y;
    fp = scores>=mean(thresh) & ~y;
    tn = scores<mean(thresh) & ~y;
    fn = scores<mean(thresh) & y;
    labeltypes = [tp, fp, tn, fn];
    cols = flipud(brewermap(4, 'RdYlGn'));
    for l = 1:size(labeltypes, 2)
        y_l = X_env(labeltypes(:, l), :);
        shadedErrorBar(x, mean(y_l, 1), nansem(y_l, 1), ...
                'lineprops', {'Color', cols(l, :), 'linewidth', 2}); hold on;
    end
    legend({'TP', 'FP', 'TN', 'FN'});
    xline(0, 'LineStyle', '--', 'HandleVisibility','off');
    box off; set(gca, 'FontSize', 15);
    xlabel('Time (s)');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% (3 min) Logistic on aud/envelope for different times preceding onset

% timing = [5, 10, 20, 50]; 
% timing = [5, 10, 20];
timing = 25;
nreps = 20;
label = {'Envelope', 'Phn', 'Spectrogram', 'Spectrogram+Envelope'};
% timing = 50; 
wrdid = 1;
AUC_all = cell(1, 2);

for i = {Dwrd, TDwrd}
    Swrd = i{1};

    % randomly sample 2000 instances to test on
    rng(1)
    envidx = randsample(find(cellfun(@(x) ~isempty(x), [Swrd.env])), 2000);
    X_env = cell2mat([Swrd.env(envidx)]);
    X_aud = cat(3, Swrd.aud{envidx});

    y = Swrd.wordOns(envidx)>0;
    
    % colormap for dimex and timit    
    if startsWith(Swrd.sentid{1}, 's')
        cols = brewermap(5, 'YlGnBu'); % dimex color map
        cols = cols(2:end, :);
    else
        cols = brewermap(5, 'YlOrRd'); % timit color map
        cols = cols(2:end, :);
    end

    if length(timing)==1
        cols = cols(4:end, :);
    end
    
    ctr = 1;
    AUC = nan(4, length(timing), nreps);
    acc = nan(4, length(timing), nreps); % test accuracy
    weights = cell(length(timing), 1); % aud weights
    for t = timing
        disp(['timing: ' num2str(t)])
        % tps = 51-t:51; % before
        tps = 51-round(t/2):51+round(t/2); % around
    
        % 1 is env
        X = X_env(:, tps);
        [fp_env, tp_env, AUC(1, ctr, :), Xenv, ~, acc(1, ctr, :)] = ...
            logistic(X, y, 0, [], tps, nreps); % pca on env
    
        % 2 is env + phn (preceding and succeeding)
        X = [Swrd.precPhn(envidx) Swrd.onsPhn(envidx)]; % Xenv 
        phnidx = X(:, 1)>0 & X(:, 2)>0 & ~any(isnan(X), 2);
        [fp_envphn, tp_envphn, AUC(2, ctr, :), ~, scores_env, acc(2, ctr, :)] = ...
            logistic(X(phnidx, :), y(phnidx), 0, [], tps, nreps); % no pca on env + phn
    
        % 3 is aud
        X = reshape(X_aud(:, tps, :), 80*length(tps), []);
        [fp_aud, tp_aud, AUC(3, ctr, :), Xaud, ~, acc(3, ctr, :), weights{ctr}] ...
            = logistic(X', y, 1, [], tps, nreps);  
        
        % 4 is aud + env
        X = [Xenv Xaud];
        [fp_audenv, tp_audenv, AUC(4, ctr, :), pcaX, ~, acc(4, ctr, :)] = ...
                logistic(X, y, 0, [], tps, nreps);
    
        fps = {fp_env fp_envphn fp_aud fp_audenv};
        tps = {tp_env tp_envphn tp_aud tp_audenv};
        for j = 1:4
            subplot(2, 2, j);
            plot(fps{j}, tps{j}, 'LineWidth', 2, 'Color', cols(ctr, :)); hold on;
    
            % formatting
            h = refline(1, 0);
            h.LineStyle = '--';
            h.Color = 'k';
            xlabel('False positive rate') 
            ylabel('True positive rate')
            title([label{j} ' ROC']);
            legend(split(num2str(timing*10)));
            box off; set(gca, 'FontSize', 15);
        end
        ctr = ctr + 1;
    end
    AUC_all{wrdid}=AUC;
    wrdid=wrdid+1;
    
    % examine model weights
    figure; 
    subplot(1, 3, 1);
    imagesc(squeeze(mean(weights{max(timing)==timing}, 1))');
    ylabel('Mel spctrm');
    xlabel('Time (10ms)');
    colormap(flipud(spectral))
    set(gca, 'YDir', 'normal', 'FontSize', 15, 'Xtick', [1 25 50], ...
        'Xticklabel', {'-0.25', '0', '0.25'});
    
    avg = {cat(3, Swrd.aud{Swrd.wordOns & ~Swrd.sentOns}); ...
        cat(3, Swrd.aud{~Swrd.wordOns & ~Swrd.sentOns})};
    for j = 1:2
        subplot(1, 3, j+1);
        imagesc(squeeze(mean(avg{j}, 3)));
        ylabel('Mel spctrm');
        xlabel('Time (s)');
        colormap(flipud(spectral))
        set(gca, 'YDir', 'normal', 'FontSize', 15, 'Xtick', [1, 50, 100], ...
            'Xticklabel', {'-0.5', '0', '0.5'});
    end
end

tmp = brewermap(5, 'YlGnBu');
tmp2 = brewermap(5, 'YlOrRd');

% Plot all timings
cols = [{tmp(2:end, :)}, {tmp2(2:end, :)}];
figure;
feats = [1 3 4];
for i = 1:3 % which acoustic feature    
    feat = feats(i);
    subplot(1, 4, i);
    for t = 1:length(timing)
        for ls = 1:2
            boxplot(squeeze(AUC_all{ls}(feat, t, :))*100, 'Position', ...
                t-0.35+0.25*ls, 'Color', cols{ls}(t, :), ...
                'BoxStyle', 'filled', 'OutlierSize', 0.01); hold on;
        end
        p = ranksum(squeeze(AUC_all{1}(feat, t, :)), squeeze(AUC_all{2}(feat, t, :)));
        text(t, 83, getSigStr(p), 'FontSize', 15);
    end
    % formatting
    xticks(1:4)
    xticklabels(split(num2str(timing*10)));
    ylim([48 85]);
    yline(50, '--k', 'LineWidth', 2);
    
    if i == 1, ylabel('AUC'); end
    xlabel('Around onset (ms)');
    set(gca, 'FontSize', 13);
    title(label(feat));
        box off;

    subplot(1, 4, 4);
    for ls = 1:2
        boxplot(squeeze(AUC_all{ls}(2, 1, :))*100, 'Position', 0.65+0.25*ls, ...
            'Color', cols{ls}(length(cols{ls})-1, :), ...
            'BoxStyle', 'filled', 'OutlierSize', 0.01); hold on;
    end
    p = ranksum(squeeze(AUC_all{1}(2, 1, :)), squeeze(AUC_all{2}(2, 1, :)));
    text(1, 83, getSigStr(p), 'FontSize', 15);
    xticks([])
    ylim([50 85]);
    yline(0.5, 'Color', 'k');
    
    xlabel('Around onset (ms)');
    set(gca, 'FontSize', 13);
    title(label(2));
    yline(50, '--k', 'LineWidth', 2);
    box off;
end

cols = [{tmp(2:end, :)}, {tmp2(2:end, :)}];
figure;
feats = [1 2 3];
for i = 1:3 % which acoustic feature    
    feat = feats(i);
    auc = nan(2, 20);
    for ls = 1:2
        auc(ls, :) = squeeze(AUC_all{ls}(feat, size(AUC_all{ls}, 2), :));
        boxplot(auc(ls, :)*100, 'Position', i-0.35+0.25*ls, 'Color', ...
            cols{ls}(length(cols{ls})-1, :), 'BoxStyle', ...
            'filled', 'OutlierSize', 0.01); hold on;
    end
    p = ranksum(auc(1, :), auc(2, :));
    text(i, 83, getSigStr(p, 2), 'FontSize', 15);    
end

%formatting
ylim([48 85]);
yline(0.5, 'Color', 'k');
yline(50, '--k', 'LineWidth', 2);
xlabel('Feature');
ylabel('AUC');    
xticks(1:3);
xticklabels({'Envelope', 'Phoneme', 'Spectrogram'})
set(gca, 'FontSize', 13);
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* X_neural;

%% (2 min) Logistic on env/phnfeat/etc for different times preceding onset

% timing = [5, 10, 20, 50]; 
timing = 20;
nreps = 20;
fields = {'env', 'formant', 'consphnfeat', 'aud'};
pcaflags = [0 0 0 1];
Swrds = {Dwrd, TDwrd};
corpus_details = {dimex_details, timit_details};

AUC_all = cell(1, 2);
for lang = 1:2
    Swrd = Swrds{lang};

    % randomly sample 2000 instances to use
    rng(1);
    trlidx = randsample(find(cellfun(@(x) ~isempty(x), [Swrd.env])), 2000);
    
    % initialize variables
    AUC = nan(length(fields), length(timing), nreps);
    acc = nan(length(fields), length(timing), nreps); % test accuracy
    weights = cell(length(timing), 1); % aud weights

    ctr = 1;
    for t = timing
        for f = 1:length(fields)
            field = fields{f};
            disp(['Performing classification with  ' field '...'])
            [~, ~, AUC(f, ctr, :), ~, ~, acc(f, ctr, :), weights{ctr}, ~, ~, ~] = ...
                    acsword_logistic(Swrd, field, trlidx, nreps, t, ...
                    pcaflags(f), corpus_details{lang});
        end
        ctr=ctr+1;
    end
    AUC_all{lang} = AUC;
    clear AUC acc weights Mdl
end

tmp = brewermap(5, 'YlGnBu');
tmp2 = brewermap(5, 'YlOrRd');
cols = [{tmp(2:end, :)}, {tmp2(2:end, :)}];
figure;
feats = 1:length(fields);
for i = feats % which acoustic feature    
    auc = nan(2, 20);
    for ls = 1:2
        auc(ls, :) = squeeze(AUC_all{ls}(i, size(AUC_all{ls}, 2), :));
        boxplot(auc(ls, :)*100, 'Position', i-0.35+0.25*ls, 'Color', ...
            cols{ls}(length(cols{ls})-1, :), 'BoxStyle', ...
            'filled', 'OutlierSize', 0.01); hold on;
    end
    [h, p] = ttest2(auc(1, :), auc(2, :));
    text(i, 83, getSigStr(p, 2), 'FontSize', 15);    
end

%formatting
ylim([40 85]);
yline(0.5, 'Color', 'k');
yline(50, '--k', 'LineWidth', 2);
xlabel('Feature');
ylabel('AUC');    
xticks(1:length(fields));
xticklabels(fields)
set(gca, 'FontSize', 13);
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* X_neural;

%% (1 min) Logistic on aud cross-tested (English-Spanish, Spanish-English)
% timing = [5, 10, 20, 50]; 
% timing = [5, 10, 20];

timing = 20;
tps = 51-round(timing/2):51+round(timing/2); % around  
disp(['timing: ' num2str(timing)])

nreps = 20;
label = {'Spectrogram'};

% timing = 50;
AUC_all = cell(1, 2);
envidx = cell(1, 2);
X_aud = cell(1, 2);
y = cell(1, 2);

% first extract instances to train and test on for both corpora
ctr=1;
for i = {Dwrd, TDwrd}
    Swrd = i{1};

    % randomly sample 2000 instances to test on
    rng(1);
    envidx{ctr} = randsample(find(cellfun(@(x) ~isempty(x), [Swrd.env])), ...
        2000);
    X_aud{ctr} = cat(3, Swrd.aud{envidx{ctr}});
    y{ctr} = Swrd.wordOns(envidx{ctr})>0;
    ctr=ctr+1;
end

ctr = 1;
for i = {Dwrd, TDwrd}
    AUC = nan(2, length(timing), nreps);
    acc = nan(2, length(timing), nreps); % test accuracy
    weights = cell(length(timing), 1); % aud weights
          
    % 1 is aud
     X = reshape(X_aud{ctr}(:, tps, :), 80*length(tps), []);
    [fp_aud, tp_aud, AUC(1, ctr, :), Xaud, ~, acc(1, ctr, :), ...
        weights{ctr}, ~, Mdl, comp, mu] = logistic(X', y{ctr}, 1, [], tps, nreps); 

    % Generate the cross language predictions
    crossctr = mod(ctr, 2)+1;
    X_cross = reshape(X_aud{crossctr}(:, tps, :), 80*length(tps), []);
    [AUC(2, ctr, :), pcaX] = crosstestlogistic(X_cross', y{crossctr}', ...
        Mdl, comp, mu, nreps);
                         
    AUC_all{ctr}=AUC;
    ctr = ctr + 1;
    
    avg = {cat(3, Swrd.aud{Swrd.wordOns & ~Swrd.sentOns}); ...
        cat(3, Swrd.aud{~Swrd.wordOns & ~Swrd.sentOns})};
end
tmp = brewermap(5, 'YlGnBu');
tmp2 = brewermap(5, 'YlOrRd');

% plot the auc values
cols = [{tmp(2:end, :)}, {tmp2(2:end, :)}];
figure; 
for ls = 1:2
    auc = nan(2, 20);
    for feat = [1, 2] 
        auc(feat, :) = squeeze(AUC_all{ls}(feat, size(AUC_all{ls}, 2), :));
        boxplot(auc(feat, :)*100, 'Position', ls-0.35+0.25*feat, 'Color', ...
            cols{ls}(length(cols{ls})-1, :), 'BoxStyle', ...
            'filled', 'OutlierSize', 0.01); hold on;
        % feat-0.35+0.25*ls
    end
    [h, p] = ttest2(auc(1, :), auc(2, :));
    text(ls, 83, getSigStr(p, 2), 'FontSize', 15);    
end

% Formatting
ylim([48 85]);
yline(0.5, 'Color', 'k', 'HandleVisibility', 'off');
yline(50, '--k', 'LineWidth', 2, 'HandleVisibility', 'off');
xlabel('Feature');
ylabel('AUC');    
set(gca, 'FontSize', 13, 'Xtick', 1:2, 'Xticklabel', ...
    {'Same-Language', 'Cross-Language' });
box off;

figure;
cdata = [median(AUC_all{1}, [2, 3], 'omitnan'), ...
        median(AUC_all{2}, [2, 3], 'omitnan');]';
err = [nansem(AUC_all{1}(:, 1, :), 3), ...
    nansem(AUC_all{2}(:, 2, :), 3);];
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
ylim([0.5 1]);
yticks([0.5 1]);
box off; 
xlabel('Trained on');
legend({'Same-language', 'Cross-language'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* X_neural;


%% (3 min) Logistic on aud binned by difficulty compared to neural data

numbins = 5;
mintrials = 30;
cols = {flipud(purples), flipud(greens)};

for s = {Dwrd, TDwrd}
    Swrd = s{1};

    % decide on number bins
    % [N, ~] = histcounts(abs(scores_aud-Swrd.wordOns), 6);
    edges = linspace(0, 1, numbins);    
   
    % figure out corpus label
    if startsWith(Swrd.sentid{1}, 's')
        corpus = 'DIMEX';
        comb = [2, 4];
    else
        corpus = 'TIMIT';
        comb = [1, 4];
    end

    % show discretization
    figure;
    histogram(Swrd.ambiguity, 'EdgeColor', 'none');
    xline(edges);
    title(corpus);
    
    % ensure both classes have at least 30 examples per bin
    y_disc = discretize(Swrd.ambiguity, edges);
    bins = find(all(crosstab(Swrd.wordOns, y_disc)>mintrials));
    
    befaft = [0.5, 0.5];
    x = -befaft(1):0.01:befaft(2);
    label = {'Syllable', 'Word'};
    debug=0;
    if debug == 1
        for i = 0:1 % syllable or word
            figure('Position', [1, 2, 800, 100]);
            for b = bins    
                subplot(1, numbins, b);
            
                % plot envelope for each bin
                ysyll =  cell2mat(Swrd.env(y_disc==b & Swrd.wordOns==i));
                cdata = normr(ysyll(1:numint, :));
                imagesc('Xdata', x, 'Cdata', smoothdata(cdata, 1)); hold on;
                colormap(cols{i+1})
                brighten(0.3);
                yticks([]);
            
                set(gca, 'FontSize', 13);
                xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
                ylim([0 numint]);
                xlim([x(1) x(end)]);
                xline(0, 'LineWidth', 1.75, 'LineStyle', '--')
                yyaxis right
                plot(x, smoothdata(mean(normr(ysyll)), 'SmoothingFactor', 0.15), ...
                    'Color', [0.5 0.5 0.5], 'linewidth', 3);    
                yticks([]);
            end
            sgtitle(corpus)
        
            figure('Position', [1, 2, 800, 100]);
            for b = bins               
                subplot(1, length(N), b);
        
                ysyll =  cat(3, Swrd.aud{y_disc==b & Swrd.wordOns==i});
                imagesc('Xdata', x, 'Cdata', squeeze(mean(ysyll, 3)));
                xline(0, 'LineWidth', 1.75, 'LineStyle', '--')
                set(gca, 'YDir', 'normal');
                colormap(inferno);
                ylim([0 80]);
                yticks([])
                
            end
            sgtitle([corpus ' ' label{i+1}])
        end
    end
    
    % split neural single trial data into bins
    time_label = '500ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
    % discretize neural data in the same way as acoustic to see AUC per bin
    AUC = nan(4, numbins);
    diffscore = nan(4, numbins, 2);
    % ls x bins x (mean, sem)
    for l = [1, 2] % monolinguals
        for bin = 1:numbins
            scores = decode_details.score{l}(y_disc==bin);
            y = decode_details.y(y_disc==bin);
            nanidx = isnan(scores);
            if sum(~nanidx)>50
                % difference in auc
                [~, ~, t, auc] = perfcurve(y(~nanidx), scores(~nanidx), 1);
                AUC(l, bin) = auc;
    
                % difference in scores
                tmp = abs(scores(~nanidx)-y(~nanidx));
                diffscore(l, bin, 1) = mean(tmp);
                diffscore(l, bin, 2) = nansem(tmp);
            end
        end
    end
       
    ls = 4;
    subj = 'bilingual';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
    % ls x bins x (mean, sem)
    for l = ls
        disp(l) 
        for bin = 1:length(N)
            scores = decode_details.score{l}(y_disc==bin);
            y = decode_details.y(y_disc==bin);
            nanidx = isnan(scores);
            if sum(~nanidx)>50
                % difference in auc
                [~, ~, t, auc] = perfcurve(y(~nanidx), scores(~nanidx), 1);
                AUC(l, bin) = auc;
    
                % difference in scores
                tmp = abs(scores(~nanidx)-y(~nanidx));
                diffscore(l, bin, 1) = median(tmp);
                diffscore(l, bin, 2) = nansem(tmp);
            end
        end
    end

    % COLLAPSE: PLOT ALL BINS
    % monolinguals +
    % bilingual vs. nonnative monolingual
    lang_combs = {[1, 2], [1, 4], [2, 4]};
    ctr = 1;
    figure('Renderer', 'Painters'); 
    for lc = lang_combs
        lang_comb = lc{1};       
        cols = getColorsCrossComp(1);
        for l = lang_comb
            subplot(1, 3, ctr)  
            errorbar(1:numbins, diffscore(l, :, 1), diffscore(l, :, 2), ...
                'color',cols(l, :), 'linestyle', '-', 'lineWidth', 1);hold on
            scatter(1:numbins, diffscore(l, :, 1), 25, cols(l, :), 'filled');
    
            %AUC doesn't look as clean
            yyaxis right
            ylabel('AUC')
            yticks([0 0.9]);
            plot(1:numbins, AUC(l, :), 'color', ...
                cols(l, :), 'linestyle', '-'); hold on;
            scatter(1:numbins, AUC(l, :), 25, 'filled'); hold on;
            yyaxis left
        end
        ylabel('abs(score - correct label)');
        yticks([0 0.9]);
        box off;
        set(gca,  'FontSize', 13, 'YDir', 'reverse');
        xlabel('Binned by acoustic difficulty');
        xlim([0.75 length(N)+0.25]);
        ylim([0.2 0.61]);

%         yyaxis right     
%         meandiff = abs(diffscore(lang_comb(1), :, 1)-diffscore(lang_comb(2), :, 1));
%         plot(x, meandiff,'color', 'k', 'linestyle', ':', 'LineWidth', 2); 
%         ylabel('average group difference');
%         yticks([]);
%         ylim([0 0.1]);
        ctr = ctr+1;
    end

    % PLOT EXTREME BINS
    lang_combs = {[1, 2], [1, 4], [2, 4]};
    ctr = 1;
    figure('Renderer', 'Painters'); 
    for lc = lang_combs
        lang_comb = lc{1};       
        cols = getColorsCrossComp(1);
        for l = lang_comb
            subplot(1, 3, ctr)
            x = 1:2;   

            % mean and error
            lowambig = squeeze(mean(diffscore(l, 1:2, :), 2));
            % mean and error
            highambig = squeeze(mean(diffscore(l, end-2:end-1, :), 2));

            errorbar(x, [lowambig(1) highambig(1)], [lowambig(2) highambig(2)], ...
                'color',cols(l, :), 'linestyle', '-', 'lineWidth', 1);hold on
            
        end
        ylabel('diff(score, label)')
        yticks([0 0.9]);
        box off;
        set(gca,  'FontSize', 13, 'YDir', 'reverse');
        xlabel('acoustic ambiguity');
        xlim([0.75 2.25]);
        ylim([0.2 0.61]);
    
%         yyaxis right     
%         meandiff = abs(diffscore(lang_comb(1), ...
%             [1, end-1], 1)-diffscore(lang_comb(2), [1, end-1], 1));
%         plot(x, meandiff,'color', 'k', 'linestyle', ':', 'LineWidth', 2); 
%         ylabel('average group difference');
%         yticks([]);
%         ylim([0 0.1]);
        ctr = ctr+1;
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* diffscore N;

%% (3 min) WORD: Logistic on aud binned by difficulty compared to neural data

% corpi = {'DIMEX', 'TIMIT'};
% Swrds = {Dwrd, TDwrd}; % TDwrd

corpi = {'TIMIT'};
Swrds = {TDwrd}; % TDwrd
acsbins = 2;
diffscore = cell(2, 4, acsbins);
sorted_neur = cell(2, 4);

for lang = 1:length(corpi)
    Swrd = Swrds{lang};
    corpus = corpi{lang};

    % Perform spectrogram classification
    tps = 51-20:51+20; % 400 ms around onset
    X_aud = cat(3, Swrd.aud{:});
    X = reshape(X_aud(:, tps, :), 80*length(tps), []);
    [~, ~, ~, ~, scores_aud, ~, ~] = logistic(X', Swrd.wordOns, 1, [], tps, 10);

    % ambiuity index is the difference between the predicted score and the
    % true label
    AAI = abs(scores_aud-Swrd.wordOns);
    
    % Decide on number of bins
    [N, ~] = histcounts(AAI, acsbins);
    edges = linspace(prctile(AAI, 5), prctile(AAI, 95), acsbins+1);    
   
    % Figure out corpus label
    if startsWith(Swrd.sentid{1}, 's')
        comb = [2, 4];       
    else
        comb = [1, 4];
    end
    
    % Ensure both classes have at least 30 examples per bin
    [~, sorted] = sort(AAI);
    y_disc = discretize(AAI, edges);
    numint = 30;
    bins = find(all(crosstab(Swrd.wordOns, y_disc)>30));

    % show high ambiguity and low ambiguity examples for English
    if strcmpi(corpus, 'TIMIT')

        figure;
        tps = -0.5:0.01:0.5;

        % low ambiguity example (spectrogram)
        subplot(2, 2, 1);

        % ensure it is a word onset
        lowambig = find(Swrd.wordOns(sorted), 1, 'first');
        imagesc(TDwrd(sorted(lowambig), :).aud{1}, 'Xdata', tps);
        set(gca, 'Ydir', 'normal', 'FontSize', 13);
        yticks([1, 80]); 
        colormap(flipud(gray));
        yticklabels({'0', '8'}); 
        ylabel('kHz');
        h = xline(0);
        h.Color = 'r';
        xlabel('Time (s)');
        title([join(num2phns(TDwrd(sorted(lowambig), :).precPhns{:}, timit_details), ' ') ...
            '|' join(TDwrd(sorted(lowambig), :).phns{:}, ' ') ]);
        

        % (envelope)
        subplot(2, 2, 2);
        plot(tps, TDwrd(sorted(lowambig), :).env{1}, ...
            'LineWidth', 2, 'Color', [0.7 0.2 0.8]);
        set(gca, 'Ydir', 'normal', 'FontSize', 13);
        ylim([-1 4.5]);
        h = xline(0);
        h.Color = 'k';
        title(['AAI = ' num2str(AAI(sorted(lowambig)))]);

        % high ambiguity example (spectrogram)
        subplot(2, 2, 3);
        highambig = find(Swrd.wordOns(sorted), 10, 'last');
        highambig = highambig(5);
        imagesc(TDwrd(sorted(highambig), :).aud{1}, 'Xdata', tps);% end-7
        set(gca, 'Ydir', 'normal', 'FontSize', 13);
        yticks([1, 80]);
        colormap(flipud(gray));
        yticklabels({'0', '8'}); 
        ylabel('kHz');
        h = xline(0);
        h.Color = 'r';
        xlabel('Time (s)');
        title([join(num2phns(TDwrd(sorted(highambig), :).precPhns{:}, timit_details), ' ') ...
            '|' join(TDwrd(sorted(highambig), :).phns{:}, ' ')])

        % (envelope)
        subplot(2, 2, 4);
        plot(tps, TDwrd(sorted(highambig), :).env{1}, ...
            'LineWidth', 2, 'Color', [0.7 0.2 0.8]);
        set(gca, 'Ydir', 'normal', 'FontSize', 13);
        ylim([-1 4.5]);
        h = xline(0);
        h.Color = 'k';
        xlabel('Time (s)');
        title(['AAI = ' num2str(AAI(sorted(highambig)))]);
    end

    % use the same number of trials per bin for neural comparison
    % min_trials = min(crosstab(y_disc)); 
    
    befaft = [0.5, 0.5];
    x = -befaft(1):0.01:befaft(2);
    label ={'Syllable', 'Word'};
    
    % ls = [1 2];
    time_label = '600ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
    
    % Discretize neural data in the same way as acoustic to see AUC per bin
    AUC = nan(4, length(N));
    % ls x bins x (mean, sem)

    for l = [1, 2, 4] % Monolinguals
        if l == 4
            subj = 'bilingual';
            filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
            load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
        end     

        trlnum = nan(length(N), 1);
        for bin = 1:length(N)
            scores = decode_details.score{l}(y_disc==bin);
            trlnum(bin) = sum(~isnan(scores));
        end
        min_trials = min(trlnum);
        disp([num2str(l) ', min trials= ' num2str(min_trials)]);
        
        % Show continuous measure of scores
        scores = decode_details.score{l};
        nanidx = isnan(scores);
        sorted(ismember(sorted, find(nanidx)))=[];
        sorted_neur{lang, l} = abs(Swrd.wordOns(sorted)-scores(sorted));

        % Look at binned scores
        for bin = 1:length(N)
            scores = decode_details.score{l}(y_disc==bin);
            y = decode_details.y(y_disc==bin);
            nanidx = isnan(scores);
            c = find(strcmp(corpi, corpus));
            if sum(~nanidx)>50
                idx = find(~nanidx);

                % Take a random sample of trials to calculate the score difference
                rng(1)
                idx = idx(randperm(length(idx), min_trials));
                [~, ~, t, auc] = perfcurve(y(idx), scores(idx), 1);
                AUC(l, bin) = auc;
    
                % difference in scores
                tmp = abs(y(idx)-scores(idx));
                diffscore{c, l, bin} = tmp;
            end
        end
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* N diffscore;

%%  Plotting figure with ambiguity for the paper

figure;
Swrds = {Dwrd, TDwrd};
% red for unfamiliar language, blue for native
colors = [0.8 0.2 0.2;0.2 0.2 0.8];
corpi = {'DIMEX', 'TIMIT'};

for cr = 1:2
    corpus = corpi{cr};
    Swrd = Swrds{strcmpi(corpus, 'TIMIT')+1};
    time_label = '600ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

    binedges = [0 0.05 0.7 1];
    ambiguity = discretize(Swrd.ambiguity, binedges);

    subplot(1, 2, strcmpi(corpus, 'TIMIT')+1);
    for bin = [1, 3]

        diffscore_stats = cell(2, 1);
        for ls = [1, 2]
            scores = decode_details.score{ls};
            y = decode_details.y;
            diffscore = abs(y-scores);

            % box plot of score differences in this bin
            % one subplot per corpus

            % get the scores for this bin
            diffscore_bin = diffscore(ambiguity==bin);
            nanidx = isnan(diffscore_bin);
            diffscore_bin(nanidx) = [];

            % save out scores for stats
            diffscore_stats{ls} = diffscore_bin;

            % display number of trials
            disp([corpus ', ' num2str(ls) ', ' num2str(bin) ', ' ...
                num2str(sum(~nanidx))]);

            % make error bars with nansem as error
            err(1) = mean(diffscore_bin);
            err(2) = nansem(diffscore_bin);
            x = bin-0.35+0.25*double(ls~=cr);

            % colors by native (ls==cr) or nonnative
            errorbar(x, err(1), err(2), 'Color', ...
                colors(double(ls==cr)+1, :), 'LineWidth', 2.5); hold on;

            % flip the axis
            set(gca, 'YDir', 'reverse', 'FontSize', 13);
            ylabel('abs(score - label)');
            ylim([0.25 0.55]);
            yticks([0.25 0.5]);
            box off;
            xticks([1 3]);
            xlim([0 4]);
            xticklabels({'low', 'high'});
            xlabel('Acoustic ambiguity');
        end
        % stats
        [h, p] = ttest2(diffscore_stats{1}, diffscore_stats{2});
        text(bin-0.25, 0.3, getSigStr(p, 2), 'FontSize', 15);

        % legend
        if bin == 1 && strcmpi(corpus, 'DIMEX')
            h = legend({'native', 'foreign'}, 'Location', 'NorthEast');
        end
    end
end

%% Example of high ambiguity and low ambiguity word boundaries

highambig = 'msdh0_si2240'; % try another
lowambig = 'mrew1_si2130'; % this bother

% show spectrogram, and envelope and text


%% (3 min) ONSET SONORITY: Logistic on aud binned by difficulty compared to neural data

corpi = {'DIMEX', 'TIMIT'};
Swrds = {Dwrd, TDwrd}; % TDwrd
details = {dimex_details, timit_details};
acsbins = 3;
diffscore = cell(2, 4, acsbins);
sorted_neur = cell(2, 4);

for lang = 1:2
    Swrd = Swrds{lang};
    corpus = corpi{lang};
    corpus_details = details{lang};

    % perform spectrogram classification
    tps = 51-20:51+20; % around
    X_aud = cat(3, Swrd.aud{:});
    X = reshape(X_aud(:, tps, :), 80*length(tps), []);

    % onset sonority
    onsPhns = arrayfun(@(x) corpus_details.phnnames(x), [Swrd.onsPhn]);
    onsSonor = cellfun(@(x) ismember(x, corpus_details.features.sonorant), onsPhns);
    [~, ~, ~, ~, scores_aud, ~, ~] = logistic(X', onsSonor, 1, [], tps, 10); 
    [~, sorted] = sort(abs(scores_aud-onsSonor));
    
    % decide on number bins
    [N, ~] = histcounts(abs(scores_aud-onsSonor), acsbins);
    edges = linspace(prctile(abs(scores_aud-onsSonor), 10), ...
        prctile(abs(scores_aud-onsSonor), 90), acsbins+1);    
   
    % figure out corpus label
    if startsWith(Swrd.sentid{1}, 's')
        comb = [2, 4];
    else
        comb = [1, 4];
    end
    
    % ensure both classes have at least 30 examples per bin
    y_disc = discretize(abs(scores_aud-onsSonor), edges);
    numint = 30;
    bins = find(all(crosstab(onsSonor, y_disc)>30));
    % use the same number of trials per bin for neural comparison
    % min_trials = min(crosstab(y_disc)); 
    
    befaft = [0.5, 0.5];
    x = -befaft(1):0.01:befaft(2);
    label ={'Vowel', 'Consonant'};
    
    % ls = [1 2];
    time_label = '600ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_onsSonority_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

    % discretize neural data in the same way as acoustic to see AUC per bin
    AUC = nan(4, length(N));
    % ls x bins x (mean, sem)
    for l = [1, 2, 4] % monolinguals + bilinguals
        if l == 4
            subj = 'bilingual';
            filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
            load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');
        end

        trlnum = nan(length(N), 1);
        for bin = 1:length(N)
            scores = decode_details.score{l}(y_disc==bin);
            trlnum(bin) = sum(~isnan(scores));
        end
        min_trials = min(trlnum);
        disp([num2str(l) ', min trials= ' num2str(min_trials)]);

        % show continuous measure of scores
        scores = decode_details.score{l};
        nanidx = isnan(scores);
%         sorted(ismember(sorted, find(nanidx)))=[];
        sorted_neur{lang, l} = abs(onsSonor(sorted)-scores(sorted));
        
        c = find(strcmp(corpi, corpus));
        for bin = 1:length(N)
            scores = decode_details.score{l}(y_disc==bin);
            y = decode_details.y(y_disc==bin);
            nanidx = isnan(scores);
            if sum(~nanidx)>50
                idx = find(~nanidx);         
                
                % rand sample
                rng(1)
                idx = idx(randperm(length(idx), min_trials));
                [~, ~, t, auc] = perfcurve(y(idx), scores(idx), 1);
                AUC(l, bin) = auc;
    
                % difference in scores
                tmp = abs(y(idx)-scores(idx));
                diffscore{c, l, bin} = tmp;
            end
        end
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% Visualize the binning method and average envelope per bin
corpi = { 'TIMIT'}; % 'DIMEX',
acsbins = 3;
diffscore = cell(2, 4, acsbins);
for s = {TDwrd} % Dwrd, 
    Swrd = s{1};

    % perform spectrogram classification
    tps = 51-20:51+20; % around
    X_aud = cat(3, Swrd.aud{:});
    X = reshape(X_aud(:, tps, :), 80*length(tps), []);
    [~, ~, ~, ~, scores_aud, ~, ~] = logistic(X', Swrd.wordOns, 1, [], tps, 10);  
    
    % decide on number bins
    [N, ~] = histcounts(abs(scores_aud-Swrd.wordOns), acsbins);
    edges = linspace(0, 1, acsbins+1);    
   
    % figure out corpus label
    if startsWith(Swrd.sentid{1}, 's')
        corpus = 'DIMEX';
        comb = [2, 4];
    else
        corpus = 'TIMIT';
        comb = [1, 4];
    end

    % show discretization
    figure;
    histogram(abs(scores_aud-Swrd.wordOns), 'EdgeColor', 'none', ...
        'FaceColor', [0.6 0.6 0.6]);
    xline(edges);
    title(corpus)
    
    % ensure both classes have at least 30 examples per bin
    y_disc = discretize(abs(scores_aud-Swrd.wordOns), edges);
    numint = 100;
    bins = find(all(crosstab(Swrd.wordOns, y_disc)>30));
    
    befaft = [0.5, 0.5];
    x = -befaft(1):0.01:befaft(2);
    label ={'Syllable', 'Word'};
    for i = 0:1
        figure('Position', [1, 2, 800, 100]);
        for b = bins    
            subplot(1, length(N), b);
        
            ysyll =  cell2mat(Swrd.env(y_disc==b & Swrd.wordOns==i));
            cdata = normr(ysyll(1:numint, :));
            imagesc('Xdata', x, 'Cdata', smoothdata(cdata, 1), ...
                'AlphaData', 0.5); hold on;
            
            if i == 0
                colormap(flipud(gray));
            else
                colormap(flipud(fpurple));
            end
            brighten(0.2);
            yticks([]);
        
            set(gca, 'FontSize', 13);
            xline(0, 'LineWidth', 2, 'LineStyle', '--', 'Color','k');
            ylim([0 numint]);
            xlim([x(1) x(end)]);
            xline(0, 'LineWidth', 1.75, 'LineStyle', '--')
            yyaxis right
            if i == 0
                plot(x, smoothdata(mean(normr(ysyll)), 'SmoothingFactor', 0.15), ...
                    'Color', [0 0 0], 'linewidth', 3);    
            else
                plot(x, smoothdata(mean(normr(ysyll)), 'SmoothingFactor', 0.15), ...
                    'Color', [0.7 0.2 0.4], 'linewidth', 3);   
            end
            ylim([0 0.15])
            yticks([]);
        end
        sgtitle(corpus)
       
        figure('Position', [1, 2, 800, 100]);
        for b = bins               
            subplot(1, length(N), b);
    
            ysyll =  cat(3, Swrd.aud{y_disc==b & Swrd.wordOns==i});
            imagesc('Xdata', x, 'Cdata', squeeze(mean(ysyll, 3)));
            xline(0, 'LineWidth', 1.75, 'LineStyle', '--')
            set(gca, 'YDir', 'normal');
            colormap(inferno);
            ylim([0 80]);
            yticks([])
            
        end
        sgtitle([corpus ' ' label{i+1}])
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% vis: plotting

% monolinguals +
% bilingual vs. nonnative monolingual
lang_combs = {[1, 2], [1, 4], [2, 4]};
cols = getColorsCrossComp(1);
corpi = {'DIMEX', 'TIMIT'};

figure('Renderer', 'Painters'); 
ctr = 1;
for corp = 1:2
    for lc = lang_combs
        lang_comb = lc{1};       
        
        for l = lang_comb
            subplot(2, 3, ctr)
            x = 1:length(N);   
            y = squeeze(cellfun(@(x) mean(x, 'omitnan'), diffscore(corp, l, :)));
            err = squeeze(cellfun(@(x) nansem(x), diffscore(corp, l, :)));
            errorbar(x, y, err, 'k', 'linestyle', '-', 'lineWidth', 1);
            hold on

            scatter(x, y, 25, cols(l, :), 'filled');
        end

        % see whether difference is significant
        for i = 1:length(squeeze(diffscore(corp, l, :)))
            x = squeeze(diffscore{corp, lang_comb(1), i});
            y = squeeze(diffscore{corp, lang_comb(2), i});
            [h, p]=ttest2(x,y);
            text(i, 0.3, getSigStr(p, 2));
        end

        ylabel('abs(score - correct label)')
        
        set(gca,  'FontSize', 13, 'YDir', 'reverse');
        xlabel('acoustic ambiguity');
        xlim([0.75 length(N)+0.25]);
        ylim([0.25 0.7]);

        sgtitle(corpi{corp});
        ctr=ctr+1;
    end
end

% PLOT EXTREME BINS
native = nan(2, size(diffscore, 3));
nonnative = nan(2, size(diffscore, 3));
biling = nan(2, size(diffscore, 3));
for corp = 1:2
    for ls = [1, 2]
        % figure out which one is native vs. nonnative
        if ls == corp % dependent on ORDERING!
            % first two bins
            native(corp, :) = squeeze(cellfun(@(x) mean(x, 'omitnan'), ...
                diffscore(corp, ls, :)));
        elseif ismember(ls, [1, 2])
            % first two bins
            nonnative(corp, :) = squeeze(cellfun(@(x) mean(x, 'omitnan'), ...
                diffscore(corp, ls, :)));
        else % bilingual
            biling(corp, :) = squeeze(cellfun(@(x) mean(x, 'omitnan'), ...
                diffscore(corp, ls, :))); 
        end     
    end
end

figure; 
bar(-1*mean(native-nonnative));
ylim([-0.02 0.08])

% figure;
% diffs = {abs(nonnative-native), abs(biling-native)};
% titles = {'Monolingual: Native vs. Non-Native', 'Native: Monolingual vs. Bilingual'};
% for comp = 1:2
%     subplot(1, 2, comp);
%     half = floor(size(diffs{comp}, 2)/2);
%     boxchart(ones(half*2, 1), reshape(diffs{comp}(:, 1:half), half*2, 1), ...
%         'BoxFaceColor', [0.5 0.5 0.5]); hold on; 
%     boxchart(ones(half*2, 1)*2, reshape(diffs{comp}(:, half+1:end), half*2, 1), ...
%         'BoxFaceColor', [0.5 0.5 0.5] ); 
%     xlabel('Acoustic difficulty');
%     xticks([1, 2]);
%     xticklabels({'min', 'max'});
%     ylabel('AUC difference');
%     yticks(0:0.5:1);
%     ylim([0 0.1]);
%     set(gca, 'FontSize', 15);
% 
%     [p, h] = ranksum(reshape(diffs{comp}(:, 1:2), 4, 1), ...
%         reshape(diffs{comp}(:, 3:4), 4, 1));
%     line([1, 2], [0.08 0.08], 'Color', 'k', 'LineWidth', 1.5)
%     text(1.35, 0.09, getSigStr(p, 2), 'FontSize', 15)
%     title(titles{comp})
% end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* N;

%% vis: single scores (not binned)

titles = {'Spanish speech', 'English speech'};
cols = getColorsCrossComp(1);
figure;
for lang = 1:2

    subplot(1, 4, (lang-1)*2+1)
    y_diff = nan(2, 3);
    for l = [1, 2]

        % number of points to average
        b = 10;

        % only include trials that are included in the bins
        tmp = cat(1, diffscore{lang, l, :});
        binned = arrayfun(@(x) mean(tmp(x:x+b)), 1:b+1:length(tmp)-b);
        x = rescale(1:length(binned), 1, 3);
        y = binned;
        scatter(x, y,  15, cols(l, :), 'filled', 'MarkerFaceAlpha', 0.5); hold on;  

        x = 1:length(N);   
        y_diff(l, :) = squeeze(cellfun(@(x) mean(x, 'omitnan'), diffscore(lang, l, :)))';
        err = squeeze(cellfun(@(x) nansem(x), diffscore(lang, l, :)));
        errorbar(x, y_diff(l, :), err, 'Color', cols(l, :), 'linestyle', '-', ...
            'lineWidth', 1.8);hold on

        ylim([0.25 0.7]);
        yticks(0.3:0.1:0.6)
        set(gca, 'YDir', 'reverse');
         
%         % line fitting
%         degree = 2;
%         [p,S] = polyfit(x,y,degree);
% 
%         % get confidence bounds
%         alpha = 0.80; % Significance level
%         [yfit,delta] = polyconf(p,x,S,'alpha', alpha);
% 
%         % plot lines
%         plot(x,yfit, 'LineWidth', 2, 'Color', cols(l, :));
%         set(gca, 'YDir', 'reverse')
%         patch([x, fliplr(x)],[yfit-delta, fliplr(yfit+delta)], ...
%             cols(l, :), 'FaceAlpha', 0.1, 'EdgeColor', 'none');
    end

    for i = 1:length(squeeze(diffscore(lang, 1, :)))
        x = squeeze(diffscore{lang, 1, i});
        y = squeeze(diffscore{lang, 2, i});
        [~, p]=ttest2(x,y);
        text(i, 0.3, getSigStr(p, 2));
    end
    h=title(titles(lang));
    h.FontWeight='Normal';
    set(gca, 'FontSize', 15);

    subplot(1, 4, (lang-1)*2+2)
    bar(abs(diff(y_diff)), 'FaceColor', [0.5 0.5 0.5], ...
        'EdgeColor', 'none');
    box off;
    ylim([0 0.1]);
    yticks([0 0.1]);
    set(gca, 'FontSize', 15);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* diffscore N;


%% ------------------ funtcions -------------------------------------------

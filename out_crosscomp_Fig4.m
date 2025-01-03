%% Set up
% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;
SIDs = [sSIDs, eSIDs, {'HS11', 'HS9', 'HS10'}];

% change this for the word onset analysis )
bef=50;
aft=50;

% eSIDs sSIDs
Dwrd = loadDwrdMulti('dimex',  bef, aft, [], dimex_details); % {'EC260', 'EC266'}
TDwrd = loadDwrdMulti('timit', bef, aft, [], timit_details); % {'EC260', 'EC266'}

window = 20;
fieldname = 'aud';
Dwrd.ambiguity = getAAI(window, Dwrd, fieldname, 20);
TDwrd.ambiguity = getAAI(window, TDwrd, fieldname, 20);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;
%% A - Histogram of AAI and example low / high AAI instances

titles = {'English', 'Spanish'};
figure;
subplot(2, 1, 2);
histogram(TDwrd.ambiguity, 'EdgeColor', 'none', 'FaceColor', ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'Normalization', 'probability'); hold on;
ylabel('Probability');
yticks([0.0 0.2 0.4]);
xticks(0:0.5:1)
xlabel('Acoustic Ambiguity Index (AAI)');
xline([0.949, 0.0322]) % AAI of examples
yl = ylim();
text([0.949, 0.0322], [yl(2), yl(2)], {'high AAI example', 'low AAI example'}, ...
    'FontSize', 13, 'HorizontalAlignment', 'center');
set(gca, 'FontSize', 13);
title(titles{1}, 'FontSize', 15, 'FontWeight', 'normal'); 
box off;

subplot(2, 1, 1);
histogram(Dwrd.ambiguity, 'EdgeColor', 'none', 'FaceColor', ...
    [0.5 0.5 0.5], 'FaceAlpha', 0.5, 'Normalization', 'probability'); hold on;
ylabel('Probability');
yticks([0.0 0.2 0.4]);
xticks([]);
set(gca, 'FontSize', 13);
xline([0.9883, 0.01214]) % AAI of examples
% 'high AAI example' 0.9883 and 'low AAI example' 0.01214
yl = ylim();
text([0.9883, 0.01214], [yl(2), yl(2)], {'high AAI example', 'low AAI example'}, ...
    'FontSize', 13, 'HorizontalAlignment', 'center');
title(titles{2}, 'FontSize', 15, 'FontWeight', 'normal'); 
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% A - High and low AAI examples in Spanish and English 

% highambig = 'msdh0_si2240'; % try another
% lowambig = 'mrew1_si2130'; % this bother
% % wrdidx = find(strcmp(Dwrd.sentid, highambig));
wrdidx = [1832, 1676];
maxnum = 7;
Swrd = TDwrd;

% % find 5 highest and 5 lowest AAI that are also word onsets
[~, idx] = sort(Swrd.ambiguity, 'descend');
wrdidx = idx(Swrd.wordOns(idx)>0);
wrdidx = wrdidx(1:maxnum);
% 
% find 5 lowest AAI that are also word onsets
[~, idx] = sort(Swrd.ambiguity, 'ascend');
wrdidx = [wrdidx; idx(Swrd.wordOns(idx)>0)];
wrdidx = wrdidx(1:maxnum*2);

xdata = -0.5:0.01:0.5;
% show spectrogram, and envelope and text
for wrd = wrdidx'
    figure;
    trial = Swrd(wrd, :);
    aud = trial.aud{1};
    env = trial.env{1};
    % combine precPhns and onsPhns with a | in the middle
    % text = [join(num2phns(trial.precPhns{1}, timit_details), '') '|' ...
    %     join(num2phns(trial.onsPhns{1}, timit_details), '')];
    % remove all 'cl'
    text = join(trial.phns{1}, '');
    % plot spectrogram on the left and envelope on the right
    subplot(1, 1, 1);
    imagesc(xdata, 1:80, aud);
    colormap(flipud(gray));
    brighten(0.5);
    xline(0, 'LineWidth', 2, 'Color', 'r');
    xlabel('Time (s)');
    xticks([-0.5 0 0.5]);
    yticks([1 80]);
    yticklabels({'0', '8'});
    ylabel('kHz');
    
    set(gca, 'YDir', 'normal', 'FontSize', 13);

    yyaxis right;
    plot(xdata, env, 'Color', 'k', 'LineWidth', 3);
    yticks([]);

    xlim([-0.2 0.2]);
    xticks([-0.2 0 0.2]);
    
    % title(text);
    % subplot(1, 2, 2);
    % plot(xdata, env, 'Color', 'k', 'LineWidth', 2);
    % xline(0, 'LineWidth', 2, 'Color', 'r');
    % set(gca, 'FontSize', 15);
    % xlabel('Time (s)');
    % xlim([-0.2 0.2])
    % xticks([-0.2 0 0.2]);
    % box off;

    sgtitle([trial.precword{1} ' | ' trial.currword{1} ' ' ...
        'AAI=' num2str(trial.ambiguity(1))]);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% B - Neural word boundary decoding split by AAI (only monolinguals)

figure;
% fieldname = 'aud';
% window = 20;
% Dwrd.ambiguity = getAAI(window, Dwrd, fieldname, 20);
% TDwrd.ambiguity = getAAI(window, TDwrd, fieldname, 20);

% red for unfamiliar language, blue for native
colors = [0.8 0.2 0.2;0.2 0.2 0.8];
corpi = {'DIMEX', 'TIMIT'};
Swrds = {Dwrd, TDwrd};

for cr = 1:2
    corpus = corpi{cr};
    Swrd = Swrds{strcmpi(corpus, 'TIMIT')+1};
    time_label = '600ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

    % analysis is sensitive to these binedges!
    binedges = [0 prctile(Swrd.ambiguity, 20) prctile(Swrd.ambiguity, 60) 1];
    
    subplot(1, 2, strcmpi(corpus, 'TIMIT')+1);
    for bin = [1, 3]
        
        diffscore_stats = cell(2, 1);
        for ls = [1, 2]
            ambiguity = discretize(Swrd.ambiguity, binedges);

            scores = decode_details.score{ls};
            yidx = decode_details.yidx{ls};
            y = decode_details.ys{ls};
            diffscore = abs(y-scores);

            % box plot of score differences in this bin
            % one subplot per corpus

            % get the scores for this bin
            ambiguity = ambiguity(yidx);
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

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% B - Neural word boundary decoding split by AAI

figure;
% fieldname = 'aud';
% window = 20;
% Dwrd.ambiguity = getAAI(window, Dwrd, fieldname, 20);
% TDwrd.ambiguity = getAAI(window, TDwrd, fieldname, 20);

% red for unfamiliar language, blue for native
colors = [0.8 0.2 0.2; 0.2 0.2 0.8];
corpi = {'DIMEX', 'TIMIT'};
titles = {'Spanish', 'English'};
Swrds = {Dwrd, TDwrd};

numbins = 5;
figure; 
AUC = nan(2, numbins);
for cr = 1:2
    corpus = corpi{cr};
    Swrd = Swrds{strcmpi(corpus, 'TIMIT')+1};
    time_label = '600ms';
    subj = 'monolingual'; % 'mandarinmono', 'bilingual'
    type = '';
    filename =  [corpus '_word_decode_' subj '_' time_label type '.mat']; 
    load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

    % analysis is sensitive to these binedges!
    % do this based on prctiles
    binedges = prctile(Swrd.ambiguity, 0:100/numbins:100);
    for ls = [1, 2]
        % get all the neural data for this bin
        pcaX = decode_details.pcaX{ls};
        yidx = decode_details.yidx{ls};
        y = decode_details.ys{ls};
        ambig = discretize(Swrd.ambiguity(yidx), binedges);
        Mdl =  decode_details.Mdl{ls};
        for bin = 1:numbins
            % get the pcaX for this bin
            pcaX_bin = pcaX(ambig==bin, :);
            y_bin = y(ambig==bin);

            % run the trained Mdl to get the AUC for the bin
            [~, scores] = predict(Mdl, pcaX_bin);
            [~, ~, ~, AUC(ls, bin)] = perfcurve(y_bin, scores(:, 2), 1);
        end
    end
    subplot(1, 2, strcmpi(corpus, 'TIMIT')+1);
    % plot the AUC difference
    if strcmpi(corpus, 'DIMEX')
        bar(AUC(1, :)-AUC(2, :), 'FaceColor',  [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on;
        % plot the line trend
        x = 1:numbins;
        y = AUC(1, :)-AUC(2, :);
        plot(x, y, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'Marker', 'o'); hold on;
    else % native - foreign
        bar(AUC(2, :)-AUC(1, :), 'FaceColor', [0.8 0.8 0.8], 'EdgeColor', 'none'); hold on;
        % plot the line trend
        x = 1:numbins;
        y = AUC(2, :)-AUC(1, :);
        plot(x, y, 'Color', [0.5 0.5 0.5], 'LineWidth', 2, 'Marker', 'o'); hold on;
    end
    % ylim([-0.1 0.1]);
    yticks(0:0.05:0.5);
    xticks([1 numbins]);
    ylabel('Native - Foreign AUC');
    xticklabels({'low', 'high'});
    xlabel('Acoustic ambiguity (AAI)');
    set(gca, 'FontSize', 15);
    box off;
    title(titles{cr}, 'FontWeight', 'normal');
end


%% C - Comparison of AAI weights

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
    clim([-3 3]);
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
        clim([-3 3]);
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


%% D - Single ERPs for high and low AAI instances

% load in the data for a single subject
% EC100 
SID = 'EC100';
el = 21;
%Dwrd = loadDwrdMulti('dimex',  bef, aft, {SID}, dimex_details); % {'EC260', 'EC266'}

window = 20;
dataf = 100;
fieldname = 'aud';
maxtrials = 40;

%Dwrd.ambiguity = getAAI(window, Dwrd, fieldname, 20);

% find high and low AAI instances
[~, idx] = sort(Dwrd.ambiguity, 'descend');
% only include word onsets
idx = idx(Dwrd.wordOns(idx)>0);
wrdidx = idx(1:maxtrials);

% find 5 highest and 5 lowest AAI that are also word onsets
[~, idx] = sort(Dwrd.ambiguity, 'ascend');
idx = idx(Dwrd.wordOns(idx)>0);
wrdidx = [wrdidx; idx(1:maxtrials)]';

% plot the ERPs for these instances on a single plot
% offset each trace by 0.5
figure;
ctr = 0;
plotctr = 0;
for wrd = wrdidx
    trial = Dwrd(wrd, :);
    ctr = ctr+1;
    if isempty(trial.(SID){1})
        continue;
    end
    cols = puor(10);
    resp = trial.(SID){1}(el, :);
    x = ((1:length(resp))/dataf)-0.5;
    if ctr>maxtrials
        plot(x, resp+plotctr, 'Color', cols(3, :), 'LineWidth', 2); hold on;
    else
        plot(x, resp+plotctr, 'Color', cols(7, :), 'LineWidth', 2); hold on;
    end
    plotctr = plotctr+4;
end
xline(0, 'LineWidth', 2, 'Color', 'r');
xlim([-0.2 0.4]);
xlabel('Time (s)');


%% ----------------------- Supplementary Figures --------------------------
%% 
%% ----------------------------- UNUSED Panels ---------------------------
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


%% ------------- FIX: ACOUSTIC CLASSIFICATION OF WORD BOUNDARY -----------------
%% Logistic on aud/envelope for different times preceding onset

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


%% ------------------ functions -------------------------------------------

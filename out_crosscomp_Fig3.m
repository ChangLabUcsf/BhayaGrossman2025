% <% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;
SIDs = [sSIDs eSIDs bSIDs];
% tps = 50:55;

% selected electrodes
timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

% before and after word time points
bef=50;
aft=50;

% load word structures
if ~exist('Dwrd', 'var')
    load("data/Figure3/Figure3_DIMEXWrd.mat");
end

if ~exist('TDwrd', 'var')
    load("data/Figure3/Figure3_TIMITWrd.mat");
end

load('data/Figure2/Figure2_WordUniVar.mat');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% A/B - Example of segmentation difficulty

% English sentence example
% find English and Spanish sentences that are roughly the same length
timin = [timit_details.sentdet.duration];
dimin = [dimex_details.sentdet.duration];
submat = minus(timin, dimin');
%[~,idxs]=mink(abs(submat(:)), 10); % to find good candidates
idxs = 229738;

% saving the example aud/envelope
corpdetails = {dimex_details, timit_details};
exidx = {[8, 4], [6, 2]};
examples = [struct(); struct()];
for idx = idxs'
    [row,col]=ind2sub(size(submat),idx);
    
    % plot the spectrogram, the word boundaries, and the syllable boundaries
    figure;
    idxs = [row, col];
    for lang = 1:2
        subplot(2, 1, lang);
        details = corpdetails{lang};
    
        % get the sentence
        sent = details.sentdet(idxs(lang));
    
        xdata = (1:size(sent.aud, 2))/sent.dataf - sent.befaft(1)-0.05;
        imagesc(xdata, 1:80, sent.aud);
        brighten(0.5);
        set(gca, 'YDir', 'normal');
        yticks([1 80]);
        yticklabels({'0', '8'});
        ylabel('kHz');
        colormap(flipud(gray));
        hold on;

        yyaxis right;
        plot(xdata, sent.loudness, 'LineWidth', 2, 'Color','k');
        %yticks([]);
    
        wrdonset = find(sent.wordOns>0);
        sylonset = find(sent.syltype>0);
        xline(xdata(wrdonset), 'Color', 'b', 'LineWidth', 2);
        % only plot syllable boundaries if they are different from word boundaries
        xline(xdata(setdiff(sylonset, wrdonset)), 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
        xlim([0 2.5]);
        xticks(0:2);
        xlim([0 2.1])
        set(gca, 'FontSize', 13);
        box off;

        % get the aud and envelope centered at the exidx
        ctr = 1;
        for ex = exidx{lang}
            bef = 0.3 * sent.dataf;
            aft = 0.3 * sent.dataf;
            onstp = sylonset(ex);
            examples(lang).aud{ctr} = sent.aud(:, onstp-bef:onstp+aft);
            examples(lang).env{ctr} = sent.loudness(:, onstp-bef:onstp+aft);
            examples(lang).wordOns{ctr} = ismember(onstp, wrdonset);
            ctr=ctr+1;
        end

        % remove underscores, '.bn' and numbers from the word list
        wordlist = regexprep(strrep(sent.wordList, '_', ' '), {'\d', '\.bn'}, '');
        title(join(wordlist, ' '), 'FontSize', 15, 'FontWeight', 'normal');
    end
    %sgtitle(idx)
end

% plot example syllable and word onsets for the two example sentences
for lang = 1:2
    figure;
    for ex = 1:2
        subplot(1, 2, ex);
        xdata = -0.3:0.01:0.3;
        imagesc(xdata, 1:80, examples(lang).aud{ex});
        brighten(0.5);
        set(gca, 'YDir', 'normal');
        yticks([1 80]);
        yticklabels({'0', '8'});
        ylabel('kHz');
        colormap(flipud(gray));
        hold on;

        yyaxis right;
        plot(xdata, examples(lang).env{ex}, 'LineWidth', 2, 'Color','k');
        yticks([]);
    
        if examples(lang).wordOns{ex}
            xline(0, 'Color', 'b', 'LineWidth', 2);
        else
            xline(0, 'Color', 'b', 'LineWidth', 2, 'LineStyle', '--');
        end
        xlim([-0.15 0.3]);
        xticks(-0.2:0.1:0.2);
        set(gca, 'FontSize', 13);
        box off;
    end
    %sgtitle(idx)
end


%% C - Acoustic word-boundary decoding (takes a minute to run)

timing = 20;
tps = 51-round(timing/2):51+round(timing/2); % around  
disp(['timing: ' num2str(timing)])

nreps = 20;
label = {'Spectrogram'};

% timing = 50;
AUC_all = nan(2, nreps);
acc_all = nan(2, nreps);
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
        1900);
    X_aud{ctr} = cat(3, Swrd.aud{envidx{ctr}});
    y{ctr} = Swrd.wordOns(envidx{ctr})>0;
    ctr=ctr+1;
end

ctr = 1;
AUC_perm = cell(2, 1);
for i = {Dwrd, TDwrd}
    acc = nan(length(timing), nreps); % test accuracy
    weights = cell(length(timing), 1); % aud weights
          
     X = reshape(X_aud{ctr}(:, tps, :), 80*length(tps), []);
     pcaflag = 1;
    [~, ~, AUC_all(ctr, :), Xaud, ~, acc_all(ctr, :), ...
        weights{ctr}, ~, Mdl, comp, mu] = logistic(X', y{ctr}, pcaflag, [], tps, nreps); 

    [AUC_perm{ctr}] = logistic_permute(X', y{ctr}, pcaflag, [], tps, nreps, 25);
    
    avg = {cat(3, Swrd.aud{Swrd.wordOns & ~Swrd.sentOns}); ...
        cat(3, Swrd.aud{~Swrd.wordOns & ~Swrd.sentOns})};
    ctr = ctr + 1;
end
AUC_perm = [AUC_perm{1}; AUC_perm{2}];

% plot the auc values with one boxplot, but scatter spanish and english as filled and empty circles
figure('Position', [100, 300, 550, 300], 'Renderer', 'painters');

% combine AUC and perm boxcharts
boxchart([ones(length(AUC_all(:)), 1); 2*ones(length(AUC_perm(:)), 1)], ...
    [AUC_all(:); AUC_perm(:)], 'BoxFaceColor','none', ...
    'BoxEdgeColor', 'k', 'HandleVisibility', 'off', 'MarkerStyle', 'none'); hold on;

% boxcharts
% boxchart(AUC_all(:), 'BoxFaceColor','none', 'BoxEdgeColor', 'k', ...
%     'HandleVisibility', 'off', 'MarkerStyle', 'none'); hold on;
disp(['Median AUC = ' num2str(median(AUC_all(:)))]);
findobj(gca, 'Tag', 'Box'); hold on;
% scatter only the spanish AUC (jittered x values)
scatter(randn(nreps, 1)*0.05+0.9, AUC_all(1, :), 25, [0.5 0.5 0.5], 'filled');
% scatter only the english AUC (jittered x values)
scatter(randn(nreps, 1)*0.05+1.1, AUC_all(2, :), 25, 'k');

% plot the permuted AUC values
% boxchart(AUC_perm(:), 'BoxFaceColor','none', ...
%     'HandleVisibility', 'off', 'MarkerStyle', 'none'); hold on;

% Formatting
ylim([0.4 1]);
yline(0.5, 'Color', 'k', 'HandleVisibility', 'off');
yline(50, '--k', 'LineWidth', 2, 'HandleVisibility', 'off');
ylabel('AUC');    
legend({'Spanish', 'English'}, 'Location', 'southeast');
box off;
xlim([0.5 2.5]);
set(gca, 'FontSize', 13);
yticks(0.5:0.25:1);
xticks([]);
ylim([0.3 1]);
xticks([1, 2])
xticklabels({'true', 'permuted'})

disp(['Spanish accuracy = ' num2str(mean(acc_all(1, :)))]);
disp(['English accuracy = ' num2str(mean(acc_all(2, :)))]);
disp(['mean accuracy = ' num2str(mean(acc_all(:)))])

[h, p] = ttest2(AUC_all(1, :), AUC_all(2, :));
disp(['p-value: ' num2str(p) ' for ttest2']);

% p test for permutation
pval = sum(AUC_perm(:) > median(AUC_all(:))) / length(AUC_perm(:));
disp(pval);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% plotWORDERP: D/E - Native Brain and Example Word-Syllable ERP

% find top 3 weighted electrodes and look at word erps
f = figure; 
pthresh = 0.01;

% Spanish
% SIDs = {'EC260', 'EC260', 'EC260', 'EC260', 'EC260', 'EC260', ...
%     'EC260', 'EC260'}; % 'EC266', 'EC266'
% els = [ 221, 192, 205, 206, 207, 208, 217, 214]; %236, 205 for EC260, 178, 175
% Ec183
% SIDs = {'EC100'};% 'EC100', 'EC100', 'EC100','EC100'};
% els = [21];% 21,135, 22, 69, 70, 71, 118, 150]; %236, 205 for EC260, 178, 175

% English
% SIDs = {'EC100'};% 'EC100', 'EC100', 'EC100','EC100'};
% els = 132;%, 56, 135, 151, 22, 70, 71, 150]; %236, 205 for EC260, 178, 175

SIDs = {'EC183'};% 
els = 71;%

% els = [185];%
% SIDs = repmat({'EC163'}, length(els), 1);
% 
% els = [42, 57, 58, 72:74, 88];%
% SIDs = repmat({'EC172'}, length(els), 1);

% els = [132]; % 132, 135, 151, 22, 70, 71, 150, 236, 205
% SIDs = repmat({'EC100'}, length(els), 1);

% els = [71]; % , 56, 135, 151, 22, 70, 71, 150, 236, 205
% SIDs = {'EC183'};

% EC186, 207

% Bilingual
% EC163 - el55
% SIDs = {'EC260'};% 'EC100', 'EC100', 'EC100','EC100'};
% els = 221;%,135, 22, 70, 71, 150]; %236, 205 for EC260, 178, 175
% 
% % in figure
% % EC100: 21, 71 EC183: 135, 71
% EC186, 207
% 
% els = [55]; %236, 205 for EC260, 178, 175
% SIDs = repmat({'EC163'}, length(els), 1);

% Dwrd = loadDwrdMulti('dimex',  bef, aft, SIDs, dimex_details);
% TDwrd = loadDwrdMulti('timit',  bef, aft, SIDs, timit_details);

plotSingleTrial = 0;
numel = length(els);
Swrds = {Dwrd, TDwrd};
tps = 30:70;
x = -0.2:0.01:0.2;
for s = [2, 1]
    Swrd = Swrds{s};

    for ctr = 1:length(els)
        SID = SIDs{ctr};    
    
        dummy = struct();
        nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
            %| (Swrd.syll<2 & ~isnan(Swrd.syll));
        dummy.(SID).resp = cat(3, Swrd.(SID){~nanidx});

        dummy.(SID).resp  = dummy.(SID).resp(:, tps, :);
        dummy.wordOns = Swrd.wordOns (~nanidx);
        dummy.syllOns = ones(sum(~nanidx), 1);

        ls = wordsurp_encoding.ls(find(strcmp(wordsurp_encoding.SID, SIDs{1}), 1));
        if s == ls   
            cols = [0 0 0; 0.1, 0.1 0.9]; 
        else
            cols = [0 0 0; 0.9, 0.1 0.1]; 
        end

        subplot(2, numel, ctr+(s-1)*numel)
        rmpath(genpath('shadederror'));
        plotWordErp(dummy, SID, els(ctr), ...
            [], f, cols, 1, 0.2, 1); hold on;
        ylabel('HFA (z)');
        set(gca, 'FontSize', 13);

        nanidx = squeeze(sum(isnan(dummy.(SID).resp(els(ctr), :, :))));
        [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
            dummy.(SID).resp(els(ctr), :, ~nanidx), dummy.wordOns(~nanidx)+1, [1, 2]);
        corrected_pval = pthresh / size(dummy.(SID).resp, 2);
        fthresh = finv(1-corrected_pval, df1, df2);  
    
        scatter(x(fvals>fthresh), 0.1*ones(1, sum(fvals>fthresh)), 45, ...
            fvals(fvals>fthresh), 'filled', 'HandleVisibility', 'off');
        cm = colormap("gray");
        colormap(flipud(cm(1:200, :)))
        % ylim([0 1.3]);
        ylim([-0.2 0.5]);
        xlim([-0.2, 0.2]);
        h=xline(0);
        h.Color = 'k';
        legend('off');

        if plotSingleTrial
            numtrials = 100;
            
            wordOnsResp = squeeze(dummy.(SID).resp(els(ctr), :, logical(dummy.wordOns)));
            syllOnsResp = squeeze(dummy.(SID).resp(els(ctr), :, ~logical(dummy.wordOns)));

            % remove all nans
            wordOnsResp = wordOnsResp(:, ~any(isnan(wordOnsResp), 1));
            syllOnsResp = syllOnsResp(:, ~any(isnan(syllOnsResp), 1));

            % find syllable trials most similar to the average
            [~, idx] = sort(arrayfun(@(x) corr(mean(syllOnsResp, 2, 'omitnan'), ...
                syllOnsResp(:, x), 'Type', 'Spearman'), 1:size(syllOnsResp, 2)));

            figure;
            subplot(1, 2, 1); 
            xdata = -0.5:0.01:0.5;
            imagesc(xdata, 1:numtrials, syllOnsResp(:,idx(1:numtrials))'); hold on;
            clim([-5 5]);
            yticks([1 100]);
            ylabel('trials');
            yyaxis right; 
            plot(xdata, mean(syllOnsResp, 2, 'omitnan'), 'Color', 'k', 'LineWidth', 2.5);
            title('Syllable');
            xline(0, 'Color', 'k', 'LineWidth', 2.5);
            xlim([-0.2, 0.4]);
            ylim([0.2 1]);
            yticks([0 0.5 1]);
            colormap(flipud(prgn));
            set(gca, 'FontSize', 13);
            
            % word trials
            [~, idx] = sort(arrayfun(@(x) corr(mean(wordOnsResp, 2, 'omitnan'), ...
                wordOnsResp(:, x), 'Type', 'Spearman'), 1:size(wordOnsResp, 2)));

            subplot(1, 2, 2);
            imagesc(xdata, 1:numtrials, wordOnsResp(:,idx(1:numtrials))'); hold on;
            yticks([1 100]);
            ylabel('trials');
            clim([-5 5]);
            % make imagesc lighter
            yyaxis right; 
            plot(xdata, mean(wordOnsResp, 2, 'omitnan'), 'Color', cols(2, :), 'LineWidth', 2.5);
            xline(0, 'Color', 'k', 'LineWidth', 2.5);
            xlim([-0.2, 0.4]);
            ylim([0.2 1]);
            yticks([0 0.5 1]);
            title('Word');
            set(gca, 'FontSize', 13);
        end
    end
    clear dummy
end

% initialize design electrode structure
fieldnames = {'Spanish', 'English'};
fields = {'sp_uv_all', 'eng_uv_all', '', 'sp_uv_all'}; 
% uv feature order
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant', 'surp', 'word'

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
        binedges = [-1 0.0005:0.005:0.01 0.015:0.015:0.1];
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
        end

        desel.(SIDs{1}).selid = els;
        ls = wordsurp_encoding.ls(find(strcmp(wordsurp_encoding.SID, SIDs{1}), 1));
        if lang == ls   
            cls = flipud(blues(8));
        else
            cls = flipud(reds(8));
        end
        desel.cols = [0 0 0;cls(3:end, :)];

        % just for native brain and coverage
        % desel.(SIDs{1}).elid=[];
        % desel.(SIDs{1}).condition=[];

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

        alpha 0.8;
        % add a pie
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
    betaInfo* *encoding* *wrd*;

%% EC186: F - Single electrode quantity of Word-Syllable difference (takes a min to run)

% need to preload all participants into Dwrd and TDwrd
% TDwrd = loadDwrdMulti('timit', bef, aft, [eSIDs sSIDs], timit_details);
% Dwrd = loadDwrdMulti('dimex',  bef, aft, [eSIDs sSIDs], dimex_details);

% load in all speech responsive electrodes into a table
electbl = loadSpeechResponsive(SIDs);
% subset electbl to only sSIDs and eSIDs
electbl = electbl(ismember(electbl.SID, [eSIDs, sSIDs]), :);

% for each electrode, save the fvals and the longest contiguous significant window
for swrd = {Dwrd, TDwrd}

    Swrd = swrd{1};
    % figure out which corpus
    if startsWith(Swrd.sentid{1}, 's')
        corpus = 'dimex';
    else
        corpus = 'timit';
    end

    % define all fields to add
    fvalfield = strcat(corpus, '_fvals');
    fthreshfield = strcat(corpus, '_fthresh');
    longestsigfield = strcat(corpus, '_contigsig');

    % for all unique SIDs
    for s = [eSIDs, sSIDs]
        SID = s{1};
        sidx = strcmp(electbl.SID, SID);

        dummy = struct();
        nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
            %| (Swrd.syll<2 & ~isnan(Swrd.syll));
        dummy.(SID).resp = cat(3, Swrd.(SID){~nanidx});
        
        dummy.wordOns = Swrd.wordOns (~nanidx);
        dummy.syllOns = ones(sum(~nanidx), 1);

        % find electrodes in electbl for this sid
        els = electbl.el(sidx);

        [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
            dummy.(SID).resp(els, :, :), dummy.wordOns+1, [1, 2]);
        corrected_pval = 0.01 / size(dummy.(SID).resp, 2);
        fthresh = finv(1-corrected_pval, df1, df2); 
        
        electbl.(fvalfield)(sidx) = mat2cell(fvals, ones(length(els), 1), ...
            size(fvals, 2));
        electbl.(fthreshfield)(sidx) = mat2cell(fvals>fthresh, ones(length(els), 1), ...
            size(fvals, 2));

        % find the longest contiguous window
        tmp = arrayfun(@(x) max(diff(find(fvals(x, :)<=fthresh))), ...
            1:length(els), 'UniformOutput', false);
        % collapse into an array, if a cell is empty replace with 0
        sig = zeros(length(els), 1);
        for i = 1:length(els)
            if ~isempty(tmp{i})
                sig(i) = tmp{i};
            end
        end
        electbl.(longestsigfield)(sidx) = sig;
    end
end
% get the native and foreign contiguous significant difference
electbl.natcontig = [electbl.dimex_contigsig(ismember(electbl.SID, sSIDs)); ...
    electbl.timit_contigsig(ismember(electbl.SID, eSIDs))];
electbl.forcontig = [electbl.timit_contigsig(ismember(electbl.SID, sSIDs)); ...
    electbl.dimex_contigsig(ismember(electbl.SID, eSIDs))];

N = histcounts2(electbl.natcontig, electbl.forcontig, 0:2:40, 0:2:40);
figure;
imagesc(N);
set(gca, 'YDir', 'normal');
clim([0 30]);
colormap([1 1 1; flipud(inferno(30))]);

% remove all electrodes with no significant windows
electbl = electbl(electbl.dimex_contigsig>5 | electbl.timit_contigsig>5, :);

figure('Position', [100, 300, 550, 300], 'renderer', 'painters');
% make boxchart for each condition
idx = electbl.natcontig>5 & electbl.forcontig>5;
boxchart(ones(size(electbl.natcontig(idx), 1), 1), [electbl.natcontig(idx)], ...
                    'BoxFaceColor', [0.5 0.5 0.9], 'JitterOutliers', 'on', ...
                    'MarkerStyle', '.', 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.5); hold on;
boxchart(ones(size(electbl.forcontig(idx), 1), 1)*2, [electbl.forcontig(idx)], ...
                    'BoxFaceColor', [0.9 0.5 0.5], 'JitterOutliers', 'on', ...
                    'MarkerStyle', '.', 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.5); hold on;
% formatting
ylim([0 40]);
xlim([0.5 2.5]);
yticks(0:20:40);
yticklabels({'0', '0.2', '0.4'});
ylabel('Contiguous signifcant difference (s)');
xticks([1 2]);
xticklabels({'Native', 'Foreign'});
set(gca, 'FontSize', 13);
box off;

% Statistical testing with linear mixed effect model
tbl=table();
% first stack the SID, el, and natcontig, then stack the SID, el, and forcontig
tbl.SID = [electbl.SID; electbl.SID];
tbl.el = [electbl.el; electbl.el];
tbl.ls = [electbl.ls; electbl.ls];
tbl.contig = [electbl.natcontig; electbl.forcontig];
tbl.native = [ones(height(electbl), 1); zeros(height(electbl), 1)];
lme = fitlme(tbl, 'contig ~ native + ls + (1|hemi) + (1|SID) + (1|el:SID)');
disp(lme);
text(1.5, 35, getSigStr(lme.Coefficients.pValue(3), 1), 'FontSize', 20);

% show all the data points                    
figure;
scatter(ones(size(electbl.natcontig, 1), 1), electbl.natcontig, 25, 'filled', 'b'); hold on;
scatter(2*ones(size(electbl.forcontig, 1), 1), electbl.forcontig, 25, 'filled', 'r'); hold on;
for i = 1:size(electbl, 1)
    plot([1 2], [electbl.natcontig(i), electbl.forcontig(i)], 'Color', [0.5 0.5 0.5]);
end
plot([1 2], [median(electbl.natcontig), median(electbl.forcontig)], 'Color', 'k', 'LineWidth', 2);

% show a pie chart with the number of electrodes that show signifcant differences for only native, only foreign, and both
figure;
p = pie([sum(electbl.natcontig>5 & electbl.forcontig<=5), ...
    sum(electbl.natcontig<=5 & electbl.forcontig>5), ...
    sum(electbl.natcontig>5 & electbl.forcontig>5)], [1 1 1]);
p(1).FaceColor = [0.3 0.5 0.9];
p(1).EdgeColor = 'none';
p(3).FaceColor = [0.9 0.5 0.5];
p(3).EdgeColor = 'none';
p(5).FaceColor = [0.7 0.6 0.9];
p(5).EdgeColor = 'none';
p(2).FontSize = 20;
p(4).FontSize = 20;
p(6).FontSize = 20;
legend({'Native', 'Foreign', 'Both'}, 'Location', 'southeast');

%% ----------------------- Supplementary Figures --------------------------
%% S9A - Neural word boundary (both Spanish and English)

figure('Renderer', 'painters');
Swrds = {Dwrd, TDwrd};
titles = {'Spanish speech', 'English speech'};
decode_type = 'word'; % 'word' 
cols = {'r', 'b'};
ctr=1;
for s = Swrds
    Swrd = s{1};
    time_label = '600ms';
    % time_label = '10ms-sliding';
    lss = [1, 2];
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
        load([datapath 'Figure3/decode/' filename], 'decode_details');
    
        styles = {':', '-', '', '-'};
        switch time_label
            case '10ms-sliding'
                accls(ls, :, :) = squeeze(decode_details.AUC(ls, :, :));
            otherwise
                accls(ls, :) = squeeze(decode_details.AUC(ls, :, :));
        end 
    end
    
    if ctr == 1 % dimex
        position = [1, 2];
    else
        position = [2, 1];
    end

    for ls = lss
        h=boxchart(ones(size(accls, 2), 1)*position(ls), ...
            accls(ls, :), 'BoxEdgeColor', 'k', 'BoxFaceColor', 'w'); 
        h.JitterOutliers = 'off';
        h.MarkerStyle = '.';
        h.MarkerColor = colors(ls, :);
        hold on; % 'BoxStyle', 'filled'
        rng(1);
        scatter(randn(size(accls, 2), 1)*0.1-0.05+position(ls), accls(ls, :), 25, 'filled', ...
            cols{(ls==ctr)+1})
    end  

    xticks([1, 2]);
    xticklabels({'Native', 'Foreign'});
    xlim([0.5 2.5]);
    legend('off');
    
    combo = [1, 2;];
    pos = nan(size(combo, 1), size(combo, 2));
    for r = 1:size(combo, 2)
        for c = 1:size(combo, 1)
            pos(c, r) = position(combo(c, r));
        end
    end
     
    [~, p] = ttest2(accls(pos(1), :), accls(pos(2), :), "Tail", "right");
    if ~isempty(getSigStr(p, 2))
        text(mean(pos)-0.25, 0.8, getSigStr(p, 2), 'FontSize', 15);
    end
    yticks();
    yline(0.5, 'Color', 'k');
   
    ylim([0.45 0.9]);
    yticks(0.5:0.2:0.9);
    ylabel('AUC');
    set(gca, 'FontSize', 15);
    title(titles{ctr}, 'Fontweight', 'normal')
    box off;
    ctr = ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% S9B - Neural word boundary (temporally resolved)
% run from inside crosslang folder
titles = {'Spanish speech', 'English speech'};
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
plotBil = 0;
for c = {'DIMEX', 'TIMIT'}
    
    corpus = c{1};
    disp(['-----------' corpus '--------------']);
    filename = [corpus '_word_decode_monolingual_10ms-sliding' subtype '.mat']; 
    load([datapath 'Figure3/decode/' filename], 'decode_details');
    x = (decode_details.startp-50)./100;

    subplot(2, 1, ctr)
    colors = [1 0 0; 0 0 1];
    data = nan(2, length(x), size(decode_details.AUC, 3));
    for ls = 1:2
        % data(ls, :, :) = smoothdata(squeeze(decode_details.AUC(ls, :, :)), ...
        %     'SmoothingFactor', 0.1);
        data(ls, :, :) = squeeze(decode_details.AUC(ls, :, :));
        y = squeeze(mean(data(ls, :, :), 3));
        err = squeeze(nansem(data(ls, :, :), 3));
        shadedErrorBar(x, y, err, ...
            'lineprops', {'color', colors((ls==ctr)+1, :), 'linewidth', 2, ...
            'linestyle', linestyle}); hold on;

        elecs = decode_details.elecs{ls};
        disp(['ls: ' num2str(ls) ', elecs:' num2str(length(elecs)) ...
            ', num subj: ' num2str(length(unique(decode_details.encoding.SID(elecs))))]);

        % print max peak +- values
        maxtp = argmax(mean(squeeze(data(ls, :, :)), 2));
        %disp(datastats(data(ls, maxtp, :)))
    end

    X = reshape(cat(2, squeeze(data(1, :, :)), squeeze(data(2, :, :))), 1, ...
        length(x), size(data, 3)*2);
    y = [ones(1, size(data, 3)) 2*ones(1, size(data, 3))];
    [fvals, ~, ~, df1, df2] = Fstat_TIMIT(X, y, [1, 2]);

    % Bonferroni corrected
    fthresh = finv(1-(0.01/size(data, 2)), df1, df2);  

    scatter(x(fvals>fthresh), 0.52*ones(1, sum(fvals>fthresh)), 45, ...
        int32(fvals(fvals>fthresh))', 'filled', 'HandleVisibility', 'off');
    cm = colormap("gray");
    colormap(flipud(cm(1:200, :)))

    ylim([0.45 0.8]);
    xlim([-0.2 0.4])
    yticks((0.5:0.2:0.7));
    set(gca, 'FontSize', 15);
    
    h=title(titles{ctr}, 'Fontweight', 'normal');
    xline(0, '-k', 'linewidth', 1, 'handlevisibility', 'off');
    yline(0.5, '-k', 'linewidth', 1, 'handlevisibility', 'off');
    xlabel('Time (s)')
    ylabel('AUC');
    ctr=ctr+1;
end
legend({'foreign', 'native'}, 'location', 'northwest');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% S9C - Neural word boundary (single subjects across participants)

timelabel = '600ms';
subj = 'monolingual';

tmp = load([datapath ...
        'Figure3/decode/DIMEX_word_decode_' subj '_' timelabel '_bysubj.mat'], ...
        'decode_details');
dimex=tmp.decode_details.tbl;

tmp = load([datapath ...
        'Figure3/decode/DIMEX_word_decode_' subj '_' timelabel '.mat'], ...
        'decode_details');
dimex_grp = tmp.decode_details;

tmp = load([datapath ...
        'Figure3/decode/TIMIT_word_decode_' subj '_' timelabel '_bysubj.mat'], ...
        'decode_details');
timit=tmp.decode_details.tbl;

tmp = load([datapath ...
        'Figure3/decode/TIMIT_word_decode_' subj '_' timelabel '.mat'], ...
        'decode_details');
timit_grp=tmp.decode_details;

joined = join(timit, dimex, 'Keys', 'SID');

% calculate median AUC for timit_auc and dimex_auc, calculate numel elecs
for i = 1:height(joined)
    joined.numel_timit(i) = length(joined.elecs_timit{i});
    joined.numel_dimex(i) = length(joined.elecs_dimex{i});
    joined.median_auc_timit(i) = median(joined.auc_timit{i});
    joined.median_auc_dimex(i) = median(joined.auc_dimex{i});
end

% plot the average AUC for each subject in the two corpora (Spanish and English)
langs = {'dimex', 'timit'};
titles = {'Spanish speakers', 'English speakers'};
grp_decode = {dimex_grp, timit_grp};
cols = [0 0 1; 1 0 0];

figure; 
% for each language group
for lss = 1:2
    subplot(2, 1, lss);
    % subject indices
    sidx = find(joined.ls_timit==lss);
    auc = nan(2, length(sidx)); % 0 is native, 1 is foreign
    sz = nan(2, length(sidx));

    if lss == 1
        % native is dimex
        auc(1, :) = joined.median_auc_dimex(sidx);
        auc(2, :) = joined.median_auc_timit(sidx);
        sz(1, :) = joined.numel_dimex(sidx)*1.5;
        sz(2, :) = joined.numel_timit(sidx)*1.5;
    else
        % native is timit
        auc(2, :) = joined.median_auc_dimex(sidx);
        auc(1, :) = joined.median_auc_timit(sidx);
        sz(2, :) = joined.numel_dimex(sidx)*1.5;
        sz(1, :) = joined.numel_timit(sidx)*1.5;
    end

    % scatter plot of AUC colored by condition, with lines connecting native and foreign for the same subject
    for i = 1:length(sidx)

        % jitter the x values
        x = [1 2] + randn(1)*0.1 - 0.05;
        scatter(x, auc(:, i), sz(:, i), cols, 'filled', ...
            'MarkerFaceAlpha', 0.6, 'LineWidth', 2);
        hold on;
        plot(x, auc(:, i), 'Color', [0.5 0.5 0.5], 'LineWidth', 0.5);
    end

    % make a box plot ontop of scatter
    boxchart([ones(1, length(sidx)), 2*ones(1, length(sidx))], ...
        [auc(1, :), auc(2, :)], 'BoxFaceColor', 'w', 'BoxEdgeColor', 'k', ...
        'JitterOutliers', 'off', 'MarkerStyle', 'none', 'MarkerColor', 'k'); hold on;
    ylim([0.45 0.75]);
    xlim([0.5 2.5]);
    yline(0.5, 'k');
    yticks([0.5 0.7]);
    ylabel('AUC')
    xticks([1, 2]);
    xticklabels({'Native', 'Foreign',});
    set(gca, 'FontSize', 13);
    title(titles{lss}, 'Fontweight', 'normal');

    % two sample ttest
    [p, ~] = signrank(auc(1, :), auc(2, :), 'Tail', 'right');
    disp(['one-tailed paired sign rank ' num2str(p) ', observations = ' num2str(size(auc, 2))]);
    tbl = table();
    tbl.auc = [auc(1, :), auc(2, :)]';
    tbl.native = [ones(1, length(sidx)), 0*ones(1, length(sidx))]';
    tbl.numel = [sz(1, :)/1.5, sz(2, :)/1.5]';
    lme2 = fitlme(tbl, 'auc ~ native + (1|numel)');
    %disp(lme2)
    disp(['From linear mixed effects: ' num2str(lme2.Coefficients.pValue(2), 2)])
    text(1.5, 0.65, getSigStr(p, 2), 'FontSize', 13);
end

% add median and range of trials and electrodes across conditions
med_dimex = median(joined.numel_dimex);
range_dimex = [min(joined.numel_dimex), max(joined.numel_dimex)];
med_timit = median(joined.numel_timit);
range_timit = [min(joined.numel_timit), max(joined.numel_timit)];
disp(['Median number of electrodes for DIMEX: ' num2str(med_dimex) ...
    ', range: ' num2str(range_dimex)]);
disp(['Median number of electrodes for TIMIT: ' num2str(med_timit) ...
    ', range: ' num2str(range_timit)]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% LOGISTIC: S10A - Acoustic cues to word boundary (boundary decoding)
% timing = [5, 10, 20, 50]; 
timing = 20;
nreps = 20;
mintrials = 1905;
fields = {'env', 'formant', 'consphnfeat'};
pcaflags = [0 0 0 1];
Swrds = {Dwrd, TDwrd};
corpus_details = {dimex_details, timit_details};

AUC_all = cell(1, 2);
for lang = 1:2
    Swrd = Swrds{lang};

    % randomly sample 2000 instances to use
    rng(1);
    trlidx = randsample(find(cellfun(@(x) ~isempty(x), ...
        [Swrd.env])), mintrials);
    
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
        boxplot(auc(ls, :), 'Position', i-0.35+0.25*ls, 'Color', ...
            cols{ls}(length(cols{ls})-2, :), 'BoxStyle', ...
            'filled', 'Symbol', ''); hold on;
    end
    [p, ~] = signrank(auc(1, :), auc(2, :));
    text(i, 0.8, getSigStr(p, 1), 'FontSize', 20, 'FontWeight', 'bold');    
end

%formatting
ylim([0.40 0.9]);
yline(0.5, 'Color', 'k');
yticks(0.4:0.2:0.8);
xlabel('Feature');
ylabel('AUC');    
xticks(1:length(fields));
xticklabels(fields)
set(gca, 'FontSize', 13);
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* X_neural;

%% S10B - Acoustic cues to word boundary (cross training)
timing = 20;
mintrials = 1905;
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
        mintrials);
    X_aud{ctr} = cat(3, Swrd.aud{envidx{ctr}});
    y{ctr} = Swrd.wordOns(envidx{ctr})>0;
    ctr=ctr+1;
end

ctr = 1;
for i = {Dwrd, TDwrd}
    AUC = nan(2, nreps);
    acc = nan(2, nreps); % test accuracy
    weights = cell(length(timing), 1); % aud weights
          
    % 1 is aud
     X = reshape(X_aud{ctr}(:, tps, :), 80*length(tps), []);
    [fp_aud, tp_aud, AUC(1, :), Xaud, ~, acc(1, :), ...
        weights{ctr}, ~, Mdl, comp, mu] = logistic(X', y{ctr}, 1, [], tps, nreps); 

    % Generate the cross language predictions
    crossctr = mod(ctr, 2)+1;
    X_cross = reshape(X_aud{crossctr}(:, tps, :), 80*length(tps), []);
    [AUC(2, :), pcaX] = crosstestlogistic(X_cross', y{crossctr}', ...
        Mdl, comp, mu, nreps);
                         
    AUC_all{ctr}=AUC;
    ctr = ctr + 1;
    
    avg = {cat(3, Swrd.aud{Swrd.wordOns & ~Swrd.sentOns}); ...
        cat(3, Swrd.aud{~Swrd.wordOns & ~Swrd.sentOns})};
end
tmp = brewermap(7, 'YlGnBu');
tmp2 = brewermap(7, 'YlOrRd');

% plot the auc values
cols = [{tmp(2:end, :)}, {tmp2(2:end, :)}];
figure; 
for ls = 1:2
    auc = nan(2, 20);
    for feat = [1, 2] 
        auc(feat, :) = squeeze(AUC_all{ls}(feat, :));
        boxplot(auc(feat, :), 'Position', ls-0.35+0.25*feat, 'Color', ...
            cols{ls}(length(cols{ls})-feat*2, :), 'BoxStyle', ...
            'filled', 'Symbol', ''); hold on;
        % feat-0.35+0.25*ls
    end
    [p, ~] = signrank(auc(1, :), auc(2, :), 'Tail', 'right');
    %[h, p] = ttest2(auc(1, :), auc(2, :), 'Tail', 'right');
    text(ls, 0.83, getSigStr(p, 1), 'FontSize', 20, 'FontWeight', 'bold');    
end

% Formatting
difflang = signrank(AUC_all{1}(1, :), AUC_all{2}(1, :));
text(1.5, 0.8, getSigStr(difflang, 1), 'FontSize', 18);
disp(difflang)
yline(0.5, 'Color', 'k', 'HandleVisibility', 'off');
xlabel('Feature');
ylabel('AUC');    
set(gca, 'FontSize', 13, 'Xtick', 1:2, 'Xticklabel', ...
    {'Same-Language', 'Cross-Language' });
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd* X_neural;
%% Set up

% Ilina Bhaya-Grossman
% 01.08.2022
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('util'))
zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;
addpath(genpath(datapath))

bef=20;
aft=50;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs, ~] = getSIDinfo();
SIDs = [sSIDs, eSIDs, {'HS8', 'HS9', 'HS10'}]; % , {'HS11', 'HS9', 'HS10'}
SIDs =[bSIDs];

% Dvow = loadDD('dimex', bef, aft, {}, datapath);
% TDvow = loadDD('timit', bef, aft, {}, datapath);

% load in all beta model versions

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
tps = 50:55;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;
%% -------------------- Electrode selection details -------------------- %%

timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

proportion = nan(length(SIDs), 1);
allelec =  nan(length(SIDs), 1);
interelec =  nan(length(SIDs), 1);
for s = 1:length(SIDs)
    SID = SIDs{s};
    if isfield(timit_elecs.allidx, SID) && isfield(dimex_elecs.allidx, SID)
        allelec(s) = length(union(timit_elecs.allidx.(SID), ...
            dimex_elecs.allidx.(SID)));
        interelec(s) = length(intersect(timit_elecs.allidx.(SID), ...
            dimex_elecs.allidx.(SID)));
        proportion(s) = (interelec(s)/allelec(s))*100;

        disp(SID)
        disp([num2str(interelec(s)) ' ' num2str(allelec(s)) ', ' ...
            num2str(proportion(s)) '%']);
    end
end
disp(median(proportion, 'omitnan'));
disp((sum(interelec, 'omitnan')/sum(allelec, 'omitnan'))*100)
disp(sum(allelec, 'omitnan'))

%% -------------------- AVERAGE SENTENCE RESPONSE -------------------------
%% -------------------- Loading sentence (and stitched) responses  ----- %%

[sent_encoding] = loadSentenceResponse(SIDs, timit_details, dimex_details, datapath);

% load elecs
dimex_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
timit_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');

% Set color scheme
binedges = -0.20:0.02:0.20;
colors = flipud(brewermap(length(binedges)-1, 'Spectral'));
cols = [colors(3, :); colors(end-2, :); colors(round(size(colors, 1)/2), :)];

maxresp = nan(2, height(sent_encoding));
maxtp = nan(2, height(sent_encoding));
type = nan(2, height(sent_encoding));
wndsz = 5; 
sliding = 65:wndsz*2:180;

swind = 150;
ewind = 100;
tempresp = nan(2, length(sliding), height(sent_encoding));
for el = 1:height(sent_encoding)

    % Aggregate all sentences, find mean 
    % (TIMIT sentences are 50ms later -- FIX!!!)
    en_resp = sent_encoding.en_sent_resp{el};
    sp_resp = sent_encoding.sp_sent_resp{el};

    % Make sure baseline period at zero
    mintp = min(size(en_resp, 2), size(sp_resp, 2)); 

    % For each sentence stitch the first 150ms + last 100ms 
    [en_resp, ~] = stitchedResp(en_resp, swind, ewind);
    [sp_resp, ~] = stitchedResp(sp_resp, swind, ewind);

    % clip to be the same length
    factor = 0.1;
    if all(isnan(en_resp)), en_mresp = nan(1, swind+ewind); 
    else
        en_mresp = [smoothdata(mean(en_resp(1, 1:swind, :), 3, 'omitnan'),...
            'gaussian', 'SmoothingFactor', factor), ...
            smoothdata(mean(en_resp(1, swind+1:end, :), 3, 'omitnan'),...
            'gaussian','SmoothingFactor', factor)];
    end
    sent_encoding.en_mresp(el) = {en_mresp};
    
    if all(isnan(sp_resp)), sp_mresp = nan(1, swind+ewind); 
    else
        sp_mresp = [smoothdata(mean(sp_resp(1, 1:swind, :), 3, 'omitnan'),...
            'gaussian', 'SmoothingFactor', factor), ...
            smoothdata(mean(sp_resp(1, swind+1:end, :), 3, 'omitnan'),...
            'gaussian', 'SmoothingFactor', factor)];
    end

    % fvals between languages
    if ~isempty(en_resp) && ~isempty(sp_resp)
        [fvals, betweenVar, withinVar, df1, df2] = ...
            Fstat_TIMIT(cat(3, en_resp, sp_resp), [ones(1, size(en_resp, 3)) ...
            2*ones(1, size(sp_resp, 3))], [1, 2]);
        fthresh = finv(1-0.001, df1, df2);
        sent_encoding.fvals(el) = {fvals};
        sent_encoding.fthresh(el) = fthresh;   
    end

    % align responses so they are as close as possible
    [~, sp_mresp] = procrustes(en_mresp',sp_mresp', "scaling", false);
    sp_mresp = sp_mresp';
    sent_encoding.sp_mresp(el) = {sp_mresp};

    % check if there is any all zero responses
    if ~(all(sp_mresp==0) || all(en_mresp==0))           

        % do sliding window over entire response
        for w = 1:length(sliding)
            wind = max(sliding(w)-wndsz , 1):min(sliding(w)+wndsz , mintp);
            tempresp(:, w, el) = mean([en_mresp(wind); sp_mresp(wind)], 2);
        end

        % find maximum window for English and Spanish response
        prom = 0.5;
        if all(isnan(en_mresp)), maxtp(1, el) = nan;
        else
            [~, loc] = findpeaks(en_mresp, 'MinPeakProminence',prom, ...
                'NPeaks',1, 'SortStr','descend');
            if ~isempty(loc)
                maxtp(1, el) = loc;
            end
%             [~, maxtp(1, el)] = max(en_mresp);           
        end
        if all(isnan(sp_mresp)), maxtp(2, el) = nan;
        else
            [~, loc] = findpeaks(sp_mresp, 'MinPeakProminence',prom, ...
                'NPeaks',1, 'SortStr','descend');
            if ~isempty(loc)
                maxtp(2, el) = loc;
            end
%             [~, maxtp(2, el)] = max(sp_mresp);
        end

        wind = nan(2, wndsz*2+1);
        if ~any(isnan(maxtp(:, el)))  
            for l = 1:2
                tmp = max(maxtp(l, el)-wndsz , 1):min(maxtp(l, el)+wndsz , swind+ewind);
                wind(l, 1:length(tmp)) = tmp;
            end
            % Note: first row is english, second row is spanish
            maxresp(:, el) = [mean(en_mresp(removeNan(wind(1, :))), 2, 'omitnan') ...
                mean(sp_mresp(removeNan(wind(2, :))), 2, 'omitnan')];
        end
            
        % find larger, and look at the ratio between the smaller response and
        % the larger response
        SID = sent_encoding.SID{el};
        if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            type(:, el) = [ismember(sent_encoding.el(el), dimex_elecs.allidx.(SID)) ...
                    ismember(sent_encoding.el(el), timit_elecs.allidx.(SID))];
            clear SID
        end
    end

    debug = 1;
    if debug
        % EC235, 249; EC161, 72; EC105, 251
        if ismember(el, 10:40:700) %94, 97, 100, 298, 316, 165 406 % 863 482
            % 967 1169 346 1049 1253 1301, 1280, 254, 785, 1405, 233, 923, 1403
            figure('Renderer','painters');
            bef = 0.5;                 
            
            xlim([-0.5 2]);
            ylim([-0.5 4.5])

            gap = 20;
            for j = 1:2
                if wind(j, 1)>swind, wind(j, :)=wind(j, :)+gap; end
            end

            sentenceERP(sp_resp, cols(1, :), bef, [swind, ewind], gap); hold on; %[0.19 0.53 0.74]
            if ~all(isnan(wind(2, :)))
                highlightERPWindow(wind(2, :)./100-bef, cols(1, :));
            end

            sentenceERP(en_resp, cols(2, :), bef, [swind, ewind], gap); %[0.83 0.24 0.30]
            if ~all(isnan(wind(1, :)))
                highlightERPWindow(wind(1, :)./100-bef, cols(2, :));
            end
            % average difference between onsets is 50ms

            title(['Language exp: ' num2str(sent_encoding.ls(el)) ...
                ', el: '  num2str(el) ' SID: ' sent_encoding.SID{el}]);     
        end
    end 
end

sent_encoding.maxresp = maxresp';
sent_encoding.maxtp = maxtp';
sent_encoding.type = type';
sent_encoding.tempresp = permute(tempresp, [3, 1, 2]);
sent_encoding.swind = repmat(swind, height(sent_encoding), 1);
sent_encoding.ewind = repmat(ewind, height(sent_encoding), 1);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Selectivity index

% native - nonnative / native + nonnative

% use area under the cure to calculate selectivity index per electrode
for el = 1:height(sent_encoding)
    ls = sent_encoding.ls(el);
    if ls == 1
        % Spanish
        unf_resp = sent_encoding.en_mresp{el};
        native_resp = sent_encoding.sp_mresp{el};
    else
        % English
        unf_resp = sent_encoding.sp_mresp{el};
        native_resp = sent_encoding.en_mresp{el};
    end
    % find the area under the curve for native

    % zero out everything that is negative
    native_resp(native_resp < 0) = 0;
    unf_resp(unf_resp < 0) = 0;
    
    native_auc = trapz(native_resp(50:200));
    unf_auc = trapz(unf_resp(50:200));

    % calculate selectivity index
    sent_encoding.selectivity(el) = (native_auc - unf_auc) / (native_auc + unf_auc);
end

% plot selectivity index
figure;
scatter(find(sent_encoding.ls == 1), ...
    sent_encoding.selectivity(sent_encoding.ls == 1), 20, 'filled'); hold on;
scatter(find(sent_encoding.ls == 2), ...
    sent_encoding.selectivity(sent_encoding.ls == 2), 20, 'filled'); hold on;
ylim([-1 1]);
h=yline(0);
h.LineWidth=2;

% use the same amount of spanish and english electrodes
sidx = find(sent_encoding.ls == 1);
eidx = find(sent_encoding.ls == 2);
sidx = sidx(randperm(length(sidx), length(eidx)));
idx = [sidx; eidx];
[h, p] = ttest(sent_encoding.selectivity(idx));
disp(['Selectivity index: ' num2str(h) ', p=' num2str(p)]);
figure; histogram(sent_encoding.selectivity(idx));



% check if this is significant using an lme
% lme = fitlme(sent_encoding, 'selectivity ~ ls + (1|SID)+(1|el:SID)');


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Visualize TIMIT vs. DIMEx mean response and maximum response time

longest = nan(1, height(sent_encoding));
for el = 1:height(sent_encoding)
    fvals = sent_encoding.fvals{el};
    meethresh = (fvals > sent_encoding.fthresh(el));
    contig = accumarray(nonzeros((cumsum(~meethresh) + 1) .* meethresh), 1);

    if ~isempty(contig)
        longest(el) = max(contig);
    end
end

% Find electrodes where the difference is longer than 200ms
idx = longest > 15;
cols = getColorsCrossComp(1);

figure;
subplot(2, 2, 1)
scatter(sent_encoding.maxtp(idx, 2) ./ 100 - bef, ...
    sent_encoding.maxtp(idx, 1) ./ 100 - bef, 25, ...
    sent_encoding.ls(idx), 'filled');
ylabel('TIMIT max timepoint');
xlabel('DIMEx max timepoint');
colormap(cols(1:2, :));
refline(1, 0);

subplot(2, 2, 2)
scatter(sent_encoding.maxresp(idx, 2), sent_encoding.maxresp(idx, 1), 25, ...
    sent_encoding.ls(idx), 'filled');
ylabel('TIMIT response');
xlabel('DIMEx response');
colormap(cols(1:2, :));
refline(1, 0);

% Calculate the differences and perform statistical tests
diffs = diff(sent_encoding.maxresp(idx, :)');
[p, ~] = ranksum(diffs(sent_encoding.ls(idx) == 1), ...
    diffs(sent_encoding.ls(idx) == 2));
disp(['Max response diff=' num2str(p)]);

subplot(2, 2, 4);
for l = 1:2
    histogram(diffs(sent_encoding.ls(idx) == l), ...
        EdgeColor='none', FaceColor=cols(l, :));
    hold on;
end
title('Response magnitude difference');

diffs = diff(sent_encoding.maxtp(idx, :)');
[p, ~] = ranksum(diffs(sent_encoding.ls(idx) == 1), ...
    diffs(sent_encoding.ls(idx) == 2));
disp(['Max tp diff=' num2str(p)])

subplot(2, 2, 3);
for l = 1:2
    histogram(diffs(sent_encoding.ls(idx) == l), ...
        EdgeColor='none', FaceColor=cols(l, :));
    hold on;
end
title('Maximum time point difference');
legend({'Spanish mono', 'English mono'});

% tmp = find(idx);
% elecs = tmp(abs(diffs) > 40);
% plotStitchedSentence(elecs, sent_encoding, sent_encoding.swind(1), ...
%     sent_encoding.ewind(1), bef);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Max peak amplitude and peak latency across subjects

longest = nan(1, height(sent_encoding));
for el = 1:height(sent_encoding)

    fvals = sent_encoding.fvals{el};
    meethresh = (fvals>sent_encoding.fthresh(el));
    contig = accumarray(nonzeros((cumsum(~meethresh)+1).*meethresh),1);

    if ~isempty(contig)
        longest(el)=max(contig);
    end
end

% find electrodes where difference is longer than 200ms
fidx = longest>-1; %20;

maxresp = sent_encoding.maxresp';
figure('Renderer', 'Painters');
ls = sent_encoding.ls;

cols = getColorsCrossComp(1);
maxls = 5;

subplot(2, 1, 1);
% response peak
for l = unique(ls)'
    idx =  ls==l & fidx';
    scatter3(maxresp(2,  idx), maxresp(1,  idx), find(idx), 10, ...
        cols(l, :), 'filled', 'MarkerFaceAlpha', 0.8); hold on;
    view(2);

    correl = corr(maxresp(2, idx)', maxresp(1, idx)', ...
    'type', 'Pearson', 'rows','complete');
    text(0, 3+l*0.2, ['r=' num2str(correl)], 'FontSize', 15);
end

% Formatting
xlabel({'Spanish peak', 'HFA'});
ylabel({'English peak',  'HFA'});
ylim([0, 3.9]);
xlim([0, 3.9]);
grid off

h=refline(1, 0);
h.Color = 'k';
h.LineWidth = 2;
h.HandleVisibility = "off";
colormap(cols);
set(gca, 'FontSize', 15);

subplot(2, 1, 2);
maxtp = (sent_encoding.maxtp')./100-bef;

% Scatter response peak latency
for l = unique(ls)'
    idx = ls==l & all(maxresp>0.75)' & fidx';
    scatter3(maxtp(2,  idx), maxtp(1,  idx), find(idx), mean(maxresp(:, idx))*10, ...
        cols(l, :), 'filled', 'MarkerFaceAlpha', 0.8); hold on;
    view(2);

    correl = corr(maxtp(2, idx)', maxtp(1, idx)', ...
    'type', 'Pearson', 'rows','complete');
    text(0, 1.5+l*0.1, ['r=' num2str(correl)], 'FontSize', 15);
end

% Formatting
xlabel({'Spanish peak latency'});
ylabel({'English peak latency'});
ylim([0.1, 2]);
xlim([0.1, 2]);
grid off

h=refline(1, 0);
h.Color = 'k';
h.LineWidth = 2;
h.HandleVisibility = "off";
colormap(cols);
set(gca, 'FontSize', 15);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Alternate visualization of peak amplitude comparison across subjects

figure;
ctr = 1;
titles = {'Spanish mono', 'English mono', 'Spanish-English bilingual'};

% Iterate over language types
for ls = [1, 2]
    idx = sent_encoding.ls == ls;
    peakresp = sent_encoding.maxresp;

    % Create subplot
    subplot(2, 3, ctr)

    % Marker type for selected from single language or both language
    type = all(sent_encoding.type,2);

    % Display the proportion of single language vs. both language
    disp(['Language type: ' num2str(ls)]);
    disp(['Both language: ' num2str(sum(idx & type)) ', ' ...
        num2str(sum(idx & type)/sum(idx)) '%']);
    disp(['Single languages: ' num2str(sum(idx & ~type)) ', ' ...
        num2str(sum(idx & ~type)/sum(idx)) '%']);
    
    % Scatter plot of peak HFA values
    scatter3(peakresp(idx & type, 2), peakresp(idx & type, 1), find(idx&type), ...
        25, [0.6 0.6 0.6], 'filled', 'Marker', 'o', 'MarkerFaceAlpha', 0.8); hold on;
    scatter3(peakresp(idx & ~type, 2), peakresp(idx & ~type, 1), find(idx&~type), ...
        25, [0 0 0], 'filled', 'Marker', '^', 'MarkerFaceAlpha', 0.6); hold on;
%     scatter3(peakresp(idx, 2), peakresp(idx, 1), find(idx), 15, [0.6 0.6 0.6], ...
%              'filled', 'Marker', 'o', 'MarkerFaceAlpha', 0.8); hold on;
    view(2);

    % Calculate correlation coefficient
    [rho, pval] = corr(peakresp(idx, 2), peakresp(idx, 1), 'rows', 'pairwise');
    text(0.5, 4, ['r=' num2str(rho, 2)], 'FontSize', 15);
    text(0.5, 3.5, ['p=' num2str(pval, 2)], 'FontSize', 15);

    ylabel({'English peak', 'HFA (z)'});
    xlabel({'Spanish peak', 'HFA (z)'});
    set(gca, 'FontSize', 15);
    ylim([0 4.5]);
    xlim([0 4.5]);
    title(titles{ctr});

    % Reference line
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    h.HandleVisibility = "off";
    %legend({'Both languages', 'Single language'}, 'Location', 'northwest');

    % Create histogram of peak HFA values
    subplot(2, 3, ctr+3);
    histogram(peakresp(idx, 2)-peakresp(idx, 1), EdgeColor='none', ...
        FaceColor = [0.4, 0.4, 0.4]); 
    xline(0);
    axis off;
    ctr = ctr + 1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Alternate visualization of peak latency comparison across subjects

figure;
ctr=1;
titles = {'Spanish mono', 'English mono', 'Spanish-English bilingual'};

% count number of significant time points difference
signum = arrayfun(@(x) sum(sent_encoding.fvals{x}>sent_encoding.fthresh(x)), ...
    1:height(sent_encoding));

% Iterate over language types
for ls = [1, 2]
    idx = sent_encoding.ls == ls & signum'>10;
    dataf = 100;
    peaktp = (sent_encoding.maxtp / dataf) - bef;

    % Create subplot
    subplot(2, 3, ctr)

    % Marker type for selected from single language or both language
    type = all(sent_encoding.type,2);

    % Display the proportion of single language vs. both language
    disp(['Language type: ' num2str(ls)]);
    disp(['Both language: ' num2str(sum(idx & type)) ', ' ...
        num2str(sum(idx & type)/sum(idx)) '%']);
    disp(['Single languages: ' num2str(sum(idx & ~type)) ', ' ...
        num2str(sum(idx & ~type)/sum(idx)) '%']);
    
    % Scatter plot of peak HFA values

    % size is the magnitude of the peak response
    sz = sent_encoding.maxresp(idx, 1) * 10;
    scatter3(peaktp(idx, 2), peaktp(idx, 1), find(idx), sz, [0.6 0.6 0.6], ...
             'filled', 'Marker', 'o', 'MarkerFaceAlpha', 0.8); hold on;
    view(2);

    % Calculate correlation coefficient
    [rho, pval] = corr(peaktp(idx, 2), peaktp(idx, 1), 'rows', 'pairwise');
    text(0.25, 1.9, ['r=' num2str(rho, 2)], 'FontSize', 15);
    text(0.25, 1.8, ['p=' num2str(pval, 2)], 'FontSize', 15);

    ylabel({'English peak time', 'from onset'});
    xlabel({'Spanish peak time', 'from onset'});
    xline(1, 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':');
    yline(1, 'LineWidth', 2, 'Color', 'k', 'LineStyle', ':');
    xline(1.5, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-');
    yline(1.5, 'LineWidth', 1.5, 'Color', 'k', 'LineStyle', '-');
    set(gca, 'FontSize', 15);
    title(titles{ctr});
    ylim([0 2]);
    xlim([0 2]);

    % Reference line
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    h.HandleVisibility = "off";
    %legend({'Both languages', 'Single language'}, 'Location', 'northwest');

    % Create histogram of peak HFA values
    subplot(2, 3, ctr+3);
    histogram(peaktp(idx, 2)-peaktp(idx, 1), EdgeColor='none', ...
        FaceColor = [0.4, 0.4, 0.4]); 
    xline(0);
    axis off;
    ctr = ctr + 1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% Compare within vs. across language single sentence response correlation

swind = 150;
ewind = 100;
d = nan(height(sent_encoding), 4);

% Calculate correlation metrics for each electrode
for el = 1:height(sent_encoding)
%     en_resp = sent_encoding.en_sent_resp{el};
%     sp_resp = sent_encoding.sp_sent_resp{el};

    en_resp = permute(mean(sent_encoding.en_rep_resp{el}, 2), [2, 1, 3]);
    sp_resp = permute(mean(sent_encoding.sp_rep_resp{el}, 2), [2, 1, 3]);
    
    if ~isempty(en_resp) && ~isempty(sp_resp)

        % Construct the onset and offset locked responses
        [en_resp, ~] = stitchedResp(en_resp, swind, ewind, 0);
        [sp_resp, ~] = stitchedResp(sp_resp, swind, ewind, 0);
    end

    if ~isempty(en_resp) && ~isempty(sp_resp)

        % Calculate single language correlation matrices
        en_corr = triu(corr(squeeze(en_resp), 'type', 'spearman'), 1);
        en_corr(en_corr == 0) = NaN;

        sp_corr = triu(corr(squeeze(sp_resp), 'type', 'spearman'), 1);
        sp_corr(sp_corr == 0) = NaN;

        d(el, 1) = mean([en_corr(:); sp_corr(:)], 'omitnan');

        % Calculate cross language correlation matrices
        crosslang_corr = corr(squeeze(en_resp), squeeze(sp_resp), 'type', 'spearman');
        crosslang_corr(crosslang_corr == 0) = NaN;
        d(el, 2) = median(crosslang_corr, [1, 2]);
        d(el, 3) = median(en_corr(:), 'omitnan');
        d(el, 4) = median(sp_corr(:), 'omitnan');
    end
end
%
% Plot the correlation metrics
figure;
ctr=1;
for ls = [1, 2]
    subplot(1, 3, ctr);
    idx = sent_encoding.ls==ls;
    scatter(d(idx, 1), d(idx, 2), 25, [0.6 0.6 0.6], ...
        'filled', 'MarkerFaceAlpha', 0.8);

    % Calculate correlation coefficient
    [rho, pval] = corr(d(idx, 1), d(idx, 2), 'rows', 'pairwise');
    text(0.1, 0.5, ['r=' num2str(rho, 2)], 'FontSize', 15);
    text(0.1, 0.47, ['p=' num2str(pval, 2)], 'FontSize', 15);

    ylim([0 0.6]);
    yticks([0 0.5]);
    xlim([0 0.6]);
    xticks([0 0.5]);

    % Reference line
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    
    xlabel({'Within language', 'correlation'});
    ylabel({'Across language', 'correlation'});
    
    set(gca, 'FontSize', 15);
    ctr=ctr+1;
end

% example electrode
figure; 
subplot(1, 3, 1); imagesc(en_corr); title('TIMIT'); caxis([0 0.6]);
subplot(1, 3, 2); imagesc(sp_corr); title('DIMEx'); caxis([0 0.6]);
subplot(1, 3, 3); imagesc(crosslang_corr); title('Cross'); caxis([0 0.6]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Plot stitched sentence response ERPs

% Asking the user for the electrode selection
% elec = input('Which electrodes?');

% Different sets of fixed electrodes
% elec = [986, 996, 946, 715, 711, 477, 237, 289, 294, 44]; 
% elec = [229, 155, 237, 142, 450];
% elec = [947, 550, 879, 712, 986, 477];
% elec = [1, 18, 52, 54, 55, 60, 77];
% elec = [698, 565, 61, 3 407, 416, 400, 636, 653]; % 
% elec = [416, 545, 400, 438]; % 
% elec = [338, 334, 223, 5]; % make threshold for one language and not the other
elec = [82, 135, 224];

% elec = find(abs(diff(sent_encoding.maxtp'))>50 ...
%     & all(sent_encoding.maxresp'>0.5) ...
%     & longest>20);

% elec = find(strcmp(sent_encoding.SID, 'EC252')' ...
%     & all(sent_encoding.maxresp'>0.75) ...
%     & ismember(sent_encoding.el, [21, 23, 52, 36])');

% find all electrodes where maxtp is very different and HG magnitude is
% large
longest = nan(1, height(sent_encoding));
for el = 1:height(sent_encoding)
    fvals = sent_encoding.fvals{el};
    meethresh = (fvals>sent_encoding.fthresh(el));
    contig = accumarray(nonzeros((cumsum(~meethresh)+1).*meethresh),1);

    if ~isempty(contig)
        longest(el)=max(contig);
    end
end

swind = sent_encoding.swind(1);
ewind = sent_encoding.ewind(1);

plotStitchedSentence(elec, sent_encoding, swind, ewind, bef);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Plotting Native brain comparing TIMIT and DIMEx response patterns (imagesc)

%'EC100', 'EC252', 'EC152', 'EC212', 'EC235', 'EC129', 'EC159', 'EC196'
nativeSIDs = {'EC214'};
% cm = spectral(20);
% cm = cm([2, 18], :);
% cm = [0.9 0.9 0.9; 0.3 0.3 0.3];

% Color by language
cm = [0 0 1; 1 0 0];
modelname={'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; 
modelfeatures  = [{'onset'}; timit_details.features.names([1:3, 8, 9, 11]); ... 
    {'peakrate'; 'F1'; 'F2'; 'F3'; 'F4'}];

bins = 15;
elecs = 203; % for EC183 , 156, 160
%elecs = 137; % for EC172
elecs = [3, 54]; % for EC214
% SID: [183, 214, 186, 195, 105]
% Elec: [72, 245, 212, 203, 63]
for l = 1:2
    binedges = linspace(prctile(sent_encoding.maxresp(:, l), 5), ...
        prctile(sent_encoding.maxresp(:, l), 100), bins);
    conds = discretize(sent_encoding.maxresp(:, l), binedges);
   
    % initialize design electrode structure
    desel=struct();
    desel.conds = unique(conds(~isnan(conds)))';
    desel.sz = (1:bins)*6;
    desel.cols = repmat(cm(l, :), bins, 1);    
    
    % split up peak rate and phonetic features again for MNI plotting
    desel.labels = [];

    % sent_encoding = sent_encoding(~strcmp(sent_encoding.SID, 'EC266'), :);
    for s=unique(sent_encoding.SID)'
        SID = s{1};
        idx = strcmp(sent_encoding.SID, SID) & sent_encoding.type(:, l)==1;
        desel.(SID).elid = sent_encoding.el(idx);
        desel.(SID).condition = conds(idx);

        desel.(SID).selid = elecs; %sent_encoding.el(idx);
    end
    [native] = plotNativeElec(nativeSIDs, desel, 1);
end

% to plot empty brain with just elec position
% desel.('EC100').elid = [];
% plotNativeElec({'EC100'}, desel, 1);

% elecs = [22, 150];
% Plot example electrodes
for sid = nativeSIDs
    SID = sid{1};   
%     elecs = find(strcmp(sent_encoding.SID, SID) & ...
%         all(sent_encoding.maxresp(:, :)>0.85, 2) & ...
%         diff(abs(sent_encoding.maxresp'))'<0.02,40)';

    els = find(strcmp(sent_encoding.SID, SID) ...
        & ismember(sent_encoding.el, elecs))';
    plotStitchedSentence(els, sent_encoding, 150, 100, 0.5, 0);
    yticks([-1 0 3]);

    for el = els
        % aggregate all sentences, find mean 
        en_resp = sent_encoding.en_sent_resp{el};
        sp_resp = sent_encoding.sp_sent_resp{el};
    
        % for each sentence use first 150ms + last 100ms 
        [en_resp, ~] = stitchedResp(en_resp, 150, 100);
        [sp_resp, ~] = stitchedResp(sp_resp, 150, 100);
        
        figure; 
        ax = subplot(2, 2, 1);
        min_sent = min([size(en_resp, 3), size(sp_resp, 3), 100]);
        imagesc(smoothdata(squeeze(sp_resp(1, :, 1:min_sent))', 2)); 
        xline(150, 'LineWidth', 2); 
        xline(50, 'LineWidth', 2, 'LineStyle', '--');
        xline(200, 'LineWidth', 2, 'LineStyle', '--');
        caxis([0 3]);
        colormap(ax, flipud(blues));
        

        ax = subplot(2, 2, 3);
        imagesc(smoothdata(squeeze(en_resp(1, :, 1:min_sent))', 2)); 
        xline(150, 'LineWidth', 2); 
        xline(50, 'LineWidth', 2, 'LineStyle', '--');
        xline(200, 'LineWidth', 2, 'LineStyle', '--');
        colormap(ax, flipud(reds));
        caxis([0 3]);
        sgtitle([SID ': ' num2str(sent_encoding.el(el)) ', ' num2str(el)])

        % imagesc the encoding models
        ax = subplot(2, 2, 2);
        [weights] = getTRFweights(SID, sent_encoding.el(el), ...
            'dimex', modelname, datapath);
        imagesc(weights(:, 1:40));
        colormap(ax, inferno);
        caxis([-1.5 1.5]);
        title('dimex');

        yticks(1:length(modelfeatures));
        yticklabels(modelfeatures);
        xlim([0.5 40])
        xticks([1 40]);
        xticklabels({'0', '-0.4'});
        xlabel('Time (s)');
        set(gca, 'FontSize', 13);
        clear strf

        ax = subplot(2, 2, 4);
        [weights] = getTRFweights(SID, sent_encoding.el(el), ...
            'timit', modelname, datapath);
        imagesc(weights(:, 1:40));
        colormap(ax, inferno);
        caxis([-1.5 1.5]);
        title('timit');

        yticks(1:length(modelfeatures));
        yticklabels(modelfeatures);
        xlim([0.5 40])
        xticks([1 40]);
        xticklabels({'0', '-0.4'});
        xlabel('Time (s)');
        set(gca, 'FontSize', 13);
        clear strf
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% Plotting MNI brain comparing TIMIT and DIMEx response patterns

% Plotting on MNI Brain the peak response difference between corpora
maxresp = sent_encoding.maxresp';
binedges = -0.20:0.02:0.20;
colors = brewermap(length(binedges)-1, 'Spectral');

nanidx = any(isnan(maxresp));
peak_diffs = nan(1, size(maxresp, 2));
peak_diffs(~nanidx) = diff(maxresp(:, ~nanidx)) / max(maxresp(~nanidx));

% Negative values indicate English greater, positive values indicate Spanish greater
conds = discretize(peak_diffs, binedges);
sids = cellfun(@(x) str2double(x(3:end)), [sent_encoding.SID]);

figure;
scatter3(1:length(peak_diffs), peak_diffs, sids, 45, sent_encoding.ls, ...
    'filled', 'MarkerFaceAlpha', 0.5);
hold on;
cm = getColorsCrossComp(1);
colormap(cm(1:4, :));
view(2);
yline(0, 'LineWidth', 2);

% Initialize design electrode structure
desel = struct();
desel.conds = unique(conds);
desel.sz = repmat(50, 1, length(desel.conds));
desel.cols = colors;

% Split up peak rate and phonetic features again for MNI plotting
desel.labels = [];
for s = unique(sent_encoding.SID)'
    SID = s{1};
    idx = strcmp(sent_encoding.SID, SID);
    desel.(SID).elid = sent_encoding.el(idx);
    desel.(SID).condition = conds(idx);
    desel.(SID).selid = []; %sent_encoding.el(idx);
end

% Plot histogram of count across y-axis (anterior to posterior)
[~] = plotMNIElec(unique(sent_encoding.SID), desel, 'lh', 0);
[~] = plotMNIElec(unique(sent_encoding.SID), desel, 'rh', 0);

nativeSIDs = {'EC252'};
for i = 1:length(nativeSIDs)
    desel.(nativeSIDs{i}).selid = 52;
end
[~] = plotNativeElec(nativeSIDs, desel, 1);

% Look at specific subjects
for s = nativeSIDs
    SID = s{1};
    figure;
    idx = strcmp(sent_encoding.SID, SID);
    scatter(sent_encoding.maxresp(idx, 2), sent_encoding.maxresp(idx, 1), ...
        25, desel.(SID).condition, 'filled', 'MarkerEdgeColor', 'k');
    hold on;
    colormap(desel.cols);
    title(SID);
    caxis([1 length(binedges)-1]);

    % Formatting
    xlabel({'Spanish peak', 'HFA (z)'});
    ylabel({'English peak', 'HFA (z)'});
    ylim([0, 3.9]);
    xlim([0, 3.9]);
    grid off;

    h = refline(1, 0);
    h.Color = 'k';
    h.LineWidth = 2;
    h.HandleVisibility = 'off';
    colormap(desel.cols);
    set(gca, 'FontSize', 15);

    % Plot the ERPs for this subject between Spanish and English
    resp = nan(450, 180, 2);
    en = cat(1, sent_encoding.en_sent_resp{idx});
    resp(sent_encoding.el(idx), :, 1) = mean(en(:, 1:180, :), 3, 'omitnan');
    sp = cat(1, sent_encoding.sp_sent_resp{idx});
    resp(sent_encoding.el(idx), :, 2) = mean(sp(:, 1:180, :), 3, 'omitnan');

    desmat = struct();
    desmat.condition = [1, 2];
    desmat.names = {'English', 'Spanish'};
end

% Legend
figure;
subs = [1:5:binedges(end)*100 binedges(end)*100];
scatter(1:length(subs), ones(1, length(subs)), ...
    135, desel.cols(subs, :), 'filled');
text((1:length(subs))-0.3, ones(1, length(subs))*0.85, ...
    split(num2str((binedges(subs)+0.01)*100)), 'FontSize', 15);
text(1-0.65, 0.65, 'En > Sp (%)', 'FontSize', 15);
text(5-0.65, 0.65, 'Sp > En (%)', 'FontSize', 15);
xlim([0 6]);
axis off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *D* maxresp;

%% plot MNI brains for Spanish only and English only electrodes

% Initialize design electrode structure
desel = struct();
desel.conds = [0, 1, 2];
desel.sz = [25, 20, 20];

binedges = -0.20:0.02:0.20;
colors = flipud(brewermap(length(binedges)-1, 'Spectral'));
desel.cols = [colors(round(size(colors, 1)/2), :); colors(end-2, :); colors(3, :)];

conds = zeros(1, height(sent_encoding));
for i = 1:height(sent_encoding)
    if sent_encoding.type(i, :)==[1, 0] 
        conds(i) = 1;
    elseif sent_encoding.type(i, :)==[0, 1] 
        conds(i) = 2;
    end
end

% Split up peak rate and phonetic features again for MNI plotting
desel.labels = [];
for s = unique(sent_encoding.SID)'
    SID = s{1};
    idx = strcmp(sent_encoding.SID, SID);
    desel.(SID).elid = sent_encoding.el(idx);
    desel.(SID).condition = conds(idx);
end

% Plot histogram of count across y-axis (anterior to posterior)
[mni_lh] = plotMNIElec(eSIDs, desel, 'lh', 0);
% light("Style","infinite","Position",[-50 -50 0]);

[mni_rh] = plotMNIElec(eSIDs, desel, 'rh', 0);
% light("Style","infinite","Position",[50 50 0]);

%% plot MNI for native electrodes and nonnative electrodes

%'EC100', 'EC252', 'EC152', 'EC212', 'EC235', 'EC129', 'EC159', 'EC196'

native = 1;

% Color by language
cm = [1 0 0; 0 0 1];

bins = 15;
binedges = linspace(0.5, 4.5, bins);
 % initialize design electrode structure
desel=struct();
desel.conds = 1:bins;
desel.sz = (1:bins)*6;

desel.cols = repmat(cm(native+1, :), bins, 1);   
% split up peak rate and phonetic features again for MNI plotting
desel.labels = [];

for ls = 3
    for s=unique(sent_encoding.SID(sent_encoding.ls==ls))'    
        SID = s{1};

        if native % timit for bilinguals
            % first row is english, so ls==1 (Spanish), native language is
            % row 2
            c =  1;% mod(ls, 2)+1;               
        else % spanish for bilinguals
            % first row is english, so ls==1 (Spanish), nonnative language is
            % row 1
            c = 2; %ls;
        end

        idx = strcmp(sent_encoding.SID, SID) & sent_encoding.type(:, c)==1;
        conds = discretize(sent_encoding.maxresp(idx, c), binedges);
        
        desel.(SID).elid = sent_encoding.el(idx);
        desel.(SID).condition = conds;
    end
end

[mni_lh] = plotMNIElec(SIDs, desel, 'lh', 0);
light("Style","infinite","Position",[-50 -50 0]);
alpha 0.5

[mni_rh] = plotMNIElec(SIDs, desel, 'rh', 0);
light("Style","infinite","Position",[50 -50 0]);
alpha 0.5

% to plot empty brain with just elec position

%% Compare Spanish / English speech responsive electrode sets

type = sent_encoding.type';
maxresp = sent_encoding.maxresp';
binedges = -0.20:0.02:0.20;
colors = flipud(brewermap(length(binedges)-1, 'Spectral'));
cols = [colors(3, :); colors(end-2, :); colors(round(size(colors, 1)/2), :)];

% Plotting the first two subplots
figure('Renderer', 'Painters'); 

subplot(1, 3, 1);
nnanidx = find(~any(isnan(type)));
both = nnanidx(all(type(:, nnanidx)));

scatter3(maxresp(2,  both), maxresp(1,  both), both, 40, cols(3, :), ...
    'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k');
view(2);
legend({'Shared response', 'DIMEx response'});

subplot(1, 3, 2);
dimexcl = nnanidx(~all(type(:, nnanidx)) & type(1, nnanidx));
timexcl = nnanidx(~all(type(:, nnanidx)) & type(2, nnanidx));

scatter3(maxresp(2,  dimexcl), maxresp(1,  dimexcl), find(dimexcl), 40, cols(1, :), ...
    'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k');
hold on;
scatter3(maxresp(2,  timexcl), maxresp(1,  timexcl), find(timexcl), 40, cols(2, :), ...
    'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'k');
view(2);
legend({'DIMEx response', 'TIMIT response'});

% Formatting for both plots
for i = 1:2
    subplot(1, 3, i);
    grid off;
    view(2);
    set(gca, 'FontSize', 15);
    ylim([0 6]);
    xlim([0 6]);
    h = refline(1, 0);
    h.Color = 'k';
    h.LineWidth = 2;
    h.HandleVisibility = 'off';
    colormap(cols);
    lgn = legend();
    lgn.Title.String = 'Electrode type';
    xlabel({'Spanish peak', 'HGA (z-score)'});
    ylabel({'English peak',  'HGA (z-score)'});
end

% Plotting the third subplot
subplot(1, 3, 3);
boxchart(mean([maxresp(2,  both); maxresp(1,  both)])', ...
    'BoxFaceColor', cols(3, :), 'XData', ones(1, length(both)), ...
    'MarkerStyle', '+', 'MarkerColor', cols(3, :), 'BoxLineColor', 'k');
hold on;
boxchart(mean([maxresp(2,  timexcl); maxresp(1, timexcl)])', ...
    'BoxFaceColor', cols(2, :), 'XData', 2 * ones(1, length(timexcl)), ...
    'MarkerStyle', '+', 'MarkerColor', cols(2, :), 'BoxLineColor', 'k');
boxchart(mean([maxresp(2,  dimexcl); maxresp(1, dimexcl)])', ...
    'BoxFaceColor', cols(1, :), 'XData', 3 * ones(1, length(dimexcl)), ...
    'MarkerStyle', '+', 'MarkerColor', cols(1, :), 'BoxLineColor', 'k');

% Formatting for the third subplot
ylim([-0.07 4]);
ylabel('mean HGA (z-score)');
set(gca, 'FontSize', 15, 'Xtick', [], 'Ytick', 1:4);            

% Separate figure looking at the single DIMEx or TIMIT electrodes to see
% whether are found more in native subjects of the language
figure;
colrs = getColorsCrossComp(1);

subplot(2, 2, 1);
scatter(maxresp(2,  timexcl), maxresp(1,  timexcl), 40, sent_encoding.ls(timexcl), ...
    'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
xlabel({'Spanish peak', 'HGA (z-score)'});
ylabel({'English peak',  'HGA (z-score)'});
title('TIMIT');
refline(1, 0);

subplot(2, 2, 2);
scatter(maxresp(2,  dimexcl), maxresp(1,  dimexcl), 40, sent_encoding.ls(dimexcl), ...
    'filled', 'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
xlabel({'Spanish peak', 'HGA (z-score)'});
ylabel({'English peak',  'HGA (z-score)'});
title('DIMEx');
colormap(colrs(1:2, :));
refline(1, 0);

subplot(2, 2, 3);
for ls = 1:2
    idx = intersect(timexcl, find(sent_encoding.ls==ls));
    histogram(maxresp(2,  idx)-maxresp(1, idx), ...
        'EdgeColor', 'none', 'FaceColor', colrs(ls, :)); 
    ylabel('Difference between peak HFA');
    colormap(colrs(1:2, :));
    xline(0, 'LineWidth', 2, 'Color', 'k');
    title('TIMIT responsive only');
    hold on;
end

subplot(2, 2, 4);
for ls = 1:2
    idx = intersect(dimexcl, find(sent_encoding.ls==ls));
    histogram(maxresp(2,  idx)-maxresp(1, idx), ...
        'EdgeColor', 'none', 'FaceColor', colrs(ls, :)); 
    ylabel('Difference between peak HFA');
    title('DIMEx responsive only');
    xline(0, 'LineWidth', 2, 'Color', 'k');
    hold on;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% Compare average sentence response between corpora across time (windowed)

tempresp = permute(sent_encoding.tempresp, [2, 3, 1]);
wndsz = 5;
cols = getColorsCrossComp(1);
sliding = 65:wndsz*2:180;
meandiff = nan(4, length(sliding));
stdiff = meandiff;

% Supplementary for all time points
figure;
for i = 1:length(sliding)
    subplot(4, length(sliding), i);
    scatter(squeeze(tempresp(2, i, :)), squeeze(tempresp(1, i, :)), 34, ...
        sent_encoding.ls, 'filled', 'MarkerFaceAlpha', 0.5);
    colormap(brewermap(4, 'Dark2'));
    box off;
    ylim([-2 4]);
    xlim([-2 4]);
    if i == 1
        ylabel('English HFA (z)');
        xlabel('Spanish HFA (z)');
    end

    yticks([0, 4]);
    xticks([0, 4]);
    title([num2str((sliding(i)-50)./100) 's center']);
    h = refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';

    ctr = 1;
    for ls = [1, 2, 4]
        subplot(4, length(sliding), i+ctr*length(sliding));
        proportdiff = diff(squeeze(tempresp(:, i, sent_encoding.ls==ls)))./...
            max(abs(squeeze(tempresp(:, i, sent_encoding.ls==ls))));
        histogram(proportdiff, 30, 'FaceColor', cols(ls, :), 'EdgeColor', 'none');
        hold on;
        ylim([0 150]);
        yticks([]);
        box off;
        [f, x] = ksdensity(proportdiff);
        yyaxis right;
        plot(x, f, 'LineWidth', 1.5, 'color', cols(ls, :));
        xline(0, 'LineWidth', 3);
        xlim([-1 1]);
        yticks([]);

        meandiff(ls, i) = mean(proportdiff, 'omitnan');
        stdiff(ls, i) = std(proportdiff, 'omitnan');

        ctr = ctr + 1;
    end
end

figure;
for ls = [1, 2, 4]
    errorbar((sliding-50)./100, meandiff(ls, :), stdiff(ls, :), '-o', ...
        'LineWidth', 1.75, 'Color', cols(ls, :));
    hold on;
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
end
ylabel('mean proportion difference');
xlabel('time following onset (s)');
set(gca, 'FontSize', 15);


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*


%% Subtraction method

fvals = cell2mat([sent_encoding.fvals(:)]);
fthresh = sent_encoding.fthresh(1);
fvals(fvals<fthresh) = NaN;
nanidx = cellfun(@(x) isempty(x), [sent_encoding.fvals]);
prop = nan(2, 2);

nodiffelecs = all(isnan(fvals(:, 150:250)), 2);
titles = {'Spanish mono', 'English mono', '', 'Bilingual'};

english = cat(1, sent_encoding.en_mresp{~nanidx});
spanish = cat(1, sent_encoding.sp_mresp{~nanidx});

% so english preference = positive, spanish pref = negative
difference = english-spanish;
maxresp = max([cellfun(@(x) max(x), sent_encoding.en_mresp), ...
    cellfun(@(x) max(x), sent_encoding.sp_mresp)], [], 2);

figure; 
ctr=1;
for startp = [1, 150]
    for ls = [1, 2, 4]
        % plotting average difference, y-axis is number of sig points
        x = median(difference(:, startp:min(250, startp+150)), 2);
        y = sum(fvals(:, startp:min(250, startp+150))>fthresh, 2);

        % only look at electrodes for which more than 10 time points are
        % different
        idx = sent_encoding.ls(~nanidx)==ls & y>1;
        prop(ls, double(ctr>2)+1) = sum(x(idx)>0)/length(x(idx));
        sz = maxresp*10;

        subplot(2, 3, ctr);    
        yyaxis right
        histogram(x(sent_encoding.ls(~nanidx)==ls & y>0), ...
            'FaceColor', [0.7 0.7 0.7], 'FaceAlpha', 0.4, ...
            'Edgecolor', 'none');
        xline(prctile(x(sent_encoding.ls(~nanidx)==ls & y>0), 50), ...
            'LineWidth', 1, 'Color', 'r');
        ylim([0 100]);
        yticks([]);
        box off;

        yyaxis left
        scatter3(x(idx), y(idx), find(idx), sz(idx), 'filled', 'k', ...
            'MarkerFaceAlpha', 0.6);
        view(2);
        title(titles{ls});
        xlim([-0.4 0.4]);
        xticks([-0.4 0.4]);
        xlabel('median HFA difference')
        xticklabels({'Spanish pref', 'English pref'});
        xline(0, 'LineWidth', 2);
        yticks([0 100]);
        ylim([0 100]);
        ylabel('# timepoints sig diff');
        set(gca, 'FontSize', 13);

        ctr=ctr+1;
    end
end



%% Compare native and unfamiliar contour area

figure;
set(gcf,'Color','w');

ha = axes;
hemi = 'rh';
hold on

imgall = load_allimgdata;

types = [1];
numsid = 0;      

for native = types

    if strcmp(hemi, 'lh')
        cortex = imgall.(SIDs{1}).img_mni.cortex;
    else
        cortex = imgall.(SIDs{4}).img_mni.cortex;
    end

    PlotBrainSurface(cortex, hemi,'lateral');
    alpha 0.9
    %light("Style","infinite","Position",[100 100 0]);

    ax_pos = get(ha,'Position');
    yy_range = get(ha,'YLim');
    zz_range = get(ha,'ZLim');
    
    % find density map of native speech, lateral side
    native_xyz = [];
    native_hga = [];
    all_xyz = [];

    for si = SIDs
        sid = si{1};
        if strcmpi(imgall.(sid).hemi,hemi)
            subject = imgall.(sid);
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix;
    
    %             anatomy = subject.img_mni.anatomy;
    %             ch_lateral = find(contains(anatomy(:,4),{'central','temporal', ...
    %                 'parietal','supramarginal','frontal','pars'}));
    %             ch_temporal = find(contains(anatomy(:,4),{'temporal'}));
    %             
                ch_sid = sent_encoding.el(strcmp(sent_encoding.SID, sid));
                ls = sent_encoding.ls(find(strcmp(sent_encoding.SID, sid), 1, 'first'));
        
                if native % spanish for bilinguals
                    if ismember(ls, [1, 2])
                        ch_sel = sent_encoding.el(sent_encoding.type(:, ls)==1);
                    else % bilingual case
                        ch_sel = sent_encoding.el(sent_encoding.type(:, 1)==1);
                    end
    
                    % for weighting
                    if ls==1
                        native_hga = [native_hga; sent_encoding.maxresp(intersect(ch_sel, ch_sid), 1)];
                    elseif ls==2
                        native_hga = [native_hga; sent_encoding.maxresp(intersect(ch_sel, ch_sid), 2)];
                    end
                else
                    if ismember(ls, [1, 2])
                        ch_sel = sent_encoding.el(sent_encoding.type(:, mod(ls, 2)+1)==1);
                    else % bilingual case
                        ch_sel = sent_encoding.el(sent_encoding.type(:, 2)==1);
                    end
    
                    % for weighting
                    if ls == 1
                        native_hga = [native_hga; sent_encoding.maxresp(intersect(ch_sel, ch_sid), 2)];
                    elseif ls==2
                        native_hga = [native_hga; sent_encoding.maxresp(intersect(ch_sel, ch_sid), 1)];
                    end
                end
    
        %         ch_use = setdiff(ch_use,ch_temporal);
                native_xyz = [native_xyz; elecmatrix(intersect(ch_sel, ch_sid),:)]; 
                all_xyz = [all_xyz; elecmatrix]; 
                numsid=numsid+1;
            end
        end
    
    end
    
    yye = min(native_xyz(:,2))-15:1:max(native_xyz(:,2)+15);
    zze = min(native_xyz(:,3))-15:1:max(native_xyz(:,3)+15); 
    ds = histcounts2(native_xyz(:,2),native_xyz(:,3),yye,zze);
    ds_all = histcounts2(all_xyz(:,2),all_xyz(:,3),yye,zze);
    ds_norm = ds./ds_all;
    ds_norm(isnan(ds_norm)) = 0;
    
%     ds = hist2w(native_xyz(:,[2, 3]),native_hga,yye,zze);
%     ds(1, :) = [];
%     ds(:, 1) = [];
%     ds(isnan(ds))=0;

    gs_kernel = fspecial('gaussian', [15, 15], 3);
    ds_sm = conv2(ds_norm,gs_kernel,'same');
    
    yy = yye(1)+diff(yye(1:2))/2:diff(yye(1:2)):yye(end);
    zz = zze(1)+diff(zze(1:2))/2:diff(zze(1:2)):zze(end);
    yyq = yy(1):0.2:yy(end);
    zzq = zz(1):0.2:zz(end); 
    ds_q = interp2(zz,yy,ds_sm,zzq,yyq','cubic');
    
    ha_ct = axes('Position',ax_pos);
    hold on
    [~,hc] = contourf(yyq,zzq,ds_q',15,'LineColor','none');
    % [~,hc] = contourf(yye,zze,ds_sm',15,'LineColor','none');
    
    %caxis([-max(max(ds_q))*0.6,max(max(ds_q))])
    %caxis([0 1])
    caxis([0 0.18])
    if strcmp(hemi, 'lh')
        set(ha_ct,'Color','None','XDir','reverse');
    else
        set(ha_ct,'Color','None');
    end
    axis equal
    axis off
    xlim(yy_range);
    ylim(zz_range);
    %colormap(ha_ct,ds_cmap{ai});
    if ~native
        colormap(ha_ct, flipud(reds));
    else
        colormap(ha_ct, flipud(blues));
    end
    drawnow;
    
    
    pause(0.1); % keep this here to be able to update transparency
    % input('Press any key...')   % keep this here to be able to update transparency
    hFills = hc.FacePrims;  % array of TriangleStrip objects
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    for i = 1:length(hFills)
        hFills(i).ColorData(4) = 180;   % default=255
    end
    hFills(1).ColorData(4) = 0;
    disp(numsid);
end

figure;
if ~native
    colormap(flipud(reds));
else
    colormap(flipud(blues));
end
caxis([0 0.2]);
axis off;
colorbar();

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% ------------------- scratch --------------------------------------------

% itzik analysis: more variability in the response when you DO NOT know the
% language 
figure;
lang_fields = {'en_rep_resp', 'sp_rep_resp'};
sents = [10, 8];
min_rep = [6, 6];
titles = {'TIMIT', 'DIMEx'};
EV = cell(1, 2);
for lang = 1:2
    field = lang_fields{lang};
    subplot(1, 2, lang)

    % all the vars
    x = nan(sents(lang), 2);
    data = nan(sents(lang), 2);
    err = nan(sents(lang), 2);
    sent_length = nan(sents(lang), 1);

    for sent = 1:sents(lang)
        pd = nan(height(sent_encoding), 1);
        for i = 1:height(sent_encoding)
            rep_resp = sent_encoding.(field){i};
            if ~isempty(rep_resp) && size(rep_resp, 2)>min_rep(lang)

                % controlling for number of trials/repeats
                % time points x repeats x sentence
                sent_tmp = find(isnan(rep_resp(:, 1, sent)), 1, 'first')-1;
                if ~isempty(sent_tmp)
                    sent_length(sent) = sent_tmp;

                    % euclidean distance between all pairs of repeats
                    pd(i) = mean(pdist(rep_resp(1:sent_length(sent), ...
                        1:min_rep(lang), sent)'));    
                else 
                    sent_length(sent)=size(rep_resp, 1);
                    pd(i) = mean(pdist(rep_resp(:, ...
                        1:min_rep(lang), sent)'));
                end
            end
        end
        x(sent, :) = [2*sent-0.25 2*sent+0.25];
        data(sent, :) = [mean(pd(sent_encoding.ls==1), 'omitnan') ...
            mean(pd(sent_encoding.ls==2), 'omitnan')];
        err(sent, :) = [nansem(pd(sent_encoding.ls==1)) nansem(pd(sent_encoding.ls==2))];
    end

    % compute explainable variance for all sentences across Spanish and English monolinguals
    % replace all 0s with nans
    represp = cellfun(@(x) replaceMat(reshape(permute(x, [1, 3, 2]), [], size(x, 2)), 0, NaN), ...
                sent_encoding.(field), 'UniformOutput', false);

    % only use cells where there are at least 10/8 repeats
    idx = cellfun(@(x) size(x, 2), represp)>min_rep(lang);

    % also remove EC219 because it has a lot of nans
    idx = idx & ~strcmp(sent_encoding.SID, 'EC219');

    mintpsent = min(cellfun(@(x) size(x, 1), represp(idx)));
    idx = idx & cellfun(@(x) size(x, 1), represp)==mintpsent;
    tmp = arrayfun(@(x) represp{x}(1:mintpsent, 1:min_rep(lang)), find(idx), ...
        'UniformOutput', false);
    represp = cat(3, tmp{:});

    % remove all nans
    nanidx = all(isnan(represp(:, :, 1)), 2);
    represp = represp(~nanidx, :, :);

    % compute explainable variance
    EV{lang} = [squeeze(explainable_variance(represp)), find(idx)];
  
    % sort by sentence length:
    [sort_length, sort_idx] = sort(sent_length);

    bar(x(sort_idx, :), data(sort_idx, :), 3.5);                    
    hold on;
    for ls = 1:2
        er = errorbar(x(sort_idx, ls), data(sort_idx, ls),err(sort_idx, ls));    
        er.Color = [0 0 0];                            
        er.LineStyle = 'none';  
        er.HandleVisibility = 'off';
    end
    text(2*(1:sents(lang)), ...
        repmat(30+(lang-1)*10, sents(lang), 1), num2str(sent_length(sort_idx)/100))  
    
    title(titles{lang});
    legend({'Spanish mono', 'English mono'});
    box off;
end

interidx = intersect(EV{1}(:, 2), EV{2}(:, 2));
figure;
x = arrayfun(@(x) EV{1}(EV{1}(:, 2)==x, 1), interidx);
y = arrayfun(@(x) EV{2}(EV{2}(:, 2)==x, 1), interidx);
lss = sent_encoding.ls(interidx);

colors = getColorsCrossComp(1);
for ls = 1:2
    scatter3(x(lss==ls), y(lss==ls), interidx(lss==ls), 45, colors(ls, :), 'filled', ...
        'MarkerFaceAlpha', 0.6); hold on;
end
view(2)
xlabel('TIMIT ev'); ylabel('DIMEx ev');
legend({'Spanish mono', 'English mono'});
set(gca, 'FontSize', 13);
ylim([0 0.7]);
xlim([0 0.7]);
h=refline(1, 0);
h.Color = 'k';
h.LineWidth = 1.5;
h.HandleVisibility = 'off';


%% plotting
ploss = {{'g', 'k'}, {'b', 'p'}, {'d', 't'}};
figure;
ctr=1;
for p = ploss
    plos = p{1};
    subplot(1, 3, ctr)
    voiced = allplos(1, :) == vidx(strcmp(vs, plos{1}));
    voiceless = allplos(1, :) == vidx(strcmp(vs, plos{2}));
    histogram(allplos(2, voiced), 'Normalization', 'probability'); hold on
    histogram(allplos(2, voiceless),'Normalization', 'probability');
    legend(plos);
    ylabel('Probability');
    xlim([-0.45 11]);
    xticks([0 5 10]);
    ylim([0 0.45]);
    xticklabels({'0', '50', '100'});
    xlabel('VOT (ms)'); 
    set(gca, 'FontSize', 15);
    title(['v: ' num2str(sum(voiced)) ', nv: ' num2str(sum(voiceless))])
    ctr = ctr+1;
end

% combining all voiced and voiceless
ploss = {{'g', 'b', 'd'},{'k', 'p', 't'}};
cols = [0.2 0.2 0.2; 0.7 0.7 0.7];
figure;
idx = cell(2, 1);
for ctr=1:2
    plos = ploss{ctr};    
    idx{ctr} = sum(allplos(1, :) == vidx(ismember(vs, plos)));
    histogram(allplos(2, logical(idx{ctr})), 'Normalization', ...
        'probability', 'EdgeColor', 'k', 'FaceColor', cols(ctr, :)); hold on
end
title(['v: ' num2str(sum(idx{1})) ', nv: ' num2str(sum(idx{2}))])
legend({'voiced', 'voiceless'});
ylabel('Probability'); 

xlim([-0.45 11]);
xticks([0 5 10]);
ylim([0 0.45]);
xticklabels({'0', '50', '100'});
xlabel('VOT (ms)'); 
set(gca, 'FontSize', 15);


%% functions


function [weights] = getTRFweights(SID, el, corpus, modelname, datapath)    
    [strf] = loadMultModelStrf(SID, modelname, corpus, datapath, 1);  
    weights = strf{1}.meanStrf(:, :, el);
end
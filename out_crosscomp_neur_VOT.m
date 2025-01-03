%% Set up

% Ilina Bhaya-Grossman
% 01.08.2022
% addpath(genpath('../../../ecog_scripts'))
% addpath(genpath('../../../plotting_scripts'))
% addpath(genpath('loadData'))
% 
% addpath brewer
% addpath(genpath('util'))
% zFolder = 'block_z'; % 'block_z'
% [datapath, dpath] = setDatapath;
% addpath(genpath(datapath))
% 
% bef=20;
% aft=50;
% 
% % Note - EC202 has no STG coverage
% [sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
% sSIDs = [sSIDs, {'EC225', 'EC203'}];
% eSIDs = [eSIDs, {'EC242', 'EC142', 'EC219'}]; %, 'EC142'
% SIDs = [sSIDs, eSIDs]; %, {'HS11', 'HS9', 'HS10'}

%colors = getColorsCrossComp(4);

% timit_details = load('out_sentence_details_timit_all_loudness.mat');
% dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
% tps = 50:55;

out_crosscomp_startup;

plos = timit_details.features.plosive;
TDcons = loadDDcons('timit', 20, 50, [eSIDs sSIDs bSIDs], [], [], ...
    timit_details.sentdet, plos);

Dcons = loadDDcons('dimex', 20, 50, [], [], [], ...
    dimex_details.sentdet, plos);

% remove all cases in which stop-consonant is not pre-vocalic
TDcons = prevocCons(TDcons, timit_details.features);
Dcons = prevocCons(Dcons, dimex_details.features);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons*;

%% -------------------- VOT distribution in timit -------------------------

% set up
struct = timit_details;
closstr = 'cl'; %'_c'; % "_c" for dimex, "cl" for timit

features = struct.features;
sentdet = struct.sentdet;
phnnames = struct.phnnames;

% closure idx
allplos = [];
cls = features.plosive(contains(features.plosive, closstr));
clidx = find(ismember(phnnames, cls));
vidx = find(ismember(phnnames, cellfun(@(x) x(1), cls, 'UniformOutput', false)));
vs = phnnames(ismember(phnnames, cellfun(@(x) x(1), cls, 'UniformOutput', false)));
for i = 1:length(sentdet)
    % find all plosives
    [phntype, tps] = find(sentdet(i).phnmatonset);
    if isfield(phnnames, sentdet(i))
        sentphn = sentdet(i).phnnames;
    else
        sentphn = arrayfun(@(x) phnnames(x), phntype, 'UniformOutput', false);
    end

    disp(['sent ' num2str(i) ': ' strjoin([sentphn{:}])]);
    for c = find(ismember(phntype, clidx))' 
       if c+1<length(sentphn) && ismember(phntype(c+1), vidx)
           % calculate time delay between plosive (eg tens of ms)
           dur = tps(c+2)-tps(c+1);
           % includes phntype, duration of VOT, sentence, starting tp
           allplos = [allplos [phntype(c+1); dur; i; tps(c)]];
%            % includes phntype, duration of VOT, sentence, starting tp
%            allplos = [allplos [phntype(c-1); dur; i; tps(c-1)]];
       end
    end 
end

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
cols = [0.8 0.1 0.1; 0.1 0.2 0.8];
figure;
idx = cell(2, 1);
for ctr=1:2
    plos = ploss{ctr};    
    idx{ctr} = sum(allplos(1, :) == vidx(ismember(vs, plos)));
    histogram(allplos(2, logical(idx{ctr})), 'Normalization', ...
        'probability', 'FaceColor', cols(ctr, :), ...
        'EdgeColor','none'); hold on
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

figure;
histogram(allplos(2, allplos(1, :)==10), ...
    'FaceAlpha',0.35, 'FaceColor','r', 'EdgeColor','none'); hold on;
histogram(allplos(2, allplos(1, :)==43),  ...
    'FaceAlpha',0.35, 'FaceColor','b', 'EdgeColor','none'); hold on;
xlim([0.5 7.5]);
% yticks([]);
xticks(1:7);

% yticks([]);
ylabel('Count');
% xticklabels(split(num2str(-15:10:60)));
xlabel('VOT (ms)');
set(gca, 'FontSize', 15);
box off;

yyaxis right
[y, xi] = ksdensity(allplos(2, allplos(1, :)==10));
plot(xi, y, 'LineWidth', 3, 'Color','r');
[y, xi] = ksdensity(allplos(2, allplos(1, :)==43));
plot(xi, y, 'LineWidth', 3, 'Color','b', 'LineStyle','-');


figure;
plotVOTFreq(TDcons, {{'b'}, {'p'}},0)


%% matrix with number of trials per category for each VOT bin

[tbl,~, ~,labels] = crosstab(arrayfun(@(x) phnnames{x}, ...
    allplos(1, :), 'UniformOutput', false), allplos(2, :));

figure;
subplot(3, 1, [1 2]);
ord = {'p', 'b', 'k', 'g', 't', 'd'};
idx = cellfun(@(x) find(strcmp(labels(1:6, 1), x)), ord);
tbl = tbl(idx, 1:6);

imagesc(tbl, 'AlphaData', 0.5, 'XData',10:10:60);
for r = 1:6
    for c = 10:10:60
        text(c, r, num2str(tbl(r, c/10)), 'FontSize', 12);
    end
end

yticks(1:1:6);
yticklabels(ord);
xticks([]);
ylabel('phoneme category');
set(gca, 'FontSize', 15);
colormap(brewermap(10, 'YlGn'));
caxis([0 100]);

subplot(3, 1, 3);
uv = [sum(tbl(ismember(ord, {'p', 't', 'k'}), :)); ...
    sum(tbl(ismember(ord, {'b', 'g', 'd'}), :)); ...
    sum(tbl)];
imagesc(uv, 'AlphaData', 0.5, 'XData',10:10:60);
for r = 1:3
    for c = 10:10:60
        text(c, r, num2str(uv(r, c/10)), 'FontSize', 12);
    end
end
yticks(1:1:3);
yticklabels({'unvoiced', 'voiced', 'all'});
xlabel('VOT (ms)');
set(gca, 'FontSize', 15);
caxis([0 250]);

%% preceding phoneme in each VOT category
ord = {'p', 'b', 'k', 'g', 't', 'd'};
idx = ismember(TDcons.phn, ord);

figure;
[tbl, ~, ~, labels] = crosstab(TDcons.nextPhn(idx));
% subset to only look at cases where next phoneme happens in min+ trials
mintrl=20;
labels(tbl<mintrl) = [];

% vot = diff(TDcons.phnTimes);
ctr = 1;
allphn = zeros(length(ord), length(labels));
for phn = ord
    idx = ismember(TDcons.phn, phn);
    [phntbl, ~, ~, phnlab] = crosstab(TDcons.nextPhn(idx));
    phntbl(~ismember(phnlab, labels)) = [];
    phnlab(~ismember(phnlab, labels)) = [];

    % to match the label ordering of next phn
    phnlabidx = cellfun(@(x) find(strcmp(labels, x)), phnlab);
    allphn(ctr, phnlabidx) = phntbl(:);
    ctr = ctr+1;
end

imagesc(allphn,'AlphaData', 0.5);
for r = 1:length(ord)
    for c = 1:length(labels)
        text(c-0.2, r, num2str(allphn(r, c)), 'FontSize', 12);
    end
end
yticks(1:1:6);
yticklabels(ord);
xlabel('next phoneme');
xticks(1:length(labels));
xticklabels(labels)
set(gca, 'FontSize', 15);
colormap(brewermap(10, 'YlGn'));
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns;

%% VOT by succeeding phn (see previous cell for counts)

figure; 
ord = { 'k',  't', 'd'};
orth = {{'c', 'k', 'ck'}, {'t'}, {'d'}};

ctr = 1;
for plos = ord
    subplot(2, 3, ctr);
    idx = strcmp(TDcons.phn, plos);
    [tbl, ~, ~, labels] = crosstab(TDcons.nextPhn(idx));
    next = labels(tbl>15);

    vot = diff(TDcons.phnTimes)-1;
    TDcons.wrd()
    for n = 1:length(next)
        %scatter(ctr-0.1+rand(1, sum(idx))*0.2, vot(idx), 45, 'filled'); hold on;
        wrds = TDcons.wrd(idx & strcmp(TDcons.nextPhn, next{n}));
        c = 1;
        pos = nan(length(wrds), 1);
        for w = wrds
            tmp = cell2mat(cellfun(@(y) strfind(w{1},y), ...
                orth{ctr}, 'UniformOutput',false));
            if ~isempty(tmp), pos(c) = min(tmp); end
            c = c+1;
        end
        subplot(2, 3, ctr);
        boxplot(vot(idx & strcmp(TDcons.nextPhn, next{n})).*10, ...
            'Position', n, 'Colors','k'); hold on

        subplot(2, 3, ctr+3);
        boxplot(pos, 'Position', n, 'Colors','k'); hold on
    end
    subplot(2, 3, ctr);
    xticks(1:length(next));
    xticklabels(next);
    xlabel('Next phoneme')
    ylabel('VOT (ms)')
    set(gca, 'FontSize', 15);
    title(plos);
    box off;

    subplot(2, 3, ctr+3);
    xticks(1:length(next));
    xticklabels(next);
    xlabel('Next phoneme')
    ylabel('phoneme word position')
    set(gca, 'FontSize', 15);
    title(plos);
    box off;
    ctr = ctr +1;
end

next = {'l', 'ix', 'ah'};


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
betaInfo* *encoding* allidx fthresh *cons* phns;



%% TODO: spectrogram for specific stop-consonant category



%% comparing unvoiced distribution in timit and dimex

figure;
colors = getColorsCrossComp(4);
y = diff(Dcons.phnTimes(:, ismember(Dcons.phn, {'p', 'k', 't'})));
[fi, x] = ksdensity(y, 'Bandwidth', 0.55);

histogram(y*10, 9, 'Normalization', 'Probability', 'FaceColor', colors(1, :), ...
    'EdgeColor','none', 'FaceAlpha', 0.6);  hold on;
plot(x*10, fi, 'LineWidth', 3, 'Color',colors(1, :), 'HandleVisibility', 'off');

y = diff(TDcons.phnTimes(:, ismember(TDcons.phn, {'p', 'k', 't'})));
[fi, x] = ksdensity(y, 'Bandwidth', 0.55);

histogram(y*10, 9, 'Normalization', 'Probability', 'FaceColor', colors(1, :)+0.6, ...
    'EdgeColor','none', 'FaceAlpha', 0.8); 
plot(x*10, fi, 'LineWidth', 3, 'Color', colors(1, :)+0.4, 'HandleVisibility', 'off');
xlim([15 120]);

legend({'Spanish unvoiced', 'English unvoiced'});

xlabel('VOT (ms)');
ylabel('Probability');
yticks(0.2:0.2:0.8)
box off;
set(gca, 'FontSize', 15);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons*;

%% -------------------- Electrode selection: VOT --------------------------
%% calculate PSI to select plosive electrodes

% load in TRF models   
varnames = {'SID', 'el', 'ls', 'base_rsq', 'maxbeta', 'voicing_beta'};
vot_encoding = array2table(zeros(0, 6), 'VariableNames', varnames);
modelnames={'onset_phnfeatonset_maxDtL'};   % 'onset_phnfeatonset_maxDtL'{'onset_maxDtL_aud'};   % 
corpus = {'timit', 'dimex'};

% determine unique variance per feature and primary encoding
% [sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
for s = [sSIDs eSIDs bSIDs] % bSIDs%
        SID = s{1}; 
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs}));
%     corpusStrf=loadMultModelStrf(SID, modelnames(1), 'timit', ...
%         datapath, 1); % for the other two modelname lists
    corpusStrf{1} = loadMultModelStrf(SID, modelnames(1:end), corpus{1}, ...
        datapath, 1);
    corpusStrf{2} = loadMultModelStrf(SID, modelnames(1:end), corpus{2}, ...
        datapath, 1);

    % timit has more subject coverage (includes mandarin patients)
    if ~any(cellfun(@(x) isempty(x), corpusStrf{1}))
        numel = length(corpusStrf{1}{1}.meanTestR);
        basersq = ones(numel, 2)*-1;
        maxbeta = ones(numel, 2)*-1;

        % models without pitch
        basersq(:, 1)   = (corpusStrf{1}{1}.meanTestR.^2)';
        % extract only betas from the phonetic feature sets
        [~,maxbeta(:, 1)]    = max(squeeze(mean(corpusStrf{1}{1}.meanStrf(2:15, :, :), 2)));
        voicing_beta = squeeze(corpusStrf{1}{1}.meanStrf(1+12, :, :));

        if ~any(cellfun(@(x) isempty(x), corpusStrf{2}))
            numel_dimex = length(corpusStrf{2}{1}.meanTestR);
            basersq(1:numel_dimex, 2)   = (corpusStrf{2}{1}.meanTestR.^2)';
            [~,maxbeta(1:numel_dimex, 2)]    = ...
                max(squeeze(mean(corpusStrf{2}{1}.meanStrf(2:15, :, :), 2)));

        end
    
        % find which electrodes meet threshold rsq for the full bin model
        els = find(sum(basersq > 0.01, 2)>0);
        
        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        
        %uvPeakrate(els)
        voicing = mat2cell(voicing_beta(:, els), 61, ones(length(els), 1));
        tmp = table(sids, els, lss, basersq(els, :), maxbeta(els, :), ...
            voicing', 'VariableNames', ...
            varnames);
        vot_encoding = [vot_encoding; tmp];
    else
        warning(['Missing subject ' SID])
    end
end

% selecting for plosive, voicing, or obstruent
% idx = sum(ismember(vot_encoding.maxbeta, [8, 12, 13]), 2)>0;
% voicing_elec = voicing_elec(idx, :);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns;

%% calculate spearmans rho for VOT

rho = nan(1, height(vot_encoding)); 
p = nan(1, height(vot_encoding)); 
% determine unique variance per feature and primary encoding
phns = {'b', 'p', 'g', 'k', 'd', 't'};
vot_max = 8;
mean_resp = nan(height(vot_encoding), vot_max-2);

% per category pair
phn_pairs = {{'b', 'p'}, {'g', 'k'}, {'d', 't'}};
rho_pairs = nan(length(phn_pairs), height(vot_encoding));
p_pairs = nan(length(phn_pairs), height(vot_encoding)); 
mean_pairesp = nan(length(phn_pairs), height(vot_encoding), vot_max-2);
maxtp = nan(1, height(vot_encoding));
for e = 1:height(vot_encoding)
    SID = vot_encoding.SID{e};    
    el = vot_encoding.el(e);

    if isfield(TDcons, SID)
        if el <= size(TDcons.(SID).resp, 1)
            vot = diff(TDcons.phnTimes);
            % vot = TDcons.dur*100;
            phnidx = ismember(TDcons.phn, phns) & vot<vot_max;
            % can subset to syllable onset, word initial, or next phn stress
    
            resp = squeeze(TDcons.(SID).resp(el, :, :));
            % max window
            [~, maxtp(e)] = max(mean(resp, 2, 'omitnan'));
            peak_resp = mean(resp(max(maxtp(e)-3, 1):min(maxtp(e)+3, ...
                bef+aft), :), 'omitnan');
            
            nanidx = isnan(peak_resp);
            if sum(nanidx)<0.8*length(phnidx)
                [rho(e), p(e)] = corr(peak_resp(phnidx&~nanidx)', ...
                    vot(phnidx&~nanidx)', 'Type', 'Spearman');
                
                % mean groups, non-downsampled
                % [means,grps] = grpstats(peak_resp(phnidx&~nanidx), ...
                % discretize(vot(phnidx&~nanidx), [0:1:5 20]), ["mean", "gname"]);
                % mean_resp(e, str2double(grps)) = means;
    
                % downsampled
                [means,grps] = grpstats(peak_resp(phnidx&~nanidx), ...
                    vot(phnidx&~nanidx), ["mean", "gname"]);
                mean_resp(e, str2double(grps)-1) = means;
    
                % per pair p and mean resp
                ctr = 1;
                for pair = phn_pairs
                    phnpair_idx = phnidx & ismember(TDcons.phn, pair{1});
                    [rho_pairs(ctr, e), p_pairs(ctr, e)] = corr(peak_resp(phnpair_idx...
                        &~nanidx)', vot(phnpair_idx&~nanidx)', 'Type', 'Spearman');
    
                    % nondownsampled
                    % [mean_pair,grps] = grpstats(peak_resp(phnidx&~nanidx), ...
                    % discretize(vot(phnidx&~nanidx), [0:1:5 20]), ["mean", "gname"]);
                    % mean_pairesp(ctr, e, str2double(grps)) = mean_pair;
    
                    % downsampled
                    [mean_pair,grps] = grpstats(peak_resp(phnidx&~nanidx), ...
                    vot(phnidx&~nanidx), ["mean", "gname"]);
                    mean_pairesp(ctr, e, str2double(grps)-1) = mean_pair;
                    clear mean_pair grps
                    ctr = ctr + 1;
                end
            end
        end
    end
end

vot_encoding.rho = rho';
vot_encoding.p = p';
vot_encoding.resp = mean_resp;
dataf = 100;
vot_encoding.maxtp =(maxtp'-bef)./dataf;

vot_encoding.p_pairs = p_pairs';
vot_encoding.rho_pairs = rho_pairs';

mean_pairesp = mat2cell(mean_pairesp, length(phn_pairs), ones(height(vot_encoding), 1), ...
    vot_max-2);
mean_pairesp = cellfun(@(x) squeeze(x), mean_pairesp, 'UniformOutput', false);
vot_encoding.resp_pairs = mean_pairesp';

% selecting for plosive, voicing, or obstruent (primary encoding)
% idx = sum(ismember(vot_encoding.maxbeta, [8, 12, 13]), 2)>0;

% the significance testing changes depending on which subjects are included
% and which threshold is used! Please note when reporting results.
vot_encoding = vot_encoding(vot_encoding.p<0.05, :); 
%vot_encoding = vot_encoding(abs(vot_encoding.rho)>0.05, :);

% selected electrodes change depending on whether you use 8, 9 or 10 VOT
% bins, selects more electrodes when fewer more dense bins are used

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% visualizing average tuning curves 

plus_idx = vot_encoding.rho>0.0;

ctr = 1;
resp_type = cell(2, 1);
stats_y = cell(2, 2);
for r = {plus_idx; ~plus_idx}'
    resp_type{ctr} = [vot_encoding.resp(r{1}, :)];
    % divide by max
    % v_plus = bsxfun(@rdivide,v_plus,max(v_plus, [], 2));

    % rescale from 0 to 1
    resp_type{ctr} = rescale(resp_type{ctr},'InputMin', ...
        min(resp_type{ctr},[],2), 'InputMax', max(resp_type{ctr},[],2));
    ctr = ctr + 1;
end

figure;
vot = 10:10:(vot_max-2)*10;

% old red / blue 
% pt_cols = getColors(3, length(vot));
% line_cols = [0 0 1; 1 0 0];
pt_cols = getColorsCrossComp(5);
line_cols = getColorsCrossComp(4);
style = {':', '-', '',  '-.'};
for ls = [1, 2, 4]
    el_type = {vot_encoding.ls(plus_idx)==ls, vot_encoding.ls(~plus_idx)==ls};
    disp(['ls: ' num2str(ls) ' minus:' num2str(sum(el_type{2})) ...
        ', plus:' num2str(sum(el_type{1}))]);

    for t = 1:2
        subplot(1, 2, t)
        y = smoothdata(resp_type{t}(el_type{t}, :), 2, ...
            'SmoothingFactor',0.01);    
        stats_y(ls, t) = {y}; 
        if ~isempty(y)
            errorbar(vot, mean(y), nansem(y),'Color', line_cols(t, :), ...
                'LineWidth', 2, 'CapSize', 0,  ...
                'LineStyle', style{ls}, 'HandleVisibility','on'); hold on;
            scatter(vot, mean(y), 75, pt_cols, 'filled', 'HandleVisibility','off');
        end
        hold on
        ylim([0.0 1.01]);
        xlim([5 vot(end)+5]);
        xlabel('VOT binned (ms)');
        ylabel('peak HFA (% of max)')
        yticks([0 0.5 1]);
        yticklabels({'0', '50', '100'});        

        set(gca, 'FontSize', 15);
        box off;
    end
    ylim([0.1 1.01]);
end

% check significance
for t = 1:2
    span = stats_y{1, t};
    eng = stats_y{2, t};

    for v = 1:vot_max-2
        [p, ~] = ranksum(span(:, v), eng(:, v));
        disp(p)
        if p < 0.01
            subplot(1, 2, t)
            text(v*10, 1.05, '*', 'FontSize', 20);
        end
    end
end

legend({'Spanish', 'English', 'Bilingual'});


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;
%% visualize average fitted curves

plus_idx = vot_encoding.rho>0.0;

ctr = 1;
resp_type = cell(2, 1);
stats_y = cell(2, 2);
for r = {plus_idx; ~plus_idx}'
    resp_type{ctr} = [vot_encoding.resp(r{1}, :)];
    % divide by max
    % v_plus = bsxfun(@rdivide,v_plus,max(v_plus, [], 2));

    % rescale from 0 to 1
    resp_type{ctr} = rescale(resp_type{ctr},'InputMin', ...
        min(resp_type{ctr},[],2), 'InputMax', max(resp_type{ctr},[],2));
    ctr = ctr + 1;
end

figure;
vot = 10:10:(vot_max-2)*10;
pt_cols = getColorsCrossComp(5);
line_cols =getColorsCrossComp(4); %[0 0 1; 1 0 0];
style = {':', '-', '-.'};
for ls = 1:2
    el_type = {vot_encoding.ls(plus_idx)==ls, vot_encoding.ls(~plus_idx)==ls};
    disp(['ls: ' num2str(ls) ' minus:' num2str(sum(el_type{2})) ...
        ', plus:' num2str(sum(el_type{1}))]);

    for t = 1:2
        subplot(2, 2, t)
        y = resp_type{t}(el_type{t}, :);    
        sids = vot_encoding.SID(el_type{t});
        stats_y(ls, t) = {y}; 

        % fitted curves
        if ~isempty(y)
            fsigm = @(param,xval) param(1)+(param(2)-param(1))./...
                    (1+10.^((param(3)-xval)*param(4)));
            y_fit=nan(size(y, 1), size(y, 2));
            for i = 1:size(y, 1)
                [param, stat] = sigm_fit(vot, y(i, :), [], [0 1 0.65 0], 0);
                if stat.rsq>0.15
                    y_fit(i, :) = fsigm(param, vot);
                else
                    y_fit(i, :) = nan(1, length(vot));
                end
            end
        end

        if ~isempty(y_fit)
            errorbar(vot, mean(y_fit, 'omitnan'), nansem(y_fit),'Color', line_cols(t, :), ...
                'LineWidth', 2, 'CapSize', 0,  ...
                'LineStyle', style{ls}, 'HandleVisibility','on'); hold on;
            scatter(vot, mean(y_fit, 'omitnan'), 75, pt_cols, 'filled', 'HandleVisibility','off');
        end
        hold on
        ylim([0.0 1.01]);
        xlim([5 vot(end)+5]);
        xlabel('VOT binned (ms)');
        ylabel('peak HG (% of max)')
        yticks([0 0.5 1]);
        yticklabels({'0', '50', '100'});        

        set(gca, 'FontSize', 15);
        ylim([0.0 1.01])
        box off;
        legend({'Spanish', 'English', 'Bilingual'});

        % showing average tuning curve per subject
        subplot(2, 2, t+2)
        for s = unique(sids)'
            SID = s{1};
%             sz = sum(strcmp(sids, SID))*0.25;
            sz = sum(abs(vot_encoding.rho(strcmp(sids, SID))));
            plot(mean(y_fit(strcmp(sids, SID), :))', ...
                'LineWidth', sz, 'Color', line_cols(t, :), ...
                'LineStyle', style{ls}, 'HandleVisibility','on'); hold on;
        end
        xlim([0.5 0.1*vot(end)+0.5]);

    end
    ylim([0.0 1.01]);
end

% check significance
for t = 1:2
    span = stats_y{1, t};
    eng = stats_y{2, t};

    for v = 1:vot_max-2
        [p, ~] = ranksum(span(:, v), eng(:, v));
        
        subplot(2, 2, t)
        text(v*10, 1.05, getSigStr(p, 2), 'FontSize', 10);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;


%% visualize single electrode tuning

plus_idx = vot_encoding.rho>0;

ctr = 1;
resp_type = cell(2, 1);
stats_y = cell(2, 2);
for r = {plus_idx; ~plus_idx}'
    resp_type{ctr} = [vot_encoding.resp(r{1}, :)];
    % divide by max
    % v_plus = bsxfun(@rdivide,v_plus,max(v_plus, [], 2));

    % rescale from 0 to 1
    resp_type{ctr} = rescale(resp_type{ctr},'InputMin', ...
        min(resp_type{ctr},[],2), 'InputMax', max(resp_type{ctr},[],2));
    ctr = ctr + 1;
end

vot = 10:10:(vot_max-2)*10;
% line_cols = [0 0 1; 1 0 0];
line_cols = getColorsCrossComp(4);

figure;
subplot(1, 2, 1);
% [~, i] = sortrows(discretize(resp_type{1}, ...
%     0:0.1:1), [6, 1, 2, 5, 4, 3]);
imagesc(resp_type{1}); % (i, :)
xticklabels(split(num2str(10:10:60)));
xlabel('VOT bin (ms)')
ylabel('Electrode Count')
set(gca, 'FontSize', 13);

subplot(1, 2, 2);
% [~, i] = sortrows(discretize(resp_type{2}, ...
%     0:0.1:1), [1, 6, 2, 5, 4, 3]);
imagesc(resp_type{2}); % (i, :)
colormap(brewermap(15, 'Greys'))
xticklabels(split(num2str(10:10:60)));
xlabel('VOT bin (ms)')
ylabel('Electrode Count')
set(gca, 'FontSize', 13);
cbh = colorbar;
ylabel(cbh, 'peak of HFA (% max)')

figure;
infl = cell(2, 2);
for ls = 1:2
    el_type = {vot_encoding.ls(plus_idx)==ls, vot_encoding.ls(~plus_idx)==ls};
    idxs = {vot_encoding.ls==ls&plus_idx, vot_encoding.ls==ls&~plus_idx};
    for t = 1:2
        subplot(2, 2, (ls-1)*2+t)
        y = resp_type{t}(el_type{t}, :); 
        wdt = vot_encoding.rho(idxs{t});
        stats_y(ls, t) = {y}; 
        if ~isempty(y)

            fsigm = @(param,xval) param(1)+(param(2)-param(1))./...
                    (1+10.^((param(3)-xval)*param(4)));
            y_fit=nan(size(y, 1), size(y, 2));
            for i = 1:size(y, 1)
                [param, ~] = sigm_fit(vot, y(i, :), [], [0 1 0.65 0]);
                y_fit(i, :) = fsigm(param, vot);
                if param(3)>0 && param(3)<100
                    infl{t, ls}(i) = param(3);
                else
                    infl{t, ls}(i) = NaN;
                end
            end
            h = plot(vot, y_fit, 'Color', line_cols(t, :));

            for i = 1:length(h)
                h(i).Color = [h(i).Color min(abs(wdt(i))*3, 1)]; 
                h(i).LineWidth = abs(wdt(i))*9; 
            end
        end
        hold on
        xlim([5 vot(end)+5]);
        xlabel('VOT binned (ms)');
        ylabel('peak HFA (% of max)')
        yticks([0 0.5 1]);
        yticklabels({'0', '50', '100'});        

        set(gca, 'FontSize', 15);
        box off;
    end
    ylim([0 1.01]);
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% visualize latency effects

% figure;
scatter(vot_encoding.maxtp, vot_encoding.rho, ...
    discretize(1-vot_encoding.p, 10)*7, vot_encoding.rho<0, 'filled', ...
    'MarkerFaceAlpha', 0.6);
colormap(getColorsCrossComp(4));
yline(0);
ylabel('rho (direction of tuning)');
xlabel('peak time point (from consonant onset)');
ylim([-0.12 0.16])
xlim([0 0.55]);
set(gca, 'FontSize', 15);

figure;
colors = getColorsCrossComp(4);
plus_idx = vot_encoding.rho>0;
scatter(vot_encoding.maxtp(plus_idx), vot_encoding.p(plus_idx), 45, ...
    colors(1, :), 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName','V-'); hold on;
scatter(vot_encoding.maxtp(~plus_idx), vot_encoding.p(~plus_idx), 45, ...
    colors(2, :), 'filled', 'MarkerFaceAlpha', 0.6, 'DisplayName','V+');
ylabel('p-val');
xlabel('HG peak relative to onset (s)');

ylim([-0.01 0.062])
xlim([0 0.55]);
yticks(0:0.02:0.06);
xticks(0:0.2:0.6);
set(gca, 'FontSize', 15);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% plot tuning curve by category type

cols = getColors(1);
vot = 10:10:(vot_max-2)*10;
pt_cols = getColors(3, length(vot));
line_cols = [0 0 1; 1 0 0];
style = {':', '-', '-.'};
labels = {'Spanish', 'English', 'Bilingual'};
phn_pairs = {{'b', 'p'}, {'g', 'k'}, {'d', 't'}};

figure;
for pair = 1:3

    sig_idx = vot_encoding.p_pairs(:, pair)<0.06;
    plus_idx = vot_encoding.rho>0 & sig_idx;

    ctr = 1;
    resp_type = cell(2, 1);
    for r = {plus_idx; ~plus_idx&sig_idx}'
        resp_type{ctr} = cell2mat(cellfun(@(x) x(pair, :), ...
            [vot_encoding.resp_pairs(r{1})], 'UniformOutput', false));
        % divide by max
        % v_plus = bsxfun(@rdivide,v_plus,max(v_plus, [], 2));
    
        % rescale from 0 to 1
        resp_type{ctr} = rescale(resp_type{ctr},'InputMin', ...
            min(resp_type{ctr},[],2), 'InputMax', max(resp_type{ctr},[],2));
        ctr = ctr + 1;
    end

    disp(['------------ ' strjoin(phn_pairs{pair}) ' -----------------'])
    for ls = 1:2
        el_type = {vot_encoding.ls(plus_idx)==ls, ...
            vot_encoding.ls(~plus_idx & sig_idx)==ls};
        
        disp([ labels{ls} ': minus:' num2str(sum(el_type{2})) ...
            ', plus:' num2str(sum(el_type{1}))]);
    
        for t = 1:2
            subplot(3, 2, t + (pair-1)*2)
            if sum(el_type{t})>0
                y = resp_type{t}(el_type{t}, :);  
                y_mean = y;
                if size(y, 1)>1
                    y_mean = mean(y, 1);
                end

                errorbar(vot, y_mean, nansem(y, 1),'Color', line_cols(t, :), ...
                    'LineWidth', 2, 'CapSize', 0,  ...
                    'LineStyle', style{ls}, 'HandleVisibility','on'); hold on;
                scatter(vot, y_mean, 75, pt_cols, 'filled', 'HandleVisibility','off');
                hold on
                ylim([0.0 1.05]);
                xlim([5 vot(end)+5]);
                yticks(0:0.5:1);
                xticks([10 vot(end)]);
                set(gca, 'FontSize', 15);
                legend('off');
                box off;
            end
        end

        if pair == 3
            xticks(vot)
            xlabel('VOT binned (ms)');
            legend('on');
        end
     
        ylim([0 1.05]);
    end
    legend(labels);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% p_vals by category type
% are the same electrodes correlated to VOT across category?
p_thresh = 0.06;
phn_pairs = {{'b', 'p'}, {'g', 'k'}, {'d', 't'}};
labels = cellfun(@(x) strjoin(x, ', '), phn_pairs, 'UniformOutput',false);

figure;
comps = [1, 2; 1, 3; 2, 3];
for c = 1:3
    subplot(2, 3, c)
    x = vot_encoding.p_pairs(:, comps(c, 1));
    y = vot_encoding.p_pairs(:, comps(c, 2));

    scatter(x, y, 45, vot_encoding.rho>0, 'filled', ...
        'HandleVisibility', 'off'); 
    caxis([0 1]);
    colormap([1 0 0; 0 0 1]);
    yline(p_thresh, 'LineWidth', 1.5); 
    xline(p_thresh, 'LineWidth', 1.5);
    
    xlabel(['(' strjoin( phn_pairs{comps(c, 1)}, ' ,') ')']);
    ylabel(['(' strjoin( phn_pairs{comps(c, 2)}, ' ,') ')'],...
        'HandleVisibility','off');
    legend('on');
    legend({'p thresh'});

    set(gca, 'Yscale', 'log',  'Xscale', 'log', 'FontSize', 15);

    % show rho/effect size
    subplot(2, 3, c+3)
    x_rho = vot_encoding.rho_pairs(:, comps(c, 1));
    y_rho = vot_encoding.rho_pairs(:, comps(c, 2));

    szs = [10, 40, 80];
    sz = arrayfun(@(x) szs(x), sum([x<p_thresh y<p_thresh], 2)+1);
    scatter(x_rho, y_rho, sz, ...
        vot_encoding.rho>0, 'filled', 'MarkerFaceAlpha', 0.6,...
        'HandleVisibility', 'off'); 
    caxis([0 1]);
    colormap([1 0 0; 0 0 1]);  
    xline(0, 'LineWidth', 1.5); 
    yline(0, 'LineWidth', 1.5); 
    refline(1, 0);
    
    xlabel(['(' strjoin( phn_pairs{comps(c, 1)}, ' ,') ')']);
    ylabel(['(' strjoin( phn_pairs{comps(c, 2)}, ' ,') ')'],...
        'HandleVisibility','off');
    legend('on');

    set(gca, 'FontSize', 15);
end

figure;
subplot(1, 3, [1 2])
[x, col] = min(vot_encoding.p_pairs, [], 2);
col(x>p_thresh)=0;
y = vot_encoding.p;
bm = [0 0 0; 0 0.8 0.6; 0.35 0 0.85; 0.85 0.42 0.81; 0.53 0.32 0.65];
%col = cell2mat(arrayfun(@(x) bm(x, :), col+1, 'UniformOutput',false));
gscatter(x, y, col, bm, ".", 25);
xlabel('min p-val for category pair');
ylabel('p-val across category pairs');
set(gca, 'FontSize', 15, 'Xscale', 'log');
xline(p_thresh, 'LineWidth', 1.5); 
legend([{'none'} labels {'p thresh'}])


subplot(1, 3, 3)
p_pairs = [vot_encoding.p_pairs];
p_pairs(p_pairs>p_thresh) = 1;
imagesc(p_pairs);
colormap(brewermap(2, 'RdGy'));
ylabel('Electrode');
xlabel('Category pair');
cbh = colorbar;
xticks(1:3);
xticklabels(labels)
cbh.Ticks = [0.25 0.75];
cbh.TickLabels = {'sig', 'n.s.'};
set(gca, 'FontSize', 15);


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% FIX: plotting single VOT erps

%meeting slides
% V+: Spanish - EC172 el219, English - EC222 el125 
% V-: Spanish - EC152 el101, English - EC186 el231

% Quals slides
% V-: Spanish - EC105 el 185, English - EC195 el137
% V+: Spanish - EC100 el70, English - EC195 el217

phns = {'g', 'k', 'b', 'p' , 'd', 't'}; % , ,  
rmpath('util/shadederror');
% % SID = 'EC266'; % 'EC222'
% % el = 134;, 122
SID = 'EC100'; % 'EC222'
el = 134;


plotVOTerp(TDcons, SID, el, 'timit', phns, 0);
legend('FontSize', 15);
yticks(0:0.5:1)
xticks(-0.2:0.2:0.4);
set(gca, 'FontSize', 15);

figure; 
votidx = find(vot_encoding.el == el & strcmp(vot_encoding.SID, SID));
vot = 10:10:(8-2)*10;
pt_cols = flipud(getColorsCrossComp(5));
y_mean = rescale(vot_encoding.resp(votidx, :),'InputMin', ...
        min(vot_encoding.resp(votidx, :),[],2), ...
        'InputMax', max(vot_encoding.resp(votidx, :),[],2));

scatter(vot, y_mean, 75, pt_cols, 'filled', 'HandleVisibility','off'); hold on;

fsigm = @(param,xval) param(1)+(param(2)-param(1))./...
                    (1+10.^((param(3)-xval)*param(4)));
[param, stat] = sigm_fit(vot, y_mean, [], [0 1 0.65 0], 0);
y_fit = fsigm(param, vot);
plot(vot, y_fit, 'Color', 'k', 'LineWidth', 1.7);

% formatting
ylabel({'peak HFA', '(% of max)'});
xlabel('VOT binned (ms)');
yticks([0 1]);
yticklabels({'0', '100'});
xticks(0:20:60)
xlim([8, 59])
set(gca, 'FontSize', 13);

% f = figure; 
% subplot(1, 3, 1);
% addpath('util/shadederror');
% plotPhonemeErp(TDcons, SID, el, 'timit', f, [1 0 0; 0 0 1], {'p'; 'b'});
% legend('FontSize', 15);
% set(gca, 'FontSize', 15);
% subplot(1, 3, 2);
% plotPhonemeErp(TDcons, SID, el, 'timit', f, [1 0 0; 0 0 1], {'k'; 'g'});
% legend('FontSize', 15);
% yticks([]);
% set(gca, 'FontSize', 15);
% subplot(1, 3, 3);
% plotPhonemeErp(TDcons, SID, el, 'timit', f, [1 0 0; 0 0 1], {'t'; 'd'});
% legend('FontSize', 15);
% yticks([]);
% set(gca, 'FontSize', 15);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% showing single TRF weights

row = 43; % 11, 43
modelnames={'onset_phnfeatonset_maxDtL', 'onset_aud'};   
% 'onset_phnfeatonset_maxDtL'{'onset_maxDtL_aud'};   % 
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};


el = vot_encoding.el(row);
SID = vot_encoding.SID{row};

figure
for l = 1:2
    corpusStrf = loadMultModelStrf(SID, modelnames, corpus{l}, ...
        datapath, 1);
    featnames = details{l}.features.names;
    
    subplot(2, 2, l);
    imagesc(squeeze(corpusStrf{1}.meanStrf(2:end-1, :, el)));
    if l == 1
        yticks(1:2+length(featnames));
        yticklabels((featnames')) % {'onset'}, , {'peakRate'}
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
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% showing frequency selectivity across VOT encoding types

modelnames={'onset_phnfeatonset_maxDtL', 'onset_aud'};   % 'onset_phnfeatonset_maxDtL'{'onset_maxDtL_aud'};   % 
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};

aud_weights = nan(height(vot_encoding), 80); % max tp frequency selectivity
feat_weights = nan(height(vot_encoding), 14);

vot_rsq = [];
rsq_idx = [];
for s = unique(vot_encoding.SID)'
    SID = s{1};
    
    corpusStrf = loadMultModelStrf(SID, modelnames, corpus{1}, ...
        datapath, 1);

    if ~isempty(corpusStrf{2})

        els = strcmp(vot_encoding.SID, SID);
        vot_rsq = [vot_rsq corpusStrf{2}.meanTestR(vot_encoding.el(els)).^2];
        rsq_idx = [rsq_idx find(els)'];

        for row = find(els)'
            el = vot_encoding.el(row);
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

            if corpusStrf{2}.meanTestR(el)^2>0.001
                aud_weights(row, :)  = y;
            end
        end
    end
end

v_minus = find(vot_encoding.rho>0 & ~isnan(mean(aud_weights, 2)));
v_plus = find(vot_encoding.rho<0 & ~isnan(mean(aud_weights, 2)));
rhos = abs(vot_encoding.rho);

[g, ~, ~, labels] = crosstab(max_feat, vot_encoding.rho>0);

x = cellfun(@(x) str2num(x), labels(:, 1));
b = bar(x, g, 'stacked');
b(1).FaceColor = [1, 0, 0];
b(2).FaceColor = [0, 0, 1];

xticks(1:14);
xlim([0.4 14.8]);
xticklabels(featnames);

ylabel('Electrode count');
xlabel('Primary Phonetic Feature');
set(gca, 'FontSize', 15);
box off;


% r-squared comparison
figure;
scatter(vot_rsq, vot_encoding.base_rsq(rsq_idx, 1), 55, ...
    vot_encoding.rho(rsq_idx)<0, 'filled', ...
    'MarkerFaceAlpha', 0.7); 
colormap([0 0 1; 1 0 0]);

ylim([0 0.42]);
xlim([0 0.42]);

h=refline(1, 0);
h.LineWidth = 1.5;
h.LineStyle = '-';
h.Color = 'k';
ylabel('phonetic feature R^2');
xlabel('mel spectrogram R^2');
set(gca, 'FontSize', 15);

% beta weight comparison (what does the average frequency selectivity look like?)
figure;   
subplot(2, 1, 1);
[~, tpidx] = max(aud_weights(v_minus, :), [], 2);
[~, sortidx] = sort(tpidx);
imagesc(aud_weights(v_minus(sortidx), :)); %.*rhos(v_minus(sortidx)));
ylabel('Electrode count');
xlabel('Frequency Bin');
set(gca, 'FontSize', 15);
title('V- population');
box off;

subplot(2, 1, 2);
[~, tpidx] = max(aud_weights(v_plus, :), [], 2);
[~, sortidx] = sort(tpidx);
imagesc(aud_weights(v_plus(sortidx), :)); %.*rhos(v_plus(sortidx)));
ylabel('Electrode count');
xlabel('Frequency Bin');
set(gca, 'FontSize', 15);
title('V+ population');
box off;

cm = brewermap(50, 'Greys');
colormap(cm(1:end-8, :));

cbh = colorbar();
ylabel(cbh, 'STRF weight (% max)');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% number of electrodes that prefer V+ or V- across subject groups

figure;
subplot(1, 2, 1)
b=bar([1, 2], [crosstab(vot_encoding.rho(vot_encoding.ls==2)<0)./length(vot_encoding.rho(vot_encoding.ls==2)), ...
    crosstab(vot_encoding.rho(vot_encoding.ls==1)<0)./length(vot_encoding.rho(vot_encoding.ls==1))]');
% , crosstab(vot_encoding.rho(vot_encoding.ls==4)<0)./length(vot_encoding.rho(vot_encoding.ls==4))
b(1).FaceColor = [0.1 0 0.8];
b(2).FaceColor = [0.8 0.1 0.1];

xticklabels({'English', 'Spanish', 'Bilingual'});
yticks([0 0.5 1])
ylim([0 1])
ylabel('% of electrodes');
xlabel('Subject group');
set(gca, 'FontSize', 15);
legend('V-', 'V+');
box off;

vplus_prct = cell(1, 4);
vplus_cnt = cell(1, 4);
for s = unique(vot_encoding.SID)'
    SID = s{1};
    sidx = find(strcmp(vot_encoding.SID, SID));
    
    ls = vot_encoding.ls(sidx(1));
    vplus_prct{ls} = ...
        [vplus_prct{ls} sum(vot_encoding.rho(sidx)<0)/length(sidx)];
    vplus_cnt{ls} = [vplus_cnt{ls}, length(sidx)];
end
xlim([0.5 2.5]);

subplot(1, 2, 2)
for ls = [1, 2]
    jitter = rand(length(vplus_prct{ls}), 1)*0.2; 
    sz = vplus_cnt{ls}*10;
    scatter(ones(length(vplus_prct{ls}), 1)*ls-0.1+jitter, ...
        vplus_prct{ls}, sz, [0.5 0.5 0.5], 'filled', 'MarkerFaceAlpha', 0.4); hold on;
    b = boxchart(ones(length(vplus_prct{ls}), 1)*ls, vplus_prct{ls});
    b.BoxFaceColor = [0.7 0.7 0.7];
    b.WhiskerLineStyle = 'none';
    b.MarkerStyle = 'none';
    disp(sz)
end
xticks([1, 2, 4]);
xticklabels({'Spanish', 'English', 'Bilingual'});
xlabel('Subject group');
ylabel('% V+ electrodes')
yticks([0 0.5 1]);
box off;
set(gca, 'FontSize', 15);
xlim([0.25 2.75]);
ylim([0 1.1])

tbl = table();
% rho direction not value
lss = [1, 2];
% tbl.rho = double(vot_encoding.rho(ismember(vot_encoding.ls, lss))>0);
tbl.rho = vot_encoding.rho(ismember(vot_encoding.ls, lss));
tbl.sid = vot_encoding.SID(ismember(vot_encoding.ls, lss));
tbl.ls = vot_encoding.ls(ismember(vot_encoding.ls, lss));
tbl.el = vot_encoding.el(ismember(vot_encoding.ls, lss));
% tbl.ls(tbl.ls==4) = 2;
lme = fitlme(tbl,'rho~ls+(1|sid)+(1|el:sid)');
disp(lme)
title(lme.Coefficients{2, 6});
% fitlme(tbl,'rho~ls')

%% MNI elecs for direction of VOT encoding

count = nan(1, 2);
colors = getColorsCrossComp(4);
for ls = 1:2
    % initialize design electrode structure
    desel=struct();
    desel.conds = 0:6;
    desel.sz = [60 60];
    desel.cols = getColorsCrossComp(4);

    % split up peak rate and phonetic features again for MNI plotting
    desel.labels = {'V-', 'V+'};
    for s=unique(vot_encoding.SID)'
        SID = s{1};
        
        idx = find(strcmp(vot_encoding.SID, SID));
    
        if ismember(vot_encoding.ls(idx(1)), ls)
            desel.(SID).elid = vot_encoding.el(idx);
            desel.(SID).condition = (vot_encoding.rho(idx)<0);
        end
    end
    
    [mni_lh] = plotMNIElec(unique(vot_encoding.SID), desel, 'lh', 0);
    title(ls)
    % add a pie
%     axes('Position',[.7 .1 .3 .3])
%     p = pie([sum(mni_lh.cond>0), sum(mni_lh.cond<1)], [1 1]); 
%     p(1).FaceColor = [1, 0, 0];
%     p(1).EdgeColor = 'none';
%     p(3).FaceColor = [0, 0, 1];
%     p(3).EdgeColor = 'none';
    
    [mni_rh] = plotMNIElec(unique(vot_encoding.SID), desel, 'rh', 0);
    title(ls)
    % add a pie
%     axes('Position',[.7 .1 .3 .3])
%     p = pie([sum(mni_rh.cond>0), sum(mni_rh.cond<1)], [1 1]); 
%     p(1).FaceColor = [1, 0, 0];
%     p(1).EdgeColor = 'none';
%     p(3).FaceColor = [0, 0, 1];
%     p(3).EdgeColor = 'none';
    count(ls, :) = [sum(mni_lh.cond>0)/length(mni_lh.cond), ...
        sum(mni_rh.cond>0)/length(mni_rh.cond)];

end

figure; 
b = bar(count');
cols = [colors(2, :)-0.15; colors(2, :)+0.15];
for i = 1:2
    b(i).FaceColor = cols(i, :);
    b(i).EdgeColor = 'none';
end
colormap();
set(gca, 'FontSize', 15, 'Xtick', 1:2, 'Xticklabel', {'LH', 'RH'});
leg = legend({'Spanish native', 'English native'}, 'Location', 'northwest   ');
htitle = get(leg,'Title');
set(htitle,'String','Subject group')

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% Native plot

% initialize design electrode structure
desel=struct();
desel.conds = 1:5;
desel.sz =  [fliplr(10:30:90) 40:30:90];
b = balanced(8);
desel.cols = [1 0 0; b(4, :); 0 0 0; b(5, :); 0 0 1];
ls = [1, 2];
% split up peak rate and phonetic features again for MNI plotting
desel.labels = 1:5; %{'V-', 'V+'};
conds = discretize(vot_encoding.rho, [-0.2 linspace(-0.1,0.1, 4) 0.2]);
for s=unique(vot_encoding.SID)'
    SID = s{1};
    
    idx = find(strcmp(vot_encoding.SID, SID));

    if ismember(vot_encoding.ls(idx(1)), ls)
        desel.(SID).elid = vot_encoding.el(idx);
        desel.(SID).condition = conds(idx); %(vot_encoding.rho(idx)<0);
    end
end

plotNativeElec(unique(vot_encoding.SID), desel, 1)

%% TO FIX: looking at curves for non-downsampled VOT

ord = {'p', 'b', 'k', 'g', 't', 'd'};
idx = ismember(TDcons.phn, ord);

figure;
for r = 1:height(vot_encoding)
    SID = vot_encoding.SID{r};
    el = vot_encoding.el(r);
    rho = vot_encoding.rho(r);

    resp = squeeze(TDcons.(SID).resp(el, :, :));

    % max window
    [~, maxtp] = max(mean(resp, 2, 'omitnan'));
    peak_resp = mean(resp(max(maxtp-3, 1):min(maxtp+3, ...
        bef+aft), :), 'omitnan');
        
    peak_resp = rescale(peak_resp,'InputMin', ...
        min(peak_resp,[],2), 'InputMax', max(peak_resp,[],2)); 
    binedges = 0:0.005:0.08;
    dur_bin = discretize(TDcons.dur, binedges);

    % scatter(TDcons.dur(idx), peak_resp_bin(idx), 25, 'filled'); hold on;
    binwidth = (binedges(2) - binedges(1))/2;

    if rho<0
        subplot(1, 2, 1);      
    else
        subplot(1, 2, 2);
    end

    bin_peak_resp  =  grpstats(peak_resp(idx), dur_bin(idx));
    bin_peak_resp = rescale(bin_peak_resp,'InputMin', ...
        min(bin_peak_resp), 'InputMax', max(bin_peak_resp)); 
    scatter(binedges(1:end-1)+binwidth, bin_peak_resp, 25, [0.5 0.5 0.5], ...
        'filled', 'MarkerFaceAlpha', abs(rho)/max(abs(vot_encoding.rho))); hold on;
    line(binedges(1:end-1)+binwidth, bin_peak_resp, '0.5'); hold on;
    

%     [y, x] = grpstats(peak_resp, discretize(TDcons.dur, 0.02:0.01:0.2));
%     [~, ord] = sort(x);
%     plot(x(ord), y(ord), 'LineWidth', 2, 'Color', 'r');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;


%% TO FIX: looking at #st vs #t responses 

% compare average responses for each electrode when [t] is preceded by an s
% versus not preceded by an s
vot = diff(TDcons.phnTimes)-1;
st_idx = strcmp(TDcons.precPhn, 's') & strcmp(TDcons.phn, 't');
t_idx = ~strcmp(TDcons.precPhn, 's') & strcmp(TDcons.phn, 't');

figure; 
subplot(1, 2, 1);
histogram(vot(st_idx)); hold on;
histogram(vot(t_idx));

rng(9);
t_idx_res = zeros(1, length(vot));
for i = unique(vot(st_idx))
   possible = ~strcmp(TDcons.precPhn, 's') & strcmp(TDcons.phn, 't') & ...
       vot == i;

   % resample from the #t distribution to match vot of #st distribution   
   t_idx_res(randsample(find(possible), ...
       min(sum(vot(st_idx)==i), sum(possible)))) = 1;
end
t_idx_res = logical(t_idx_res);

subplot(1, 2, 2);
histogram(vot(st_idx)); hold on;
histogram(vot(t_idx_res));

% thresholding for rho of 0.05
vot_encoding_max = vot_encoding(abs(vot_encoding.rho)>0.05, :);
resp = cell(height(vot_encoding_max), 2);

for i = 1:height(vot_encoding_max)
    el = vot_encoding_max{i, 'el'};
    sid = vot_encoding_max{i, 'SID'}{1};
    maxtp = vot_encoding_max{i, 'maxtp'};

    resp{i, 1} = squeeze(TDcons.(sid).resp(el, bef+uint8(maxtp*100), st_idx));
    resp{i, 2} = squeeze(TDcons.(sid).resp(el, bef+uint8(maxtp*100), t_idx_res));
    
    [~, p(i)]=ttest2(resp{i, 1}, resp{i, 2});
    
    if p(i)<0.01
        figure;
        subplot(1, 2, 1);
        boxplot(resp{i, 1}); hold on; 
        boxplot(resp{i, 2}, 'Position', 2);
        xticks([1, 2]);
        xticklabels({'#st', '#t'})

        subplot(1, 2, 2);
        scatter(10:10:60, vot_encoding_max.resp(i, :), 'k', 'filled');
        line(10:10:60, vot_encoding_max.resp(i, :), 'Color', 'k', 'LineWidth', 1.5);
    end
end

% compare average responses for each electrode when [t-d] continuum is 
% preceded by an s versus not preceded by an s

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max;

%% ------------- NEURAL CLASSIFICATION OF VOICING -------------------------
%% LDA for VOT category discrimination

% vot_encoding = vot_encoding(vot_encoding.p<0.01, :);
% single time frame LDA dimensionality reduction
Scons = TDcons;
corpus = 'timit';
Scons.corpus = corpus;
max_comp = 500;
nfolds = 20;
voicing = {{'b', 'g', 'd'}, {'p', 'k', 't'}};
phns = {'p', 'b', 'k', 't','g', 'd'}; % 'k', 't','g', 'd'
% language background
subjs = [1, 2]; %4; %[1, 2]
idx = ismember(Scons.phn, phns) & diff(Scons.phnTimes)<8;

% load in all COT discriminating electrode responses
% create matrix -- el x time (onset:onset+200) x trial
% load in data from pool of subjects (English vs. Spanish monolinguals)
varnames = {'SID', 'el', 'lang', 'neuralresp'};
consresp = array2table(zeros(0,4), 'VariableNames', varnames);
    % n2pe/ps --> phonetic english, phonetic spanish

% aggregate neural data
for s = [eSIDs sSIDs] % bSIDs %, 'EC100', 'EC105', 'EC163', 'EC214'
    SID = s{1};  

    % 1 - spanish, 2 - english, 3 - bilingual
    ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs}));
    
    % use neural window from onset to 300 ms after onset & baseline    
    tps = bef+10:bef+30;
    startp = bef+10;
    %tps = 1:7;
    ctr = 1;
    els = vot_encoding.el(strcmp(vot_encoding.SID, SID));
    
    if ~isempty(els)
        allr = nan(length(els), length(tps), length(Scons.phn(idx)));
        for e = els'
            allr(ctr, :, :) = squeeze(Scons.(SID).resp(e, tps, idx)); 
            
            % add to line to table
            t2 = table({SID}, e, ls, {allr(ctr, :, :)}, ...
                'VariableNames', varnames);
            consresp  = [consresp ; t2];
            
             % turn into average phoneme pair representation
            ctr = ctr+1;        
        end
    end
    clear allr e ctr t2 
end

AUC = nan(length(subjs), nfolds);
acc = cell(length(subjs), 1);
% electrode weighting through the logistic regression and PCA 
comp_weighted = cell(length(subjs));
elecs = cell(length(subjs));
scores = cell(length(subjs));
yidx = cell(length(subjs), 1);
for ls = subjs
    % using more than simple p/b contrast
    y = Scons.phn(idx);
    % 1 is voiced, 2 is unvoiced
    y = arrayfun(@(x) ismember(x, voicing{1}), y)+1; 
    elshape = reshape(repmat(find(consresp.lang==ls), 1, length(tps)), ...
        sum(consresp.lang==ls)*length(tps), []);
    
    % run pca
    % A = reshape(cell2mat([consresp.('neuralresp')]), ...
    %     height(consresp)*length(tps), []);
    A = squeeze(cell2mat([consresp.('neuralresp')(consresp.lang==ls)]));
    [X_mod, nanrow, nancol, y] = makeDataMatrix(A, y, ...
        consresp.SID(consresp.lang==ls));
    tmp = find(idx);
    yidx{ls} = tmp(~nancol);
    elshape(nanrow) = [];
    elecs(ls) = {unique(elshape)};

    [~, ~, AUC(ls, :), pcaX, tmp_score, acc{ls}, weights] = logistic(X_mod', ...
        y'-1, 1, [], tps, nfolds);
    scores{ls} = nan(sum(idx), 1);
    scores{ls}(~nancol) = tmp_score(:, 1);
    comp_weighted(ls) = {weights};
    
    % calculate silhouette score
    disp('---------------------- CONS NEUR LDA -------------------------------');
    disp(['accuracy = ' num2str(mean(acc{ls}))]);
    disp(['AUC = ' num2str(mean(AUC(ls, :)))]);
    disp(['n trials = ' num2str(length(y))]);
end

decode_details = struct();
decode_details.acc = acc; % accuracy for both subject group types
decode_details.startp = startp;
decode_details.tps = tps;
decode_details.AUC = AUC;
decode_details.elecs = elecs;
decode_details.encoding = consresp; % for threshold 
decode_details.score = scores;

% 1 is voiced, 2 is unvoiced
decode_details.y = arrayfun(@(x) ismember(x, voicing{1}), Scons.phn(idx))+1;
decode_details.yidx = yidx;
decode_details.weights = comp_weighted;

if isscalar(subjs) && subjs == 3
    subj = 'mandarinmono';
elseif isscalar(subjs) && subjs == 4
    subj = 'bilingual';
else
    subj = 'monolingual';
end

subtype = '';
filename = [corpus '_voicing_decode_' subj '_' subtype '.mat']; % dimex filename
save([datapath 'ecog_decode/VOT/' filename], 'decode_details');

disp(['saving ..... ' filename])

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons*;

%% plotting decoding results

Scons = TDcons;
corpus = 'timit';
lss = [1, 2, 4];
type = '';

figure('Position', [100, 300, 150, 300], 'renderer', 'painters');
position = [1, 2, 0, 3];
for ls = lss

    if ls == 3
        subj = 'mandarinmono';
    elseif ls == 4
        subj = 'bilingual';
    else
        subj = 'monolingual';
    end
    
    filename =  [corpus '_voicing_decode_' subj '_.mat']; 
    load([datapath 'ecog_decode/VOT/' filename], 'decode_details');

    styles = {':', '-', '', '-'};
%     subjs = {'Spanish', 'English'};
%     startp = decode_details.startp;
    accls(ls, :) = squeeze(decode_details.AUC(ls, :, :));
end


colors = getColors(1);
for ls = lss
    boxplot(accls(ls, :)*100, 'Position', position(ls), 'Color', colors(ls, :), ...
            'BoxStyle', 'filled'); hold on;
    xlabel('Subject Group');
end

%formatting
xticklabels({'Spanish', 'English', 'Bilingual' });
xticks([1 2 3]);
xlim([0.5 3.5]);
legend('off');

pos = [1, 3; 2, 3; 1, 2];
combo = [1, 4; 2, 4; 1, 2];
for i = 1:3
    
    [p, h] = ranksum(accls(combo(i, 1), :), accls(combo(i, 2), :));
    disp(p)
    if ~isempty(getSigStr(p))
        plot([pos(i, 1)-0.02 pos(i, 2)+0.02], ...
            [70+5*i 70+5*i], 'LineWidth', 1.5, 'Color', 'k');
        text(mean(pos(i, :))-0.25, 70+5*i+1, getSigStr(p), 'FontSize', 15);
    end
end

% for sliding window
% labels = startp+round(decode_details.timing(1)/2)-startp(1);
% xticklabels(split(num2str(labels*10)));
yline(0.5, 'Color', 'k');
yticks();
ylim([30 70]);
yticks([30, 70])
ylabel('AUC')    
set(gca, 'FontSize', 15);
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% plotting decoding native vs. nonnative


corp = {'dimex', 'timit'};
Sconss = {Dcons, TDcons};
native = [];
nonnative=[];
for s = 1:2

    Scons = Sconss{s};
    corpus = corp{s};
    lss = [1, 2];
    type = '';
    
    figure('Position', [100, 300, 150, 300], 'renderer', 'painters');
    position = [1, 2, 0, 3];
    for ls = lss
    
        if ls == 3
            subj = 'mandarinmono';
        elseif ls == 4
            subj = 'bilingual';
        else
            subj = 'monolingual';
        end
        
        filename =  [corpus '_voicing_decode_' subj '_.mat']; 
        load([datapath 'ecog_decode/VOT/' filename], 'decode_details');
    
        styles = {':', '-', '', '-'};
    %     subjs = {'Spanish', 'English'};
    %     startp = decode_details.startp;
        accls(ls, :) = squeeze(decode_details.AUC(ls, :, :));
    end

    if s==ls && ismember(ls, [1, 2])
            native = [native; accls(ls, :)];
    elseif ismember(ls, [1, 2])
            nonnative = [nonnative; accls(ls, :)];
    end  
    
    colors = getColors(1);
    for ls = lss
        boxplot(accls(ls, :)*100, 'Position', position(ls), 'Color', colors(ls, :), ...
                'BoxStyle', 'filled'); hold on;
        xlabel('Subject Group');
    end
    
    %formatting
    xticklabels({'Spanish', 'English', 'Bilingual' });
    xticks([1 2 3]);
    xlim([0.5 3.5]);
    legend('off');
    
    pos = [1, 3; 2, 3; 1, 2];
    combo = [1, 4; 2, 4; 1, 2];

    combo = [1, 2];
    pos = [1, 2];
    for i = 1:1
        
        [p, h] = ranksum(accls(combo(i, 1), :), accls(combo(i, 2), :));
        disp(p)
        plot([pos(i, 1)-0.02 pos(i, 2)+0.02], ...
            [70+5*i 70+5*i], 'LineWidth', 1.5, 'Color', 'k');
        text(mean(pos(i, :))-0.25, 70+5*i+1, getSigStr(p, 2), 'FontSize', 15);
    end
    
    % for sliding window
    % labels = startp+round(decode_details.timing(1)/2)-startp(1);
    % xticklabels(split(num2str(labels*10)));
    yline(0.5, 'Color', 'k');
    yticks();
    ylim([30 70]);
    yticks([30, 70])
    ylabel('AUC')    
    set(gca, 'FontSize', 15);
    box off;
end

figure;
boxchart(ones(numel(native), 1), native(:), 'BoxFaceColor', [0.5 0.5 0.5]); hold on;
boxchart(ones(numel(nonnative), 1)*2, nonnative(:), 'BoxFaceColor', [0.5 0.5 0.5]);
ylabel('AUC');
xlabel('Monolingual Subject Type');
xticks([1, 2]);
ylim([0.5 0.9]);
yticks([0.3 0.8])
set(gca, 'FontSize', 15);
xticklabels({'native', 'non-native'});

[p, ~] = ranksum(native(:), nonnative(:));
plot([1 2], ...
    [0.67 0.67], 'LineWidth', 1.5, 'Color', 'k'); hold on;
text(mean(pos(i, :))-0.25, 0.7, getSigStr(p, 2), 'FontSize', 15);


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd*;

%% ------------------------- FUNCTIONS ------------------------------------

function [colors] = getColors(type, numcol)

    if nargin<2, numcol = 5; end
    switch type
        case 1 % language type
            % spanish, english, mandarin
            colors = brewermap(4, 'Dark2');
        case 2 % primary feature encoding map
            colors = brewermap(3, 'Dark2');
        case 3 % vot map
            colmap = mkcolormap([1 0 0],[0.8 0 0.8],[0 0 1]);
            colors = flipud(colmap(round(linspace(1, 200, numcol)), :));
    end
end

function plotVOTFreq(TDcons, phns, syllOns)
    % psychometric function based on transcribed category information
    % phns should be grouped into two cells, ordered voiced and un-voiced
    % e.g. {{'b', 'g'}, {'p', 'k'}}
    
    % calculate VOT
    vot = diff(TDcons.phnTimes)-1;
    maxvot = 5; % 50 ms
    if ~syllOns
        phnidx = ismember(TDcons.phn, phns{1}) | ismember(TDcons.phn, phns{2}) ...
            & vot<=maxvot;
    else
        phnidx = ismember(TDcons.phn, phns{1}) | ismember(TDcons.phn, phns{2}) ...
            & vot<=maxvot & TDcons.syllOns;
    end
    
    reps = 50;
    y = nan(maxvot, reps);
    for v = 1:maxvot        
        % sample with replacement (50 times)
        for s = 1:reps
            % sample idx
            sidx = randsample(find(phnidx), round(sum(phnidx)*0.3));

            % percent voiced from sample
            unvoiced = ismember(TDcons.phn(sidx), phns{2});
            y(v, s) = sum(vot(sidx)==v & unvoiced)/sum(vot(sidx)==v);
        end
    end

    x = (1:maxvot)*10;
    errorbar(x, mean(y, 2), std(y, [], 2, 'omitnan'), 'o', 'CapSize',0, 'Color', 'k'); 
    hold on;
    scatter(x,mean(y, 2), 60, getColors(3, maxvot), 'filled');    
    
    % fit sigmoidal line
    [param,~]=sigm_fit(x,mean(y, 2));
    fsigm = @(param,xval) param(1)+...
        (param(2)-param(1))./(1+10.^((param(3)-xval)*param(4)));

    sigmx = 10:1:50;
    plot(sigmx, fsigm(param, sigmx), 'LineWidth', 2,'LineStyle', ...
        '--',  'Color', 'k');

    ylabel('% labeled unvoiced');
    yticks([0 0.5 1]);
    yticklabels({'0', '50', '100'});
    xticks(10:10:50);
    xlabel('VOT (ms)');
    set(gca, 'FontSize', 15);
end

function plotVOTerp(TDcons, SID, els, ~, phns, syllOns)
    % calculate VOT
    vot = diff(TDcons.phnTimes);
    if ~syllOns
        phnidx = ismember(TDcons.phn, phns) & vot<7;
    else
        phnidx = ismember(TDcons.phn, phns) & vot<7 & TDcons.syllOns;
    end

    desmat.names = {'10', '20', '30', '40', '50'};
    desmat.condition = vot(phnidx)-1;

    disp(crosstab(desmat.condition));
    plotResp(TDcons.(SID).resp(:, :, phnidx), desmat,els, SID, ...
        [], flipud(getColorsCrossComp(5)), [], [], 1); % getColors(3, length(desmat.names))

end

% remove all non-prevocalic consonants in Dcons struct
function [Dcons] = prevocCons(Dcons, features)
    clostr = 'cl'; % not true for DIMEx

    % modify such that preceding phoneme is the phoneme before closure
    % add word onset that readjusts to ignore the closure label
    for i = find(ismember(Dcons.phn, {'p', 'b', 'k', 'g', 'd', 't'}))

        currPhn = Dcons.phn(i);
        % preceding is closure
        if i>1
            precPhn = strcmp(Dcons.phn(i-1), join([currPhn clostr], '')); 
            % preceding is in the same sentence and temporally adjacent
            sameSent = Dcons.sent(i-1)==Dcons.sent(i);
            adj = abs(Dcons.phnTimes(2, i-1)-Dcons.phnTimes(1, i))<3;
            
            % if sentence is the same and the onset 
            if precPhn && sameSent && adj
                Dcons.precPhn(i) = Dcons.precPhn(i-1);
                Dcons.wordOns(i) = Dcons.wordOns(i-1);
                Dcons.syllOns(i) = Dcons.syllOns(i-1);
                Dcons.sentOns(i) = Dcons.sentOns(i-1);

                % stress
                Dcons.precPhnStress(i) = Dcons.precPhnStress(i-1);
            end
        end
    end

    vows = features.sonorant;
    phnidx = ismember(Dcons.nextPhn, vows);

    numel = length(Dcons.phn);
    for field = fieldnames(Dcons)'
        if size(Dcons.(field{1}), 2)==numel
            Dcons.(field{1}) = Dcons.(field{1})(:, phnidx);
        elseif size(Dcons.(field{1}), 3)==numel
            Dcons.(field{1}) = Dcons.(field{1})(:, :, phnidx);
        elseif startsWith(field, 'EC')
            Dcons.(field{1}).resp = Dcons.(field{1}).resp(:, :, phnidx);
        end
    end
end

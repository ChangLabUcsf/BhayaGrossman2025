%% Set up
% Ilina Bhaya-Grossman
% 05.21.2025
out_crosscomp_startup;
if ~exist('Dcons', 'var')
    load("data/ExtendedFigures/ExtFigure5_DIMEXCons.mat");
end

if ~exist('TDcons', 'var')
    load("data/ExtendedFigures/ExtFigure5_TIMITCons.mat");
end

if ~exist('vot_encoding', 'var')
    load("data/ExtendedFigures/ExtFigure5_NeuralVOT.mat");
end

if ~exist('imgall', 'var')
    load([datapath '/Figure1/Figure1_ImgData.mat'], 'imgall');
end

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons* phnnames vs vidx;

%% A - VOT distribution in timit

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

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons* phnnames vs vidx;

%% UNUSED - Additional information...

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
plotVOTFreq(TDcons, {{'b'}, {'p'}}, 0);

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons* phnnames vs vidx;

%% UNUSED - Matrix with number of trials per category for each VOT bin

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
clim([0 100]);

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
clim([0 250]);

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons* phnnames vs vidx;

%% B - Plotting example VOT erps

% meeting slides
% V+: Spanish - EC172 el219, English - EC222 el125 
% V-: Spanish - EC152 el101, English - EC186 el231

% Quals slides
% V-: Spanish - EC105 el 185, English - EC195 el137
% V+: Spanish - EC100 el70, English - EC195 el217

phns = {'g', 'k', 'b', 'p' , 'd', 't'}; % , ,  
rmpath('../util/shadederror');
% SID = 'EC266'; % 'EC222'
% el = 134;, 122
% SID = 'EC100'; % 'EC222'
% el = 134;

% SID = 'EC195'; % 'EC222'
% el = 217;

SID = 'EC186'; % 'EC222'
el = 239;

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


%% C - Calculate PSI to select plosive electrodes and spearman rho for VOT
% the significance testing changes depending on which subjects are included
% and which threshold is used! Please note when reporting results.
vot_max = 8;
vot_encoding = vot_encoding(vot_encoding.p<0.05, :); 

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

for t = 1:2
    span = stats_y{1, t};
    eng = stats_y{2, t};

    for v = 1:vot_max-2
        [p, ~] = ranksum(span(:, v), eng(:, v));
        disp(p)
        if p < 0.05
            subplot(1, 2, t)
            text(v*10, 1.05, getSigStr(p, 1), 'FontSize', 20);
        end
    end
end
legend({'Spanish', 'English', 'Bilingual'});

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* phns vot_max phnnames vs vidx;
%% UNUSED - Plot tuning curve by category type

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

%% D - Proportion of V+ electrodes across participants

figure;
subplot(1, 2, 1)
b=bar([1, 2], [crosstab(vot_encoding.rho(vot_encoding.ls==2)<0)./length(vot_encoding.rho(vot_encoding.ls==2)), ...
    crosstab(vot_encoding.rho(vot_encoding.ls==1)<0)./length(vot_encoding.rho(vot_encoding.ls==1))]');
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
    rng(1);
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

clearvars -except *all* subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh *cons* phnnames vs vidx vot_max;

%% E - MNI elecs for direction of VOT encoding

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
    
    [mni_lh] = plotMNIElec(unique(vot_encoding.SID), desel, 'lh', 0, 1, imgall);

    title(ls)
    
    [mni_rh] = plotMNIElec(unique(vot_encoding.SID), desel, 'rh', 0, 1, imgall);
    title(ls)
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

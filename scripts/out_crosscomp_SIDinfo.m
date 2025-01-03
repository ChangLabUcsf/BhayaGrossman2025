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
SIDs = [sSIDs, eSIDs]; % , {'HS11', 'HS9', 'HS10'}

%% -------------------- PROFICIENCY ---------------------------------------
%% -------------------- Hemisphere stats ----------------------------------

[imgall] = load_allimgdata;

hemis = cellfun(@(x) imgall.(x).hemi, sSIDs, 'UniformOutput', false);
disp('Spanish mono')
disp(['Left: ' num2str(sum(strcmp(hemis, 'lh'))) ', ' ...
        'Right: ' num2str(sum(strcmp(hemis, 'rh'))), ', ' ...
        'Total: ' num2str(length(hemis))]);

hemis = cellfun(@(x) imgall.(x).hemi, eSIDs, 'UniformOutput', false);
disp('English mono')
disp(['Left: ' num2str(sum(strcmp(hemis, 'lh'))) ', ' ...
        'Right: ' num2str(sum(strcmp(hemis, 'rh'))), ', ' ...
        'Total: ' num2str(length(hemis))]);

hemis = cellfun(@(x) imgall.(x).hemi, bSIDs, 'UniformOutput', false);
disp('Bilingual')
disp(['Left: ' num2str(sum(strcmp(hemis, 'lh'))) ', ' ...
        'Right: ' num2str(sum(strcmp(hemis, 'rh'))), ', ' ...
        'Total: ' num2str(length(hemis))]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% -------------------- Monolingual proficiency ------------------------ %%

varnames = {'SID', 'ls', 'spanish', 'english', 'mandarin', 'other'};
proficiency =  array2table(zeros(0, 6), 'VariableNames', varnames);

% All Spanish subjects with proficiency in all tested languages
% missing language info for: EC100, EC105, EC152
sSIDs = {'EC172', 'EC163', 'EC105', 'EC100', 'EC203', 'EC214', 'EC152', 'EC225', 'EC252'}';
sprof = {5, 5, 5, 5, 5, 5, 5, 5, 5}';
eprof = {0, 0, 0, 0, 1, 0, 0, 0, 0}'; 
mprof = {0, 0, 0, 0, 0, 0, 0, 0, 0}'; 
other = {0, 0, 0, 0, 0, 0, 0, 0, 0}'; 
lang = repmat({'s'}, length(sSIDs), 1);
proficiency = [proficiency; [sSIDs, lang, sprof, eprof, mprof, other] ];

% All English subjects with proficiency in all tested languages
% double checked
eSIDs={'EC196', 'EC195', 'EC183', 'EC212', 'EC186', 'EC219', 'EC221', 'EC222', 'EC235', 'EC242'}';
eprof = {5; 5; 5; 5; 5; 5; 5; 5; 5; 5};
sprof = {1; 1; 0; 0; 3; 2; 2; 1; 3; 0};
mprof = {0; 0; 0; 0; 0; 0; 0; 0; 0; 0};
other = {1; 1; 0; 2; 0; 0; 0; 0; 1; 0}; 
lang = repmat({'e'}, length(eSIDs), 1);
proficiency = [proficiency; [eSIDs, lang, sprof, eprof, mprof, other] ];

% All Mandarin subjects with proficiency in all tested languages
mSIDs={'HS11', 'HS10', 'HS9', 'HS8'}'; % 'HS14',
mprof = {5; 5; 5; 5}; % 5
eprof = {4; 0; 0; 0}; % 2 other language proficiency
sprof = {0; 0; 0; 0}; % 0
other = {0; 0; 0; 0}; % 0
lang = repmat({'m'}, length(mSIDs), 1);
proficiency = [proficiency; [mSIDs, lang, sprof, eprof, mprof, other] ];

% Plotting proficiencies
colors = brewermap(3, 'Dark2');
lss = {'s', 'e', 'm'};
figure('Position',[500, 80, 200, 600])
for i = 1:3
    subplot(4, 1, i);
    ls = lss{i};
    idx = find(strcmp(proficiency.ls, ls));
    numid = length(idx);
    rs = (rand(numid, 1)*0.5)-0.25;
    x = [zeros(numid, 1)+rs ones(numid, 1)+rs ...
        ones(numid, 1)*2+rs ones(numid, 1)*3+rs];
        
    y = [proficiency.spanish(idx) ...
        proficiency.english(idx) ...
        proficiency.mandarin(idx) ...
        proficiency.other(idx)];
    
    scatter(x(:), y(:), ...
        65, colors(i, :), 'filled', 'MarkerFaceAlpha', 0.7); hold on;
    
    plot(x', y', 'LineWidth', 0.5, ...
        'Color', colors(i, :), 'HandleVisibility', 'off');

    % Formatting 
    xlim([-0.5 3.5])
    ylim([-0.5 5.5]);
    ylabel('Proficiency');
    xticks([]);
    
    text(3, 4.5, ['n=' num2str(numid)]);
    set(gca, 'FontSize', 13)
end
xticks(0:3);
xticklabels({'Spanish', 'English', 'Mandarin', 'Other'});
xtickangle(45);
xlabel('Language');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons;

%% -------------------- Bilingual proficiency -------------------------- %%

figure;
[proficiency] = getBilingProf();

% Plot proficiency
% subplot(4, 1, 4);
tmp = getColorsCrossComp(1);
colors = [tmp(1:2, :); 0.5 0.5 0.5];
L1s = {'S', 'E', 'O'};
for i = 1:3
    idx = find(strcmp(proficiency.L1, L1s{i}));
    numid = length(idx);
    rs = (rand(numid, 1)*0.24)-0.12;
    x = [zeros(numid, 1)+rs ones(numid, 1)+rs ...
        ones(numid, 1)*2+rs]; %  
        %ones(numid, 1)*3+rs
        
    y = [proficiency.span_prof(idx) ...
        proficiency.eng_prof(idx) ...
        proficiency.other_prof(idx)]; %         proficiency.mand_prof(idx) ...
    
    scatter(x(:), y(:), 65, ...
        colors(i, :), 'filled', 'MarkerFaceAlpha', 0.7); hold on;
    
    plot(x', y', 'LineWidth', 0.5, ...
        'Color', colors(i, :), 'HandleVisibility', 'off');
    
end
text(3, 4.5, ['n=' num2str(numid)]);

% Formatting
xlim([-0.25 2.25])
ylim([-0.5 5.5]);
ylabel('Proficiency');


xticks(0:3);
xtickangle(0);
xticklabels({'Spanish', 'English', 'Other'}); % 'Mandarin', 
yticks([1, 3, 5]);
yticklabels({'single words', 'fluent', 'native'});
ytickangle(40);

leg = legend({'Spanish', 'English', 'Other'});
title(leg,'L1')
xlabel('Language');
set(gca, 'FontSize', 13)

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons;

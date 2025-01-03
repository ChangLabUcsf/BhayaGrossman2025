%% set up
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath brewer
addpath(genpath('util'))
zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;
addpath(genpath(datapath))

bef=20;
aft=50;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs] = getSIDinfo();
% eSIDs = {'EC219' 'EC183' 'EC212' 'EC222', 'EC195', 'EC186', 'EC235'}; % 	
% sSIDs = {'EC172' 'EC163' 'EC105' 'EC214', 'EC100', 'EC203', 'EC252'}; % 	
SIDs=eSIDs;

% load all acoustic data
addpath(genpath('util'))
Dvow = loadDD('dimex', bef, aft, {}, datapath, 0, [], 0);
addpath(genpath('util'))
TDvow = loadDD('timit', bef, aft,{}, datapath, 0, [], 0);
Dvow.corpus = 'dimex'; 
TDvow.corpus = 'timit';

% load in all beta model versions
% 
timit_vow = {'aa', 'ae', 'ao', 'ah', 'ey', 'eh', 'ih', 'iy', 'ow'};
timit_details = load('stim_info/out_sentence_details_timit_all_loudness.mat');
dimex_details = load('stim_info/out_sentence_details_dimex_all_loudness.mat');
tps = 50:55;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;


%% load in phonetic feature RDMs
arpadef = {'aa', 'ae', 'ah', 'ao', 'aw', 'ax', 'ay', 'eh', ...
    'er', 'ey', 'ih', 'ix', 'iy', 'ow', 'oy', 'uh', 'uw'};
arpaphn = vowelfeat();

% phonetic encoding
pe = nan(length(TDvow.vowel), size(arpaphn.features, 2)-2);
ctr = 1;
for i = TDvow.vowel    
    if ismember(i, arpaphn.phnnames)
        pe(ctr, :) = table2array(arpaphn.features(...
            find(strcmp(arpaphn.features.name, i)), 3:end));  
    end
    ctr = ctr + 1;
end
TDvow.phoneticRDM{1} = corr(pe', 'Type', 'Spearman');
TDvow.phoneticRDM{2} = corr(pe(:, 1:4)', 'Type', 'Spearman');

% spectral encoding
filename = 'mds_timit_spectro50_stress.mat';
if isfile(filename)
    sp = load(filename);
else
    allaud = [TDvow.aud];
    A = reshape(allaud(:, tps, :), [], size(allaud, 3));    
    sp = corr(A, 'Type', 'Spearman');
    
    % dtw
    
end
TDvow.spectralRDM = sp;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn;

%% visualizing phonetic and spectral RDMs
figure
subplot(1, 3, 1);
load('util/arpabet2ipa.mat');
ipanames = ipa(cellfun(@(i) find(strcmp(arpabet, i)), arpaphn.phnnames));

imagesc(corr(table2array(arpaphn.features(:, 3:end))'));
yticks(1:13); xticks(1:13);
yticklabels(ipanames);
xticklabels(ipanames);
set(gca, 'FontSize', 15);
title({'Vowel Distances'}, {'by English features'})

arpaphn = vowelfeat();
subplot(1, 3, 2);
imagesc(corr(table2array(arpaphn.features(:, 3:6))'));
yticks(1:13); xticks(1:13);
yticklabels(ipanames);
xticklabels(ipanames);
title({'Vowel Distances'}, {'by Spanish features'})
set(gca, 'FontSize', 15);

sp = nan(13, 13);
for i = 1:length(arpaphn.phnnames)
%     sp(i, i) = 1;
    for j = i:length(arpaphn.phnnames)
        iid = strcmp(TDvow.vowel,arpaphn.phnnames(i));
        jid = strcmp(TDvow.vowel,arpaphn.phnnames(j));

        allaud = [TDvow.aud];
        idx = TDvow.stress == 1 & ismember(TDvow.vowel, arpaphn.phnnames);
        alli = reshape(allaud(:, tps, idx&iid), [], sum(idx&iid));  
        allj = reshape(allaud(:, tps, idx&jid), [], sum(idx&jid));
        
        sp(i, j) = mean(corr(alli, allj), [1, 2]);
        sp(j, i) = mean(corr(alli, allj), [1, 2]);       
    end
end

subplot(1, 3, 3);
imagesc(sp);
yticks(1:13); xticks(1:13);
yticklabels(ipanames);
xticklabels(ipanames);
title({'Vowel Distances'}, {'by Spectral features'})
set(gca, 'FontSize', 15);

colormap(flipud(brewermap(30, 'Spectral')));

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn;




%% ---------------------- Bilingual analysis ------------------------------
%% scatter of correlation of STRF betas between English and Spanish corpora
modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_aud'};
SIDs = eSIDs;

figure;
bimola = struct();
for s = SIDs
    SID = s{1};
    corpusStrf{1}=loadMultModelStrf(SID, modelname, 'dimex', datapath);   
    corpusStrf{2}=loadMultModelStrf(SID, modelname, 'timit', datapath);  
    
    if ~isempty(corpusStrf{1}{1}) && ~isempty(corpusStrf{2}{1})
        els = 1:min(length(corpusStrf{1}{1}.Els), length(corpusStrf{2}{1}.Els));
        bimola.(SID).corr = arrayfun(@(e) corr(reshape(corpusStrf{1}{1}.meanStrf(:, :, e), ...
            [], 1), reshape(corpusStrf{2}{1}.meanStrf(:, :, e), [], 1)), els);       
        bimola.(SID).rsq(1, els) = corpusStrf{1}{1}.meanTestR(els).^2;
        bimola.(SID).rsq(2, els) = corpusStrf{2}{1}.meanTestR(els).^2;
        clear corpusStrf

        idx = all(bimola.(SID).rsq>0.01, 1); % els(idx), repmat(str2num(SID(3:end)), 1, sum(idx))
        scatter3(bimola.(SID).rsq(1, idx), bimola.(SID).rsq(2, idx),repmat(str2num(SID(3:end)), 1, sum(idx)) , 55, ...
            bimola.(SID).corr(idx), 'filled', 'MarkerFaceAlpha', 0.7, ...
            'MarkerEdgeColor', 'k'); hold on;
        view(2);
    else
        disp(['Missing STRF for ' SID]);
    end
end

ylim([0 0.35]);
xlim([0 0.35]);
yticks(0:0.1:0.3);
xticks(0:0.1:0.3);
h = refline(1, 0);
h.Color = 'k';
caxis([0 1]);
colormap(flipud(brewermap(20, 'PuOr')));

ylabel('TIMIT TRF R^2');
xlabel('DIMEX TRF R^2');
set(gca, 'FontSize', 15);

cbh = colorbar();
ylabel(cbh, 'Beta Correlation');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn;

% comparison of correlation between English, Spanish, Bilingual

modelname = {'onset_maxDtL_maxDtLOnset_vowelOnset_aud'};

varnames = {'SID', 'el', 'lang', 'corr', 'srsq', 'ersq'};
bimola = array2table(zeros(0,6), 'VariableNames', varnames);
for SIDs = {sSIDs, eSIDs, bSIDs}
    for s = SIDs{1}
        SID = s{1};
        corpusStrf{1}=loadMultModelStrf(SID, modelname, 'dimex', datapath);   
        corpusStrf{2}=loadMultModelStrf(SID, modelname, 'timit', datapath);  

        if ~isempty(corpusStrf{1}{1}) && ~isempty(corpusStrf{2}{1})
            els = 1:min(length(corpusStrf{1}{1}.Els), length(corpusStrf{2}{1}.Els));
            
            rsq = nan(2, length(els));
            correl = arrayfun(@(e) corr(reshape(corpusStrf{1}{1}.meanStrf(:, :, e), ...
                [], 1), reshape(corpusStrf{2}{1}.meanStrf(:, :, e), [], 1)), els);       
            rsq(1, els) = corpusStrf{1}{1}.meanTestR(els).^2;
            rsq(2, els) = corpusStrf{2}{1}.meanTestR(els).^2;

            idx = all(rsq>0.1, 1);
            ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
            
            t2 = table(repmat({SID}, sum(idx), 1), els(idx)', repmat(ls, sum(idx), 1), ...
                correl(idx)', rsq(1, idx)', rsq(2, idx)', 'VariableNames', varnames);
            bimola = [bimola; t2];            
        else
            disp(['Missing STRF for ' SID])
        end
    end
end

figure;
%boxplot(bimola.corr, bimola.lang-1, );
cols = [0.6000 0.4392 0.6706; 0.3529 0.6824 0.3804; 0.8784 0.5098 0.0784];
violin({bimola.corr(bimola.lang==1), ...
    bimola.corr(bimola.lang==2), bimola.corr(bimola.lang==3)}, 'facecolor', ...
    cols, 'medc', []);
ylim([-0.01 1]);
xticks(1:3);
xticklabels({'Spanish', 'English', 'Bilingual'});
set(gca, 'FontSize', 15);
legend('off');
yticks(0:0.5:1);
ylabel({'Correlation between','TRF weights'});

SID = 'EC129';
corpusStrf{1}=loadMultModelStrf(SID, modelname, 'dimex', datapath);
corpusStrf{2}=loadMultModelStrf(SID, modelname, 'timit', datapath);
figure; 
subplot(1, 2, 1); 
imagesc(corpusStrf{1}{1}.meanStrf(:, :, 183)); 
yticks([5 84]);
yticklabels({'1' ,'8'});
xticks([1 60]); xticklabels({'0', '-0.6'}); ylabel('Time (s)')
ylabel('Frequency (kHz)');
set(gca, 'FontSize', 15);

subplot(1, 2, 2);
imagesc(corpusStrf{2}{1}.meanStrf(:, :, 183));
yticks([5 84])
yticklabels({'1' ,'8'});
xticks([1 60]); xticklabels({'0', '-0.6'}); ylabel('Time (s)')
ylabel('Frequency (kHz)');
set(gca, 'FontSize', 15);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn bimola;

%%

figure;
imagesc(flipud(TDvow.aud(:, 25:75, 2))); 
colormap(flipud(bone)); 
xline(25, 'LineWidth', 4); xticks([1 25 50]); 
xticklabels({'-0.25', '0', '0.25'}); 
yticks([1 80]); 
yticklabels({'8' '1'});
set(gca, 'FontSize', 15);
ylabel('Frequency (kHz)');
xlabel('Time (s)');

figure
imagesc(flipud(TDvow.aud(:, 25:75, 12))); 
colormap(flipud(bone)); 
xline(25, 'LineWidth', 4); xticks([1 25 50]); 
xticklabels({'-0.25', '0', '0.25'}); 
yticks([1 80]); 
yticklabels({'8' '1'});
set(gca, 'FontSize', 15);
ylabel('Frequency (kHz)');
xlabel('Time (s)');

%% scratch work

% /i/ vs. /I/, /aa/ vs. /ao/, /uh/ vs. /uw/
vowpair = {{'ih', 'iy'}, {'aa', 'ao'}, {'uh', 'uw'}};
% eSIDs(9) = [];
TDvow = addtoDD(TDvow, 'timit', bef, aft, [sSIDs eSIDs]);
load('out_elecs_voweltypeftest_bychan.mat')
varnames = {'SID', 'el', 'lang', cell2mat(vowpair{1}), cell2mat(vowpair{2}), ...
    cell2mat(vowpair{3}), 'fval'};
sensestruct = array2table(zeros(0,7), 'VariableNames', varnames);
for s = {'EC222'} %[sSIDs eSIDs]
    SID = s{1};
    idx = TDvow.stress==1&TDvow.meanf0<170;
    ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
    for el = allidx.(SID)
        % add for each electrode significant time points (& magnitude?) for
        % which there is a difference between the pairs
        tps = nan(3, 1);
        for v = 1:length(vowpair)
            grp{1} = squeeze(TDvow.(SID).resp(el, :, idx&strcmp(TDvow.vowel, vowpair{v}{1})));
            grp{2} = squeeze(TDvow.(SID).resp(el, :, idx&strcmp(TDvow.vowel, vowpair{v}{2})));
            
            debug = 1;            
            tps(v) = tpttest(grp, 0.05, debug);
        end     
        t2 = table({SID}, el, ls, tps(1), tps(2), tps(3), fvals.(SID)(el), ...
            'VariableNames', varnames);
        sensestruct = [sensestruct ; t2];
    end
end

figure;
subplot(1, 3, 1);
sidx = sensestruct.lang==1;
eidx = sensestruct.lang==2;
y = ksdensity(sensestruct.ihiy(sidx)); 
plot(y, 'LineWidth', 2, 'Color', 'r'); hold on;
y = ksdensity(sensestruct.ihiy(eidx)); 
plot(y, 'LineWidth', 2, 'Color', 'b');
% scatter(sensestruct.ihiy(sidx)-0.5+rand(sum(sidx), 1)*0.5, ...
%     sensestruct.fval(sidx), 35, 'filled', 'MarkerFaceAlpha', 0.5); hold on;
% scatter(sensestruct.ihiy(eidx)-0.5+rand(sum(eidx), 1)*0.5, ...
%     sensestruct.fval(eidx), 35, 'filled', 'MarkerFaceAlpha', 0.5);

subplot(1, 3, 2);
y = ksdensity(sensestruct.aaao(sidx)); 
plot(y, 'LineWidth', 2, 'Color', 'r'); hold on;
y = ksdensity(sensestruct.aaao(eidx)); 
plot(y, 'LineWidth', 2, 'Color', 'b');
% scatter(sensestruct.aaao(sidx)-0.5+rand(sum(sidx), 1)*0.5, ...
%     sensestruct.fval(sidx), 35, 'filled', 'MarkerFaceAlpha', 0.5); hold on;
% scatter(sensestruct.aaao(eidx)-0.5+rand(sum(eidx), 1)*0.5, ...
%     sensestruct.fval(eidx), 35, 'filled', 'MarkerFaceAlpha', 0.5);

subplot(1, 3, 3);
y = ksdensity(sensestruct.uhuw(sidx)); 
plot(y, 'LineWidth', 2, 'Color', 'r'); hold on;
y = ksdensity(sensestruct.uhuw(eidx)); 
plot(y, 'LineWidth', 2, 'Color', 'b');
% scatter(sensestruct.uhuw(sidx)-0.5+rand(sum(sidx), 1)*0.5, ...
%     sensestruct.fval(sidx), 35, 'filled', 'MarkerFaceAlpha', 0.5); hold on;
% scatter(sensestruct.uhuw(eidx)-0.5+rand(sum(eidx), 1)*0.5, ...
%     sensestruct.fval(eidx), 35, 'filled', 'MarkerFaceAlpha', 0.5);
legend({'spanish', 'english'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn bimola *struct;




%% -------------------------- RSA on Vowel sounds -------------------------
%% PCA on spectrogram of DIMEX and TIMIT vowels

figure;
ctr = 1;
tps = 50:70;
Svows = {Dvow, TDvow};
for c = {'dimex', 'timit'}
    corpus = c{1};

    addpath(genpath(datapath))
    
    % Load the data for the current corpus
    Svow = Svows{ctr};

    % Construct the filename for saving the results
    filename = ['mds_' corpus '_spectro50_stress.mat'];
    
    if isfile(filename)
        load(filename, 'Y');
    else
        % Concatenate all the audio data from Svow and reshape the audio
        % data
        allaud = [Svow.aud];
        A = reshape(allaud(:, tps, :), [], size(allaud, 3));

        % Compute the correlation matrix using Spearman correlation coefficient
        B = corr(A, 'Type', 'Spearman');

        % Remove all non-stressed vowel instances
        B(Svow.stress ~= 1, :) = [];
        B(:, Svow.stress ~= 1) = [];
        
        % Apply multidimensional scaling to obtain a 2D representation
        Y = mdscale(abs(B), 2);

        % Save the results in the specified filename for future use
        save(filename, 'B', 'Y');
    end
   
    subplot(2, 1, ctr);
    
    % Create a scatter plot with colored markers for stressed vowels
    scatter(Y(:, 1), Y(:, 2), 35, Svow.vowelType(Svow.stress == 1), 'filled');
    
    % Set the colormap to 'Dark2' from ColorBrewer colormap collection
    colormap(brewermap(5, 'Dark2'));
    
    % TODO: Add a legend with corresponding vowel labels
    legend(unique(Svow.vowel(Svow.stress == 1)), 'Location', 'best');
    ctr = ctr + 1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections;

%% scratch
%% ---------------------- RSA on Neural Data ------------------------------
%% time-course RSA table construction

ueSIDs = {'EC196','EC195','EC183','EC212','EC186','EC219', 'EC221', 'EC222'};
load('out_elecs_voweltypeftest_bychan.mat');
vowrsa = struct();

% load in data from pool of subjects (English vs. Spanish monolinguals)
TDvow = addtoDD(TDvow, 'timit', bef, aft, [sSIDs ueSIDs bSIDs]);
varnames = {'SID', 'el', 'lang', 'neural RDM', 'n2s', 'n2ps', 'n2pe'};
vowresp = array2table(zeros(0,7), 'VariableNames', varnames);
    % n2pe/ps --> phonetic english, phonetic spanish

for s = [ueSIDs sSIDs bSIDs]
    SID = s{1};  
    
    % 1 - spanish, 2 - english, 3 - bilingual
    ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
    
    % use neural window from onset to 100 ms after onset & baseline
    idx = TDvow.stress == 1 & ismember(TDvow.vowel, arpaphn.phnnames);
     
    allr = nan(length(allidx.(SID)), length(TDvow.vowel(idx)), ...
        length(TDvow.vowel(idx)));
    ctr = 1;
    for e = allidx.(SID)
        feat = [];
        
        allr(ctr, :, :) = corr(squeeze(TDvow.(SID).resp(e, bef:bef+10, idx))); 
        B = triu(squeeze(allr(ctr, :, :)));
        B = B(:);
        B(B==0) = NaN;

        % spectral RDM
        A = triu(TDvow.spectralRDM.B(idx, idx));
        A = A(:);
        A(A==0) = NaN;

        vowrsa.(SID).spectral(ctr) = corr(A, B, 'rows','complete');

        A = triu(TDvow.phoneticRDM{1}(idx, idx));
        A = A(:);
        A(A==0) = NaN;

        vowrsa.(SID).phonetic_english(ctr) = corr(A, B, 'rows','complete');

        A = triu(TDvow.phoneticRDM{2}(idx, idx));
        A = A(:);
        A(A==0) = NaN;

        vowrsa.(SID).phonetic_spanish(ctr) = corr(A, B, 'rows','complete');
                
        % add to line to table
        t2 = table({SID}, e, ls, {allr(ctr, :, :)}, vowrsa.(SID).spectral(ctr), ...
            vowrsa.(SID).phonetic_spanish(ctr), vowrsa.(SID).phonetic_english(ctr), ...
            'VariableNames', varnames);
        vowresp = [vowresp; t2];
        
         % turn into average phoneme pair representation
        ctr = ctr+1;
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn;

%% Subject pool information --> replace with idx used above

cols = [0.6000 0.4392 0.6706; 0.3529 0.6824 0.3804; 0.5020 0.4510 0.6745];
% MNI brain with speech responsive electrodes
elecs = containers.Map;
elecs('EC100') = 71;

% initialize design electrode structure
desel=struct();
desel.conds = 1;
desel.sz = [55];%[0 0]; %
desel.cols = [0.6000    0.4392    0.6706];
desel.labels = split(num2str(1:length(desel.conds)));

for s=[sSIDs eSIDs bSIDs]
    %find rows corresponding to subject
    SID = s{1};
    
    desel.(SID).elid = [allidx.(SID)];
    % set up conditions
    conds = ones(1, length(allidx.(SID)));       
    desel.(SID).condition=conds;
end

for key = elecs.keys()
    % selected electrodes
    key = key{1};
    desel.(key).selid = elecs(key);
end

plotMNIElec(sSIDs, desel, 'lh');
plotMNIElec(sSIDs, desel, 'rh');

desel.cols = [0.3529    0.6824    0.3804];
plotMNIElec(eSIDs(1:8), desel, 'lh');
% plotMNIElec(eSIDs(1:8), desel, 'rh');

desel.cols = [0.8784 0.5098 0.0784];
plotMNIElec(bSIDs, desel, 'lh');
plotMNIElec(bSIDs, desel, 'rh');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections arpaphn;

%% visualize with scatter and box plot
figure;
ueSIDs = {'EC196','EC195','EC183','EC212','EC186','EC219', 'EC221', 'EC222'};
h=subplot(2, 2, 1);
els=scatter(vowresp.n2ps(vowresp.lang==1), vowresp.n2pe(vowresp.lang==1), 45, ...
    vowresp.n2s(vowresp.lang==1), 'filled', 'MarkerEdgeColor', 'k'); hold on;
disp(length(els.YData));
ylabel({'English Phonetic','Feature Correlation'});
xlabel({'Spanish Phonetic','Feature Correlation'});
ax = ancestor(h, 'axes');
xticks(0:0.01:0.02); yticks(0:0.01:0.02);
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xtickformat('%.2f');
title(['Spanish subjects (n = ' num2str(length(sSIDs)) ')']);

ylim([-0.005 0.02]); xlim([-0.005 0.02]);
h = refline(1, 0); hold on;
h.LineWidth = 1.5;
h.Color = 'k';
colormap([1 1 1; 1 1 1; brewermap(5, 'OrRd')]);
caxis([0 0.07]);

h=subplot(2, 2, 2);
els = scatter(vowresp.n2ps(vowresp.lang==2), vowresp.n2pe(vowresp.lang==2), 45, ...
    vowresp.n2s(vowresp.lang==2), 'filled', 'MarkerEdgeColor', 'k'); hold on;
disp(length(els.YData));
ylabel({'English Phonetic','Feature Correlation'});
xlabel({'Spanish Phonetic','Feature Correlation'});
ax = ancestor(h, 'axes');
xticks(0:0.01:0.02); yticks(0:0.01:0.02);
ax.XAxis.Exponent = 0;
ax.YAxis.Exponent = 0;
xtickformat('%.2f');
title(['English subjects (n = ' num2str(length(ueSIDs)) ')']);
cbh = colorbar();
ylabel(cbh, 'Spectral Correlation');

ylim([-0.005 0.02]); xlim([-0.005 0.02]);
h = refline(1, 0); hold on;
h.LineWidth = 1.5;
h.Color = 'k';
colormap([1 1 1;  brewermap(6, 'OrRd')]);
caxis([0 0.07]);

ctr = 1;
for s = {sSIDs, ueSIDs}
    
    SIDs= s{1};
    allsid_phn = {[], [], []};
    for s = SIDs
        SID = s{1};
        allsid_phn{1} = [allsid_phn{1} vowrsa.(SID).phonetic_english];
        allsid_phn{2} = [allsid_phn{2} vowrsa.(SID).phonetic_spanish];
        allsid_phn{3} = [allsid_phn{3} vowrsa.(SID).spectral];
    end
    subplot(2, 2, 2+ctr);
    boxchart([allsid_phn{1}; allsid_phn{2}; allsid_phn{3}]', ...
        'MarkerStyle', 'none', 'BoxFaceColor', [0.5 0.5 0.5]); 
    ylabel('corr with neural RDM');
    xticklabels({'English', ...
        'Spanish', 'Spectral'});   
    yline(0);
    ylim([-0.02 0.05]);
    yticks(0:0.02:0.1);

    set(gca, 'FontSize', 13);
    ctr = ctr + 1;
end

%% REDO 
SID = 'EC172';
ne = nan(length(allidx.(SID)), ...
    length(arpaphn.phnnames), length(arpaphn.phnnames));

ctr = 1; 

for e = allidx.(SID)   
    for i = 1:length(arpaphn.phnnames)
        for j = i+1:length(arpaphn.phnnames)
            iid = strcmp(TDvow.vowel,arpaphn.phnnames(i));
            jid = strcmp(TDvow.vowel,arpaphn.phnnames(j));

            alli = squeeze(TDvow.(SID).resp(e, bef:bef+10, idx&iid));
            alli(:, isnan(alli(1, :))) = [];

            allj = squeeze(TDvow.(SID).resp(e, bef:bef+10, idx&jid));
            allj(:, isnan(allj(1, :))) = [];
            ne(ctr, i, j) = mean(corr(alli, allj), [1, 2]);
            ne(ctr, j, i) = mean(corr(alli, allj), [1, 2]);
        end
    end
    ctr = ctr + 1;
end

figure;
imagesc(squeeze(median(ne)));
yticks(1:13); xticks(1:13);
yticklabels(arpaphn.phnnames);
xticklabels(arpaphn.phnnames);
colormap(flipud(brewermap(30, 'Spectral')));
title({'Vowel Distances'}, {'by Neural features'})
set(gca, 'FontSize', 15);


%% plot different formants


f=figure;
subplot(3, 1, 1)
plotVariance(Dvow.formantVals, Dvow.vowel,crest(5), 0, 1, f);
[h, xedges, yedges] = histcounts2(Dvow.formantVals(1, :), ...
    Dvow.formantVals(2, :));
contour(yedges(2:end), xedges(2:end), h, [20 20], ...
    'Color', 'k', 'LineWidth', 2)
legend('off');

subplot(3, 1, 2)
excl = contains(unique(TDvow.vowel), 'x') | ...
    contains(unique(TDvow.vowel), 'm') | ...
    cellfun(@(x) length(x)~=2, unique(TDvow.vowel));
timit_vow2 = unique(TDvow.vowel);
timit_vow2(excl)=[];

plotVariance(TDvow.formantVals(:, ismember(TDvow.vowel,timit_vow2)), ...
    TDvow.vowel(ismember(TDvow.vowel,timit_vow2)), crest(20), 0, 1, f);
[h, xedges, yedges] = histcounts2(TDvow.formantVals(1, :), ...
    TDvow.formantVals(2, :));
contour(yedges(2:end), xedges(2:end), h, [15 15], ...
    'Color', 'k', 'LineWidth', 2);
legend('off');

subplot(3, 1, 3)


bar([length(unique(Dvow.vowel)), sum(~excl)]);
ylabel('Number of categories');



%% NEUR: Construct neural confusion matrices v1 (single trial) 

% setup 
SIDs = sSIDs;

% was set at 35
max_comp = 200;

idxs = {find(Dvow.stress==1); find(Tvow.stress==1 ...
    & ismember(Tvow.vowel, timit_vow))};
test_idx = cell(2, 1);
Svows = {Dvow, TDvow};

% load in all vowel discriminating electrode responses
% create matrix -- el x time (onset:onset+200) x trial
% same electrodes for both corpus representations
load('select_elec/out_elecs_voweltypeftest_bychan.mat');

conf_neur = cell(2, 1);
for i = 1:2 
    % single time frame LDA dimensionality reduction
    Svow = Svows{i};
    idx = idxs{i};       
    
    % load in data from pool of subjects (English vs. Spanish monolinguals)
    %TDvow = addtoDD(TDvow, 'timit', bef, aft, SIDs);
    %Dvow = addtoDD(Dvow, 'dimex', bef, aft, SIDs);
    varnames = {'SID', 'el', 'lang', 'neuralresp', 'fval'};
    vowresp = array2table(zeros(0,5), 'VariableNames', varnames);
    
    % aggregate neural data
    for s = SIDs
        SID = s{1};  
    
        % 1 - spanish, 2 - english, 3 - bilingual
        ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
        
        % use neural window from onset to 300 ms after onset & baseline    
        tps = bef+10:bef+30;
        wind = 21;
        allr = nan(length(allidx.(SID)), wind, ... % length(tps)
            length(Svow.vowel(idx)));
        ctr = 1;
        for e = find(fvals.(SID)>10)
            [~, tp] = max(mean(Svow.(SID).resp(e, tps, idx), 3, 'omitnan'));
            %allr(ctr, :, :) = squeeze(Svow.(SID).resp(e, tps(tp)-2:tps(tp)+2, idx));
            allr(ctr, :, :) = squeeze(Svow.(SID).resp(e, tps, idx)); 
            
            % add to line to table
            t2 = table({SID}, e, ls, {allr(ctr, :, :)}, fvals.(SID)(e), ...
                'VariableNames', varnames);
            vowresp = [vowresp; t2];
            
             % turn into average phoneme pair representation
            ctr = ctr+1;        
        end
        clear allr e ctr t2  
    end
    
    % run pca
    resp = cell2mat([vowresp.('neuralresp')]);
    weights = reshape(repmat(vowresp.fval, 1, wind), ...
        height(vowresp)*wind, []);
   
    [A, ~, ~, y] = makeDataMatrix(resp, Svow.vowelType(idx), ...
        vowresp.SID, 1000);
    assert(~isempty(A));
    
    % find all acoustic test trials s.t. trial overlap in subjects is maximized
    overlapTrl = find(~isnan(mean(A)));
    rng(1);
    test_idx{i} = sort(overlapTrl(randperm(length(overlapTrl), ...
         round(length(idxs{i})*0.25))));
    A_test = A(:, test_idx{i});
    
    % with pca
    % weights = max(squeeze(mean(resp, 2)), [] , 2, 'omitnan');
    [~, score_test, ~, ~, exp] = pca(A_test', 'NumComponents', min(max_comp, ...
        size(A_test, 1))); % 'VariableWeights', weights

    % see how each pc correlates with each formant (and choose each pc that correlates)
%     form_corr = nan(3, size(score_test, 2));
%     for f = 1:2
%         form_corr(f, :) = corr(Svow.formantVals(f, test_idx{i})', ...
%             score_test);        
%     end
%     form_corr(3, :) = corr(Svow.meanf0(test_idx{i})', ...
%             score_test);
%     [~, comps] = maxk(max(abs(form_corr(1:2, :))), 10);
    n_comp = find(diff(cumsum(exp)>95));
    
    % dissimilarity matrix
    %conf_neur{i} = squareform(pdist(score_test(:, 1:n_comp), 'mahalanobis'));
    %conf_neur{i} = squareform(pdist(score_test(:, comps), 'mahalanobis'));
    [means, sem]=grpstats(squeeze(vowresp.neuralresp{1}(:, :, test_idx{i}))', ...
        Svow.vowelType(idxs{i}(test_idx{i})));
    
    conf_neur{i} = squareform(pdist(A_test', 'spearman'));
    figure; 
    corpus_idx = idxs{i}(test_idx{i});
    subplot(1, 2, 1);
    errorbar(means', sem'); legend(unique(Svow.vowel(corpus_idx)));
    subplot(1, 2, 2);
    [~, idx] = sort(Svow.vowelType(corpus_idx)); 
    imagesc(conf_neur{i}(idx, idx));

    clear Svow cb tmp
end

% neural space check
% for i = 1:2
%     figure;
%     y = mdscale(conf_neur{i}, 2);
%     scatter(y(:, 1), y(:, 2), 35, ...
%         Svows{i}.vowelType(idxs{i}(test_idx{i})), 'filled', ...
%         'MarkerFaceAlpha', 0.5); hold on;
%     if i==1
%         cols = brewermap(5, 'Dark2');
%         cols(5, :) = [0.2 0.6 0.9];        
%     else
%         cols=(brewermap(length(unique(Svows{i}.vowelType(test_idx{i}))), ...
%             'Spectral'));
%     end
%     yticks([]); xticks([])
%     xlabel('MD1'); ylabel('MD2');
%     set(gca, 'FontSize', 15);
%     colormap(cols);
% end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr idxs conf_neur test_idx max_comp;


%% NEUR: Construct neural confusion matrices v2 (cross validated rep trials) 

% setup
initSIDs = eSIDs;
% SIDs will be subselected based on largest trial overlap

% was set at 35
max_comp = 200;

test_idx = cell(2, 1);
Svows = {Dvow, TDvow};

% idxs = {find(Dvow.stress==1); ...
%     find(Tvow.stress==1 & ismember(Tvow.vowel, timit_vow))};
idxs = {find(Dvow.stress>-1); ...
    find(ismember(Tvow.vowel, timit_vow))};

% find overlapping repeated indices as test indices 
SIDs = cell(2, 1);
for i = 1:2
    Svow = Svows{i};
    reps = nan(length(initSIDs), length(Svow.vowel));
    for subj = 1:length(initSIDs)
        SID = initSIDs{subj};
        reps(subj, :) = Svow.(SID).repVows;
    end
    
    % find set with most overlap
    [idxs{i}, SIDs{i}] = findOverlap(reps, initSIDs, [], 2);
    if strcmp(Svow.corpus, 'timit')
        idxs{i}(~ismember(Tvow.vowel, timit_vow)) = 0;
    end
    disp(['Total trials: ' num2str(sum(idxs{i}))]);
end

% load in all vowel discriminating electrode responses
% create matrix -- el x time (onset:onset+200) x trial
% same electrodes for both corpus representations
load('out_elecs_voweltypeftest_bychan.mat');
conf_neur = cell(2, 1);
for i = 1:1
    % single time frame LDA dimensionality reduction
    Svow = Svows{i};
    idx = idxs{i};  
    
    % load in data from pool of subjects (English vs. Spanish monolinguals)
    %TDvow = addtoDD(TDvow, 'timit', bef, aft, SIDs);
    %Dvow = addtoDD(Dvow, 'dimex', bef, aft, SIDs);
    varnames = {'SID', 'el', 'lang', 'neuralresp', 'fval'};
    vowresp = array2table(zeros(0,5), 'VariableNames', varnames);
    
    % aggregate neural data
    for s = SIDs{i}
        SID = s{1};  
    
        % 1 - spanish, 2 - english, 3 - bilingual
        ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
        
        % use neural window from onset to 300 ms after onset & baseline    
        tps = bef+5:bef+50;
        wind = 46;
        reps = 10;
        allr = nan(length(allidx.(SID)), wind, ... % length(tps)
            length(Svow.vowel(idx)), reps);
        ctr = 1;
        for e = allidx.(SID)
            [~, tp] = max(mean(Svow.(SID).resp(e, tps, idx, :), [3 4], 'omitnan'));
            allr(ctr, :, :, :) = squeeze(Svow.(SID).resp(e, tps, idx, 1:reps)); 
            
            % add to line to table
            t2 = table({SID}, e, ls, {allr(ctr, :, :, :)}, fvals.(SID)(e), ...
                'VariableNames', varnames);
            vowresp = [vowresp; t2];
            
             % turn into average phoneme pair representation
            ctr = ctr+1;        
        end
        clear allr e ctr t2  
    end
    
    % run pca
    resp = cell2mat([vowresp.('neuralresp')]);
    weights = vowresp.fval;
%     weights = reshape(repmat(vowresp.fval, 1, wind), ...
%         height(vowresp)*wind, []);

    % A is electrode/time x trials x repetitions
    A = squeeze(mean(resp, 2));
    %A = reshape(resp, height(vowresp)*wind, [], reps);
    
    % remove all rows missing

    % assert no columns or rows where data is missing
    assert(all(~(sum(isnan(A))>0), "all"));
    assert(all(~(sum(isnan(A), 2)>0), "all"));
    assert(~isempty(A));

    % find all acoustic test trials     
    test_idx{i} = idxs{i};

    % leave-one-out cross-validated RSA
    tmp = nan(reps, sum(test_idx{i}), sum(test_idx{i}));
    for rep = 1:reps
        % with pca
        % weights = max(squeeze(mean(resp, 2)), [] , 2, 'omitnan');
        a = 1:10;
        A_train = squeeze(mean(A(:, :, (a~=rep)), 3));
        A_test = squeeze(A(:, :, rep));
        [~, score_test, ~, ~, exp] = pca(A_train', 'NumComponents', min(max_comp, ...
            size(A_train, 1)), 'VariableWeights',weights);
        n_comp = find(diff(cumsum(exp)>95));
        
        %% dissimilarity matrix
        %conf_neur{i} = squareform(pdist(score_test(:, 1:n_comp), 'mahalanobis'));
        %conf_neur{i} = squareform(pdist(score_test(:, comps), 'mahalanobis'));
%         [means, sem]=grpstats(squeeze(vowresp.neuralresp{1}(:, :, test_idx{i}, ...
%             (a~=rep)))', Svow.vowelType(idxs{i}(test_idx{i})));
        
        tmp(rep, :, :) = squareform(pdist(score_test(:, 1:n_comp), 'seuclidean'));        
        corpus_idx = idxs{i}; % test_idx{i}
        [~, idx] = sort(Svow.vowelType(corpus_idx)); 
        if rep == 1, figure; imagesc(squeeze(tmp(rep, idx, idx))); end

%         subplot(1, 2, 1);
%         errorbar(means', sem'); legend(unique(Svow.vowel(corpus_idx)));
%         subplot(1, 2, 2);        
    end
    conf_neur{i} = squeeze(mean(tmp));
    clear Svow cb tmp
end

% neural space check
for i = 1:1
    figure;
    y = mdscale(conf_neur{i}, 2);
    scatter(y(:, 1), y(:, 2), 35, ...
        Svows{i}.vowelType(idxs{i}(test_idx{i})), 'filled', ...
        'MarkerFaceAlpha', 0.5); hold on;
    if i==1
        cols = brewermap(5, 'Dark2');
        cols(5, :) = [0.2 0.6 0.9];        
    else
        cols=(brewermap(length(unique(Svows{i}.vowelType(test_idx{i}))), ...
            'Spectral'));
    end
    yticks([]); xticks([])
    xlabel('MD1'); ylabel('MD2');
    set(gca, 'FontSize', 15);
    colormap(cols);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr idxs conf_neur test_idx max_comp;
%% ACS: Construct acoustic confusion matrices from LDA classification

% setup
neur_corr = {[], []; [], []};
acs_corr = cell(2, 1);
idxs = {1:length(Dvow.vowel); find(ismember(Tvow.vowel, timit_vow))};
% idxs = {find(Dvow.stress==1); ...
%     find(Tvow.stress==1 & ismember(Tvow.vowel, timit_vow))};

allaud = {Dvow.aud(:, :, idxs{1}); Tvow.aud(:, :, idxs{2})};
y =  {Dvow.vowelType(idxs{1})'; Tvow.vowelType(idxs{2})'};

reps = 50;

train_idx = cell(2, 1);
cent_dist = nan(2, 2, reps);
labels = {'Spanish', 'English'};

Svows = {Dvow, Tvow};

% repetitions to get distribution across different cross-validated sets 
for rep = 1:reps 
    disp(['--------------------- rep ' num2str(rep) ' ---------------------']);
    A_train = cell(2, 1);
    A_test = cell(2, 1);
    for i = 1:2 % set training and test data
        % use random 80% of the data to train the LDA classifier
        % rng(1);
        % train_idx{i} = randperm(size(allaud{i}, 3), round(size(allaud{i}, 3)*0.75));
        % test_idx{i} = find(~ismember(1:size(allaud{i}, 3), train_idx{i}));

        % use 80% of the data that is not the shared 20% of neural data to train the LDA classifier
        test_idx{i} = find(ismember(idxs{i}, find(Svows{i}.('EC195').repVows)));
        train_idx{i} = find(~ismember(1:size(allaud{i}, 3), test_idx{i}));
    
        % preprocessing and pca   
        A_train{i} = reshape(allaud{i}(:, 50:55, train_idx{i}), [], ...
            size(train_idx{i}, 2))';
        A_test{i} = reshape(allaud{i}(:, 50:55, test_idx{i}), [], ...
            size(test_idx{i}, 2))';
    end
    
    conf_acs = cell(2, 2);
    mapx_acs = cell(2, 2);
    for i = 1:2

        % with pca
        [coeff, score, ~, ~, exp] = pca(A_train{i}, 'Algorithm', 'eig', ...
            'NumComponents', max_comp);
        n_comp = find(diff(cumsum(exp)>95));
        
        % run lda with kfold cross-validation
        kfolds = 5;
        [y_pred, mappedX, Mdl] = LDAmap(score(:, 1:n_comp), y{i}(train_idx{i}), kfolds);

        %     mappedX = squeeze(mean(mappedX, 'omitnan'));        
        %     visualizing LD transofmration for training set
        %     visLDA(mappedX, idx{i}(train_idx), Svow{i}, y_pred, y{i}(train_idx));
        %     title('Training Data');
    
        for j = 1:2
            % applying LD basis transformation for test set (same/cross language)
            mu = mean(A_test{j});
            score_test = (A_test{j}-mu)*coeff(:,1:n_comp);
            kfolds = length(Mdl.Trained);
            [y_pred, mappedtestX] = testmapLDA(Mdl,score_test(:, 1:n_comp), kfolds);
            mappedtestX = squeeze(mean(mappedtestX, 'omitnan'));
            
            % visualize with correct labels or predicted labels
            true_y = y{j}(test_idx{j});
            if rep==1           
                if j ~= i   % true labels don't exist for cross-comparison case
                    true_y = y_pred';
                end
                visLDA(mappedtestX, idxs{i}(test_idx{i}), Svows{i}, y_pred, true_y, ...
                        j==i);
            end
            % train on i, test on j
            conf_acs{i, j} = squareform(pdist(mappedtestX, 'cosine'));
            mapx_acs{i, j} = mappedtestX;
            cent_dist(i, j, rep) = mean(silhouette(mappedtestX, true_y));
        end
    end

    % calculate the clustering index per repetition
    [sort_idx]=findOptSet(conf_acs, mapx_acs, 2, {y{1}(test_idx{1}), ...
        y{2}(test_idx{2})}, {Svows{1}.formantVals(1:3, idxs{1}(test_idx{1})), ...
        Svows{2}.formantVals(1:3, idxs{2}(test_idx{2}))});
    %[sort_idx]=findOptSet(conf_acs, mapx_acs, 1, []);

    % visualize which tokens used
    for i = 1:2
        figure;
        scatter(Svows{i}.formantVals(2, test_idx{i}), ...
            Svows{i}.formantVals(1, test_idx{i}), 35, ...
            Svows{i}.vowelType(test_idx{i}), 'filled', 'MarkerFaceAlpha', 0.2); hold on;
        colormap(getColors(i+1, length(timit_vow)));
        scatter(Svows{i}.formantVals(2, test_idx{i}(sort_idx{i})), ...
            Svows{i}.formantVals(1, test_idx{i}(sort_idx{i})), 35, 'black');
        ylim([200 800]);
        xlim([500 3000]);
        set(gca, 'XDir', 'reverse', 'YDir', 'reverse');

        % pitch distribution of selected vs all tokens

        % category distribution of selected vs all tokens
    end
   
    % correlate
    for i = 1:2
        neur = conf_neur{i};
        acs = conf_acs(:, i);   
        
        for j = 1:2
            % row is listening to 1 - DIMEX, 2 - TIMIT
            % column is trained on 1 - Spanish, 2 - English
            neur_corr{i, j}(rep) = corr(vectorRDM(neur(sort_idx{i}, sort_idx{i})), ...
                vectorRDM(acs{j}(sort_idx{i}, sort_idx{i})), ...
                    'Type', 'Spearman');
        end
    end

    clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr y allaud idxs max_comp test_idx ...
    labels conf_neur acs_corr cent_dist;
end

%% NEUR-ACS: Visualizing correlation results
% listening to 1 - DIMEX, 2 - TIMIT  
ctr = 1;
x = [1 2 4 5];
cols = getColors(1);
g = [];
for i = 1:2
    % trained on 1 - Spanish, 2 - English
    for j = 1:2
        scatter(x(ctr)-0.1+rand(1, length(neur_corr{i, j}))*0.2, neur_corr{i, j}, ...
            45, cols(j, :), 'filled', 'MarkerFaceAlpha', 0.5); hold on;
        ctr = ctr+1;
    end
end
xticks([1.5 4.5]);
yticks(0:0.1:0.4);
yline(0, 'LineWidth', 1.5, 'LineStyle', '--');
xticklabels({'Spanish', 'English'});
xlabel('Listening to...');
legend({'trained on Spanish', 'trained on English'});
ylabel('neural and acoustic correlation');

set(gca, 'FontSize', 15);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr;

%% ------------------ Plain LDA on neural data -----------------------------

Svows = {Dvow, TDvow};
native = [];
nonnative = [];
for vow = 1:2

    % single time frame LDA dimensionality reduction
    Svow = Svows{vow};
    corpus = Svow.corpus;
    switch corpus
        case 'dimex'
            idx = Svow.stress==1; % 
        case 'timit'
            idx = Svow.stress==1 & ismember(Svow.vowel, timit_vow); % 
    end
    
    % load in all vowel discriminating electrode responses
    % create matrix -- el x time (onset:onset+200) x trial
    load(['select_elec/out_elecs_voweltypeftest_bychan_' corpus '_all.mat']);
    
    % load in data from pool of subjects (English vs. Spanish monolinguals)
    %TDvow = addtoDD(TDvow, 'timit', bef, aft, [sSIDs ueSIDs bSIDs]);
    %Dvow = addtoDD(Dvow, 'dimex', bef, aft, [sSIDs ueSIDs]);
    varnames = {'SID', 'el', 'lang', 'neuralresp'};
    vowresp = array2table(zeros(0,4), 'VariableNames', varnames);
        % n2pe/ps --> phonetic english, phonetic spanish
    SIDs = {sSIDs, eSIDs, [], bSIDs};
    disp(['------------------ ' upper(corpus) '------------------------']);
    nfolds = 12;
    AUC = nan(4, nfolds);
    % aggregate neural data
    pos = [1, 2, 0, 3];
    figure;
    for ls = [1:2]
        disp(['-----------------  Language type: ' num2str(ls) ' ------------']);
        % 1 - spanish, 2 - english, 3 - bilingual    
        for s = SIDs{ls} % , 'EC100', 'EC105', 'EC163', 'EC214'
            SID = s{1};     
            % use neural window from onset to 300 ms after onset & baseline    
            tps = bef+10:bef+30;
            allr = nan(length(allidx.(SID)), length(tps), ...
                length(Svow.vowel(idx)));
            ctr = 1;
            for e = allidx.(SID)
                allr(ctr, :, :) = squeeze(Svow.(SID).resp(e, tps, idx)); 
                
                % add to line to table
                t2 = table({SID}, e, ls, {allr(ctr, :, :)}, ...
                    'VariableNames', varnames);
                vowresp = [vowresp; t2];
                
                 % turn into average phoneme pair representation
                ctr = ctr+1;        
            end
            clear allr e ctr t2 
        end
        
        y = Svow.vowelType;
        X = cell2mat([vowresp.('neuralresp')]);

        rng(3)
        randidx = sort(randsample(size(X, 1), 100));
        X = X(randidx, :, :);
        [A, ~, ~, y] = makeDataMatrix(X, y(idx), vowresp.SID(randidx),  1000);
        [~, ~, AUC(ls, :), pcaX, scores, acc, mappedX] = ...
            lda(A', y, 1, [], tps, nfolds);
           
        colors = getColors(1);
        boxchart(ones(nfolds, 1)*pos(ls), AUC(ls, :),  ...
             'BoxFaceColor', colors(ls, :)); hold on;
        
        % conf_neur = pdist2(mappedX, mappedX);
        % [~, sorted] = sort(y);
        % Svow.corpus = corpus;
        % visLDA(mappedX, trls, Svow, y_pred, y);
        
        disp(['accuracy = ' num2str(mean(acc))]);
        disp(['AUC = ' num2str(mean(AUC(ls, :)))]);
        disp(['n trials = ' num2str(length(y))]);
        disp(['n classes = ' num2str(length(unique(y)))])

        if vow==ls
            native = [native; AUC(ls, :)];
        else
            nonnative = [nonnative; AUC(ls, :)];
        end
    end
    
    pos = [1, 2; 2, 3; 1, 3];
    combo = [1, 2; 2, 4; 1, 4];
    for i = 1:3    
        [~, p] = ttest2(AUC(combo(i, 1), :), AUC(combo(i, 2), :));
        if ~isempty(getSigStr(p, 2))
            plot([pos(i, 1)-0.02 pos(i, 2)+0.02], ...
                [0.7+5*i 0.7+5*i], 'LineWidth', 1.5, 'Color', 'k'); hold on;
            text(mean(pos(i, :))-0.25, 0.7+5*i+5, getSigStr(p, 2), 'FontSize', 15);
        else
            disp(p)
        end
    end
    
    %formatting
    xticklabels({'Spanish', 'English' , 'Bilingual'});
    xticks([1 2 3]);
    xlim([0.5 3.5]);
    legend('off');
    yline(100/length(unique(Svow.vowelType(idx))), 'Color', 'k','LineStyle', '--');
    yticks();
    ylim([0.5 0.9]);
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
yticks([0.5 0.9])
set(gca, 'FontSize', 15);
xticklabels({'native', 'non-native'});

[p, ~] = ranksum(native(:), nonnative(:));
plot([1 2], ...
    [0.87 0.87], 'LineWidth', 1.5, 'Color', 'k'); hold on;
text(mean(pos(i, :))-0.25, 0.9, getSigStr(p, 2), 'FontSize', 15);


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr;

%% LDA on neural data for single subjects

Svow = Dvow;
corpus = 'dimex';
Svow.corpus = corpus;
switch corpus
    case 'dimex'
        idx = Svow.stress==1;
    case 'timit'
        idx = Svow.stress==1 & ismember(Svow.vowel, timit_vow);
end

% load in all vowel discriminating electrode responses
% create matrix -- el x time (onset:onset+200) x trial
load('select_elec/out_elecs_voweltypeftest_bychan.mat');

% load in data from pool of subjects (English vs. Spanish monolinguals)
%TDvow = addtoDD(TDvow, 'timit', bef, aft, [sSIDs ueSIDs bSIDs]);
%Dvow = addtoDD(Dvow, 'dimex', bef, aft, [sSIDs ueSIDs]);
varnames = {'SID', 'el', 'lang', 'neuralresp'};
vowresp = array2table(zeros(0,4), 'VariableNames', varnames);
    % n2pe/ps --> phonetic english, phonetic spanish

% aggregate neural data
for s = [sSIDs eSIDs] % , 'EC100', 'EC105', 'EC163', 'EC214'
    SID = s{1};  

    % 1 - spanish, 2 - english, 3 - bilingual
    ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
    
    % use neural window from onset to 300 ms after onset & baseline    
    tps = bef+10:bef+30;
    allr = nan(length(allidx.(SID)), length(tps), ...
        length(Svow.vowel(idx)));
    ctr = 1;
    for e = allidx.(SID)
        allr(ctr, :, :) = squeeze(Svow.(SID).resp(e, tps, idx)); 
        
        % add to line to table
        t2 = table({SID}, e, ls, {allr(ctr, :, :)}, ...
            'VariableNames', varnames);
        vowresp = [vowresp; t2];
        
         % turn into average phoneme pair representation
        ctr = ctr+1;        
    end
    clear allr e ctr t2 
end

ypred_neur = cell(2, 1);
y_neur = cell(2, 1);
nfolds = 5;
AUC = cell(2, 1);
for ls = 1:2
    disp(['---------------------- LS ' num2str(ls) '-------------------------------']);
    ctr = 1;
    for s = unique(vowresp.SID(vowresp.lang==ls))'

        sidx = strcmp(vowresp.SID, s{1}); 
        B = cell2mat([vowresp.('neuralresp')(sidx)]);

        y = Svow.vowelType;
        [A, ~, ~, y] = makeDataMatrix(B, y(idx), vowresp.SID(sidx), 600);
        [fp, tp, AUC{ls}(ctr, :), pcaX, scores, acc, ~] = ...
            lda(A', y, 1, [], tps, nfolds);
        
        ctr = ctr+1;
    end
end
%% visualization
figure;
cols = getColors(1);
for ls = 1:2
    scatter(ls - 0.1 + rand(size(acc{ls}, 2), 1)*0.2, mean(acc{ls}, 1), 40, ...
        cols(ls, :), 'LineWidth',2);
    hold on;
    boxplot(mean(acc{ls}, 1), 'Colors','k', 'Positions', ls, 'widths', 0.45);
end

% chance
yline(1/length(unique(y)), 'LineStyle', '--', 'LineWidth', 2);
xticks(1:2);
xticklabels({'Spanish', 'English'});
set(gca, 'FontSize', 15);
ylabel('Classifier Accuracy (%)');
box off;

% acoustic accuracy
% ceil_acc = acc_acs(1);
% yline(ceil_acc, 'LineStyle', '--', 'LineWidth', 2);

%ylim([1/length(unique(y))-0.05 ceil_acc+0.05]);
xlim([0.65 2.35]);
yticks(0.2:0.2:1);
yticklabels(split(num2str(20:20:100)));

figure;
clear C
for ls = 1:2
    subplot(1, 2, ls);
    y_all = cellfun(@(x) mean(x), y_neur(ls, :), 'UniformOutput', false);
    y_pred = cellfun(@(x) mean(x), ypred_neur(ls, :), 'UniformOutput', false);
    vowels = unique(Svow.vowel);

    for s = 1:length(y_pred)
        if ~isnan(y_all{s})
            C(s, :, :) = confusionmat(y_all{s},round(y_pred{s}));
        end
    end    
    cm = confusionchart(squeeze(round(mean(C))),vowels, ...
        'Normalization', 'row-normalized');
%     cm = confusionchart(squeeze(round(mean(C))),timit_vow, ...
%         'Normalization', 'column-normalized');
    cm.RowSummary = 'row-normalized';
    set(gca, 'FontSize', 15);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr *acc* y;

%% ------------------ Plain LDA on neural data (REPEATS) -----------------------------

% single time frame LDA dimensionality reduction
Svow = TDvow;
corpus = 'timit';
max_comp = 500;
switch corpus
    case 'dimex'
        idx = Svow.stress==1;
    case 'timit'
        idx = Svow.stress==1 & ismember(Svow.vowel, timit_vow);
end

% load in all vowel discriminating electrode responses
% create matrix -- el x time (onset:onset+200) x trial
load('out_elecs_voweltypeftest_bychan.mat');
vowrsa = struct();

% load in data from pool of subjects (English vs. Spanish monolinguals)
%TDvow = addtoDD(TDvow, 'timit', bef, aft, [sSIDs ueSIDs bSIDs]);
%Dvow = addtoDD(Dvow, 'dimex', bef, aft, [sSIDs ueSIDs]);
varnames = {'SID', 'el', 'lang', 'neuralresp'};
vowresp = array2table(zeros(0,4), 'VariableNames', varnames);
    % n2pe/ps --> phonetic english, phonetic spanish

% aggregate neural data
for s = eSIDs % , 'EC100', 'EC105', 'EC163', 'EC214'
    SID = s{1};  

    % 1 - spanish, 2 - english, 3 - bilingual
    ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
    
    % use neural window from onset to 300 ms after onset & baseline    
    tps = bef+10:bef+30;
    reps = 10;
    allr = nan(length(allidx.(SID)), length(tps), ...
        length(Svow.vowel(idx)), reps);
    ctr = 1;
    if ~isempty(Svow.(SID).resp)
        for e = allidx.(SID)
            
            allr(ctr, :, :, :) = squeeze(Svow.(SID).resp(e, tps, idx, :)); 
            
            % add to line to table
            t2 = table({SID}, e, ls, {allr(ctr, :, :, :)}, ...
                'VariableNames', varnames);
            vowresp = [vowresp; t2];
            
             % turn into average phoneme pair representation
            ctr = ctr+1;        
        end
        clear allr e ctr t2 
    end
end

%% run pca
A = reshape(cell2mat([vowresp.('neuralresp')]), ...
    height(vowresp)*length(tps), [], reps);

% remove columns and rows where data is missing
% trials for which over electrodes do not have data
% remanining electrodes with missing trials
% A is electrode/time x trials
nancol = sum(isnan(A), [1 3])>0.95*size(A, 1);
A(:, nancol, :) = [];
nanrow = sum(isnan(A), [2 3])>0;
A(nanrow, :, :) = [];

% find which trials retained
tmp = find(idx);
trls = tmp(~nancol);

y = Svow.vowelType(trls)';

% with pca
[~, score, ~, ~, exp] = pca(mean(A, 3)', 'NumComponents', min(max_comp, ...
        size(A, 1))); 
n_comp = find(diff(cumsum(exp)>95));

% without pca
% coeff = A';
% n_comp = size(A, 1);

% run lda with kfold cross-validation 
kfolds = 5;
[y_pred, mappedX, ~] = LDAmap(score(:, 1:n_comp), y, kfolds);
mappedX = squeeze(mappedX(1, :, :));
conf_neur = pdist2(mappedX, mappedX);
[~, sorted] = sort(y);

Svow.corpus = corpus;
visLDA(mappedX, trls, Svow, y_pred, y);

% calculate silhouette score
s = silhouette(mappedX,y);
ri = RandIndex(y, y_pred);
disp('---------------------- VOWEL NEUR LDA -------------------------------');
disp(['Number of components used = ' num2str(n_comp)]);
disp(['RandIndex = ' num2str(ri)]);
disp(['SI = ' num2str(mean(s))]);
disp(['accuracy = ' num2str(sum(y_pred==y)/length(y))]);
disp(['n trials = ' num2str(length(y))]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* inflections neur_corr;
%% ------------------ Plain LDA on ACS data --------------------------------

% was set at 35
max_comp = 400;
Svows = {Dvow, TDvow};
corpora = {'dimex', 'timit'};
acc_acs = [];
for ctr = 1:2
    Svow = Svows{ctr};
    corpus = corpora{ctr};

    %idx = trls; % from above, matched with neural, Svow.stress==1;
    switch corpus
        case 'timit'
            idx = Svow.stress==1 & ismember(Svow.vowel, timit_vow);
        case 'dimex'
            idx = Svow.stress==1;
    end

    allaud = Svow.aud(:, :, idx);
    y =  Svow.vowelType(idx)';
    
    % preprocessing and pca
    A = reshape(allaud(:, 50:55, :), [], size(allaud, 3));
    
    % with pca
    [coeff, score, ~, ~, exp] = pca(A', 'Algorithm', 'eig', 'NumComponents', max_comp);
    n_comp = find(diff(cumsum(exp)>95));
    
    % no pca
    % coeff = A';
    % n_comp = size(A, 1);
    
    % run lda with kfold cross-validation
    kfolds = 5;
    y = arrayfun(@(x) find(x==unique(y)), y);
    [y_pred, mappedX, Mdl] = LDAmap(score(:, 1:n_comp), y, kfolds);
    mappedX = squeeze(mean(mappedX(:, :, :), 'omitnan'));
    
    Svow.corpus = corpus;
    visLDA(mappedX, idx, Svow, y_pred, y);
    conf_acs = pdist2(mappedX(:, 1), mappedX(:, 2));
    [~, sorted] = sort(y);
    
    % calculate silhouette score
    s = silhouette(mappedX,y);
    ri = RandIndex(y, y_pred);
    disp('---------------------- ACS LDA -------------------------------');
    disp(['RandIndex = ' num2str(ri)]);
    disp(['mSI = ' num2str(mean(s))]);
    disp(['accuracy = ' num2str(sum(y_pred==y)/length(y))]);
    disp(['n trials = ' num2str(length(y))]);

    acc_acs(ctr) = sum(y_pred==y)/length(y);
end

figure;
scatter(1:2, acc_acs.*100, 105, 'k', 'LineWidth', 3); hold on;
plot([0.75 1.25], [0.2 0.2].*100, 'Color', 'k', 'LineWidth', 2);
plot([1.75 2.25], [1/length(timit_vow) 1/length(timit_vow)].*100, ...
    'Color', 'k', 'LineWidth', 2);
ylim([0 100]);
set(gca, 'FontSize', 15);
ylabel('Classifier Accuracy (%)');
yticks(0:50:100);
xticks([1:2]);
xticklabels({'Spanish vowels', 'English vowels'});


%% scratch: cross time (sliding window) frame LDA dimensionality reduction
Svow = Dvow;
corpus = 'dimex';

% load in all vowel discriminating electrode responses
% create matrix -- el x time (onset:onset+200) x trial
ueSIDs = {'EC196','EC195','EC183','EC212','EC186','EC219', 'EC221', 'EC222'};
load('out_elecs_voweltypeftest_bychan.mat');
vowrsa = struct();

% load in data from pool of subjects (English vs. Spanish monolinguals)
%TDvow = addtoDD(TDvow, 'timit', bef, aft, [sSIDs ueSIDs bSIDs]);
%Dvow = addtoDD(Dvow, 'dimex', bef, aft, [sSIDs ueSIDs bSIDs]);
varnames = {'SID', 'el', 'lang', 'neuralresp'};
vowresp = array2table(zeros(0,4), 'VariableNames', varnames);
    % n2pe/ps --> phonetic english, phonetic spanish

for s = [eSIDs] % [eSIDs sSIDs bSIDs]
    SID = s{1};  
    
    % 1 - spanish, 2 - english, 3 - bilingual
    ls=find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, bSIDs}));
    
    % use neural window from onset to 100 ms after onset & baseline
    idx = Svow.stress==1;
    
    tps = bef-10:bef+40;
    allr = nan(length(allidx.(SID)), length(tps), ...
        length(Svow.vowel(idx)));
    ctr = 1;
    for e = allidx.(SID)
        allr(ctr, :, :) = squeeze(Svow.(SID).resp(e, tps, idx)); 
        
        % add to line to table
        t2 = table({SID}, e, ls, {allr(ctr, :, :)}, ...
            'VariableNames', varnames);
        vowresp = [vowresp; t2];
        
         % turn into average phoneme pair representation
        ctr = ctr+1;        
    end
    clear allr e ctr t2 
end

max_comp = 200;

ctr = 1;
wind = 2:2:length(tps)-2;
acc = nan(1, length(wind));
MdlLinear = cell(length(wind), 1);
coeffs = cell(length(wind), 1);
mappedX = [];
for t = wind
    % run pca
    tp_frame = t-1:t+1;
    a = cell2mat(vowresp.('neuralresp'));
    A = reshape((a(:, tp_frame, :)), height(vowresp)*length(tp_frame), []);    

    % remove columns and rows where data is missing
    % trials for which over electrodes do not have data
    % remanining electrodes with missing trials
    nancol = sum(isnan(A))>size(A, 1)/2;
    A(:, nancol) = [];
    nanrow = sum(isnan(A), 2)>0;
    A(nanrow, :) = [];
    % running 

    % run lda
    tmp = find(idx);
    y = Dvow.vowelType(tmp(~nancol));
    kfolds = 5;
    [coeff, ~, ~, ~, exp] = pca(A, 'Algorithm', 'eig', 'NumComponents', max_comp);  
    n_comp = find(diff(cumsum(exp)>95));    
    
    [y_pred, tmpX, MdlLinear{ctr}] = LDAmap(coeff(:, 1:n_comp), y, kfolds); 
    mappedX(ctr, :, :) = squeeze(tmpX(1, :, :));
    coeffs{ctr} = coeff;
    
    acc(ctr)=sum(y_pred==y')/length(y);
    ctr = ctr + 1;
end


%% find maximum timepoint for mapping

figure;
for t = 1:length(wind)
    % use model at maximum time point to remap vowel instances at each time point 
    
    [~, maxtp] = max(acc);
    mdl = MdlLinear{t};
        
    [W, LAMBDA] = eig(mdl.Trained{1}.BetweenSigma, ...
        mdl.Trained{1}.Sigma); 
    lambda = diag(LAMBDA);
    [~, SortOrder] = sort(lambda, 'descend');
    W = W(:, SortOrder);
    mappedX = coeffs{t}(:, 1:size(W, 1))*W(:, 1:2);       
    
    subplot(1, 2, 1);
    scatter(mappedX(:, 1), mappedX(:,  2), 35, y, 'filled', ...
    'MarkerFaceAlpha', 0.6);
    cb = colorbar;

    if strcmp(corpus, 'dimex')
        cb.Ticks = 1.35:0.8:5;
        colormap(brewermap(5, 'Dark2'));
    else
        cb.Ticks = unique(Svow.vowelType(idx))+0.7;
        colormap(brewermap(20, 'Spectral'));
    end

    cb.TickLabels=unique(Svow.vowel(idx));
    ylabel('LD 2'); xlabel('LD 1');
    yticks([]); xticks([]);
    set(gca, 'FontSize', 15);
    
    windtp = (tps(wind)-bef)/100;
    subplot(1, 2, 2);
    plot(windtp, acc, 'LineWidth', 2);
    yline(0.25);
    rectangle('Pos', [windtp(t)-1/100, 0.1, 2/100, 0.7], 'FaceColor', [0.7 0.7 0.7 0.6], ...
        'EdgeColor', 'none');
    ylim([0.1 0.8]);
    
    pause(0.5);
    
end

%% ---------------------- Cluster Evaluation Functions --------------------


% Two sample ttest over time
function [signtp] = tpttest(data, pthresh, debug)

    % data format each cell is a group with dim: time x reps
    assert(size(data{1}, 1)==size(data{2}, 1));
    tps = 1:size(data{1}, 1);    
    p = nan(1, length(tps));
    for t = tps
        [~, p(t)] = ttest2(data{1}(t, :), ...
        data{2}(t, :));
    end
    signtp = p<pthresh;
    
    % remove all single tps
    singles = strfind(signtp, [0 1 0]);
    for s = singles
        signtp(s+1) = 0;
    end
    
    % Plot ERPs for debugging purposes
    if debug
        figure;
        cols = brewermap(2, 'Set1');
        for p = 1:2
            x = -0.2:0.01:0.5;
            y = mean(data{p}, 2, 'omitnan');
            sem = std(squeeze(data{p}), 0, 2, 'omitnan')/...
                sqrt(size(squeeze(data{p}), 2));
            shadedErrorBar(x, y, sem, 'lineProps', {'markerfacecolor', cols(p, :)} ); 
            hold on;
            scatter(x(signtp), repmat(-0.2, sum(signtp), 1), 35, 'filled', 'r');
        end
    end
    signtp = sum(p<pthresh);
end


% rand index
function [ri, ari] = RandIndex(y, y_pred)
    sC = y == y';
    sC(logical(tril(ones(length(y))))) = [];
    
    sK = y_pred == y_pred';
    sK(logical(tril(ones(length(y_pred))))) = [];
    
    % pairs of elements same sets in C and in the same set in K
    a = sum(sK & sC); 
    
    % pairs of elements different sets in C and in differnt set in K
    b = sum(~sK & ~sC); 
    ri = (a + b)/nchoosek(length(y),2);
    
    % adjusted Ri
    % ari = (ri - exp_ri)
end

% calculate accuracy over varying number of components retained from PCA
function [acc, Y_pred] = accOverComp(A, max_comp, y)

    acc = nan(max_comp, 1);
    Y_pred = nan(size(A, 2), max_comp);
    coeff = pca(A, 'Algorithm', 'eig', 'NumComponents', max_comp);    
    for n = 1:max_comp       
        MdlLinear = fitcdiscr(coeff(:, 1:n), y, 'DiscrimType', 'linear');
        % turns into partitioned model
        [Y_pred(:, n),~, ~] = predict(MdlLinear, coeff(:, 1:n));
        acc(n) = sum(Y_pred(:, n)==y)/length(y);
    end
end

% assumes no valid entry is zero
function [vec] = vectorRDM(rdm)
% process RDM so can be used for correlation analysis
    % ensure it is a dissimilarity matrix
    assert(all(diag(rdm)==0));
    % take only below triangular 
    mask = logical(tril(ones(size(rdm, 1), size(rdm, 2)), -1));
    vec = rdm(mask);
end

function [y_pred, mappedX] = testmapLDA(mdl, X, kfolds)
% uses a pre-trained LDA classifier to predict labels for test set
% and map test set to trained LD space

    dims = 3;
    mappedX = nan(kfolds, size(X, 1), dims);
    allW = nan(kfolds, size(X, 2), dims);
    for k = 1:length(mdl.Trained)   
        [W, LAMBDA] = eig(mdl.Trained{k}.BetweenSigma, ...
            mdl.Trained{k}.Sigma); 
        lambda = diag(LAMBDA);
        [~, SortOrder] = sort(lambda, 'descend');
        W = W(:, SortOrder);
        allLDs = X*W;
        mappedX(k, :, :) = allLDs(:, 1:dims);
        % return corresponding weights for LDs
        allW(k, :, :) = W(:, 1:dims);
    end
    % find average mapping between trained models
    y_pred = ensemPredict(mdl, X, kfolds);

end

function [y_pred] = ensemPredict(partitionMdl, X, kfold)
    for k = 1:kfold
        mdl = partitionMdl.Trained{k};
        [tmp, ~, ~] = mdl.predict(X);
        y(k, :) = tmp;
    end
    y_pred = mode(y);
end

function [colors] = getColors(type, numel)

    switch type
        case 1 % language type
            % spanish, english, mandarin
            colors = brewermap(4, 'Dark2');
        case 2 % spanish language vowels
            colors = brewermap(5, 'Dark2');
            colors(5, :) = [0.2 0.6 0.9];
        case 3 % english language vowels
            colors=(brewermap(numel, 'Dark2'));
            colors(9, :) = [0.2 0.8 0.9];
    end
end

function [trls, selSIDs] = findOverlap(pres, SIDs, minTrl, minSubj)
    % pres as subjects x trls (1 indicates trial was presented to subject)

    varnames = {'SID set', 'trials', 'total'};
    overlapresp = array2table(zeros(0,3), 'VariableNames', varnames);
    % brute force, find as many overlaps as possible
    for k = minSubj:length(SIDs)
        C = nchoosek(1:length(SIDs), k);
        for r = 1:size(C, 1)
            % for this combination of subjects find overlap
            ovrlp = sum(pres(C(r, :), :));
            t2 = table({C(r, :)}, {ovrlp==k}, sum(ovrlp==k), ...
                'VariableNames', varnames);
            overlapresp = [overlapresp; t2];
        end
    end

    maxidx = find(overlapresp.total==max([overlapresp.total]), 1, 'last');
    selSIDs = SIDs(overlapresp.("SID set"){maxidx});
    trls = overlapresp.trials{maxidx};
end

function [dists] = calcClustDist(X, category)
    dists = nan(size(X, 1), 1);
    for cat = unique(category)'
        idx = find(category == cat);
        dists(idx) = pdist2(X(idx, :), median(X(idx, :)));
    end
end

function [C] = mat_prepro(B)
% preprocessing a matrix of form el x time x trials
    
    % find window with greatest variability
    tps = 3:size(B, 2)-2;
    for ctr = 1:length(tps)
        wind = tps(ctr)-2:tps(ctr)+2;
        tp_mean(ctr) = mean(B(:, wind, :), [1, 2, 3], 'omitnan');
    end

    [~, idx] = max(tp_mean);
    B = squeeze(mean(B(:, tps(idx)-2:tps(idx)+2, :), 2, 'omitnan'));
    C = rescale(B,'InputMin', min(B,[],1), 'InputMax', max(B,[],1));
end

function [sort_idx]=findOptSet(conf_acs, coord_acs, type, category, formants)

% finds the optimal set of tokens for acoustic dissimilarity analysis 
% (between English and Spanish)
% sorted_idx is the list of token indices

    % sort acoustic RDM by single col correlation to find set that will make up
    % lowest correlated acoustic RDMs
    cs = cell(2, 1);
    sort_idx = cell(2, 1);
    acs_corr = cell(2, 1);
    labels = {'Spanish', 'English'};
    for i = 1:2
        debug = 1;
        switch type
            case 1 % correlation method, finding least correlated columns of matrix                
                [~, sort_idx{i}] = sort(diag(corr(conf_acs{1, i}, conf_acs{2, i}, ...
                    'Type', 'Spearman')), 'ascend');

            case 2 % procrustes on both acoustic spaces
                y1 = coord_acs{1, i};
                y2 = coord_acs{2, i};
                [~, y2_trans] = procrustes(y1, y2);
                
                if debug
                    disp(['Original correlations: ' mat2str(diag(corr(y1, y2)))]);
                    disp(['New correlations: ' mat2str(diag(corr(y1, y2_trans)))]);
%                     figure;
%                     for x = 1:3
%                         subplot(1, 3, x);
%                         scatter(y1(:,x), y2_trans(:,x), 45, y2(:,x), 'filled');
%                     end

                    
                    ys = {y1, y2_trans};
                    %ys = {y1_form, y2_form};
                    cats = cell(2, 1);
                    figure;
                    for x = 1:2
                        ax(x) = subplot(1, 2, x);
                        cats{x} = arrayfun(@(x) find(x==unique(category{i})), category{i});
                        scatter(ys{x}(:,1), ys{x}(:,2), 45, cats{x}, 'filled');
                        colormap(getColors(i+1, length(unique(category{i}))));     
                        set(gca, 'FontSize', 15);
                    end
                    linkaxes(ax);

                    % comparing to formant space
                    [~, y1_form] = procrustes(formants{i}', y1);
                    [~, y2_form] = procrustes(formants{i}', y2);
                    corr(y1_form, formants{i}');
                    corr(y2_form, formants{i}');

                    figure;
                    violin([silhouette(ys{1}, cats{1}) ...
                        silhouette(ys{2}, cats{2})], 'medc', [], ...
                        'facecolor', [0.7 0.7 0.7], 'bw', [0.2 0.2]);
                    ylabel('silhouette score');
                    xticks([1 2]);
                    xticklabels({'cross-language', 'same language'});
%                     xticklabels({'same language', 'cross-language'});
                    set(gca, 'FontSize', 15);

                    

                end
    
%                 [~, sort_idx{i}]=sort(diag(pdist2(y1, ...
%                     y2_trans)), 'ascend');   
                [~, sort_idx{i}]=sort(mean(pdist2(y1, ...
                    y2_trans)), 'ascend');
        end

        % remove first xx rows to see how much correlation drops
        rd1 = conf_acs{1, i}(sort_idx{i}, sort_idx{i});
        rd2 = conf_acs{2, i}(sort_idx{i}, sort_idx{i});           
        for j = 30:size(conf_acs{1, i}, 1)        
            cs{i}(j) = corr(vectorRDM(rd1(1:j, 1:j)), ...
                vectorRDM(rd2(1:j, 1:j)), ...
                'Type', 'Spearman');      
        end        
        sort_idx{i} = sort_idx{i}(1:find(cs{i}<0.3, 1, 'last'));
    
        % calculate correlation of new test set idx
        rd1_ct = rd1(1:length(sort_idx{i}), 1:length(sort_idx{i}));
        rd2_ct = rd2(1:length(sort_idx{i}), 1:length(sort_idx{i}));
        acs_corr{i} = corr(vectorRDM(rd1_ct), vectorRDM(rd2_ct), ...
            'Type', 'Spearman');
    
        % Acoustic correlation of this repetition
        disp([labels{i} ' tokens correlation: ' num2str(acs_corr{i})]); 
        disp(['num trials used: ' num2str(length(sort_idx{i}))]);  
    
        if debug % visualization of analysis step
            figure;     
            subplot(1, 3, 1); imagesc(rd1);
            xline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--'); 
            yline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--');
    
            subplot(1, 3, 2); imagesc(rd2);
            xline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--'); 
            yline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--');
    
            subplot(1, 3, 3); scatter(1:length(cs{i}), cs{i}, ...
                20, 'k', 'filled', 'MarkerFaceAlpha', 0.7); hold on
            ylabel('cumulative correlation between matrices');
    
            yyaxis right;
            scatter(1:length(cs{i}), max([var(rd1); var(rd2)]), ...
                 30, 'r', 'filled', 'MarkerFaceAlpha', 0.5);
            xlim([20 length(cs{i})]);
            ylabel('variance in column');
            xlabel('tokens used');
            xline(length(sort_idx{i}), 'LineWidth', 2, 'LineStyle', '--');
        end
    end
end
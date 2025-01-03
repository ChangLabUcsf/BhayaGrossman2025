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
[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
SIDs = [sSIDs, eSIDs, bSIDs, mSIDs]; 

% remove EC225 % 'EC225', 'EC212', 'EC221', 
SIDs(ismember(SIDs, {'EC222'})) = [];

% Dvow = loadDD('dimex', bef, aft, {}, datapath);
% TDvow = loadDD('timit', bef, aft, {}, datapath);

% load in all beta model versions

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
asccd_details = load('out_sentence_details_accsd_loudness.mat');
tps = 50:55;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% load in model beta weights

% load in model weights from strf

modelnames={'onset_phnfeatonset_maxDtL', 'onset_aud',  'onset_maxDtL'};  
thresh = 0.05;
wind = 1:30;

corpus = {'timit', 'dimex', 'asccd'};
details = {timit_details, dimex_details, asccd_details};

% make a table with electrode, corpus, model, and weights
varnames = {'SID', 'el', 'ls', 'timit_rsq', 'dimex_rsq', 'asccd_rsq', ...
    'timit_wts', 'dimex_wts', 'asccd_wts'};
fsitable = table('Size', [0, length(varnames)], 'VariableTypes', ...
    {'string', 'double', 'double', 'double', 'double', 'double', 'cell', 'cell', 'cell'}, ...
    'VariableNames', varnames);

for s = SIDs
    SID = s{1};
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs})); 

    uv_phnfeat = cell(length(corpus), 1);
    for c = 1:length(corpus)
        corpusStrf{c} = loadMultModelStrf(SID, modelnames, corpus{c}, ...
            datapath, 1);
        if ~isempty(corpusStrf{c}{1}) && ~isempty(corpusStrf{c}{3})
            uv_phnfeat{c} = corpusStrf{c}{1}.meanTestR.^2 ...
                            - corpusStrf{c}{3}.meanTestR.^2;
        end
    end

    % find idx of corpus that are not empty
    corpusidx = find(~cellfun(@(x) isempty(x{2}), corpusStrf));
    
    % get unique variance for phonetic features
    
    % find electrodes with significant encoding
    elidx = cellfun(@(x) x{2}.meanTestR.^2 > thresh, ...
        corpusStrf(corpusidx), 'UniformOutput',false);

    % find electrodes with significant encoding in any corpus
    % only use minimum electrodes across all corpora
    minumel = min(cellfun(@(x) length(x), elidx));
    elidx = cellfun(@(x) x(1:minumel), elidx, 'UniformOutput', false);
    elidx = find(any(cell2mat(elidx') & uv_phnfeat{1}(1:minumel)>0.0));
    rsq = cell(2,3);
    wts = cell(2,3);
    
    for c = 1:length(corpus)
        for m = 1:length(modelnames)
            if ~isempty(corpusStrf{c}{m})
                rsq{m,c} = corpusStrf{c}{m}.meanTestR(elidx).^2;
                wts{m,c} = squeeze(num2cell(corpusStrf{c}{m}.meanStrf(:,wind, ...
                    elidx), [1, 2]));
            else % fill with nans
                rsq{m,c} = nan(1, length(elidx));
                wts{m,c} = cell(length(elidx), 1);
            end
        end
    end

    if length(elidx)>1
        % add to table
        fsitable = [fsitable; table(repmat(SID, length(elidx), 1), ...
            elidx', repmat(ls, length(elidx), 1), cat(1, rsq{:,1})', ...
            cat(1, rsq{:,2})', cat(1, rsq{:,3})', ...
            cat(2, wts{:,1}), cat(2, wts{:,2}), cat(2, wts{:,3}), ...
            'VariableNames', varnames)];
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* fsitable;

%% make a hierarchical clustering of model weights for phonetic features

corpus = {'timit', 'dimex', 'asccd'};
wind = 1:30;
model = 1; % phonetic features
pcaflag = 0;

% aggregate all weights into a single matrix
% each row is a model, each column is an electrode
allwts = cell(length(corpus), 1);
for c = 1:length(corpus)
    if model==1
        % electrode x feature x time
        allwts{c} = nan(height(fsitable), 14, length(wind));
    elseif model==2
        % electrode x frequence x time
        allwts{c} = nan(height(fsitable), 80, length(wind));
    end
end

for el = 1:height(fsitable)
    for c = 1:length(corpus)
        if ~isempty(fsitable{el, [corpus{c} '_wts']}{model})
            if model==1 && ismember(c, [1, 2])
                % only use the weights 2-15
                allwts{c}(el, :, :) = fsitable{el, [corpus{c} '_wts']}{model}(2:15, :);
            elseif model==1 && ismember(c, 3)
                allwts{c}(el, 1:12, :) = fsitable{el, [corpus{c} '_wts']}{model}(2:13, :);
            elseif model==2
                allwts{c}(el, :, :) = fsitable{el, [corpus{c} '_wts']}{model}(2:end, :);
            end
        end
    end
end

% do a comparison between hierarchical clustering
clusters = [5, 5];
corp_details = {timit_details, dimex_details, asccd_details};
%figure; 
ctr = 1;
corps = [1, 2];

% find elidx for which both clusters have weights
elidx = find(~isnan(sum(allwts{corps(1)}, [2, 3])) ...
        & ~isnan(sum(allwts{corps(2)}, [2, 3])) ...
        & fsitable.timit_rsq(:, 2)>0.1 ...
        & fsitable.dimex_rsq(:, 2)>0.1 );

for c = corps
    % cluster the weights
    details = corp_details{c};

    % find elidx for which both clusters have weights
    elidx = find(~isnan(sum(allwts{corps(1)}, [2, 3])) ...
            & ~isnan(sum(allwts{corps(2)}, [2, 3])) ...
            & fsitable.([corpus{c} '_rsq'])(:, 2)>0.1);

    if model==1
        % remove syllabic, obstruent, sonorant
        syllidx = find(contains(details.features.names, ...
            {'syllabic', 'obstruent', 'sonorant'}));    
        
        % use absolute value of weights
        wts = abs(allwts{c});
        wts(:, syllidx, :) = [];
    else
        wts = allwts{c};
    end

    % vectorize weights so electrode x (feature x time)
    wts = wts(elidx, :, :);
    vecwts = reshape(wts, size(wts, 1), []);

    % normalize by electrode (?)
%     vecwts = normalize(vecwts, 1, "zscore");
    vecwts = rescale(vecwts);

%     % normalize weights and threshold
%     fthresh = prctile(vecwts, 0);
%     vecwts(vecwts<fthresh) = 0;
   
    
    if pcaflag
        % run pca to reduce dimensionality
        [coeff, score, latent, tsquared, explained, mu] = pca(vecwts);
        numcomp = min(45, size(coeff, 2));
        vecwts = score(:, 1:numcomp);
        disp(sum(explained(1:numcomp)));
        pcaflag=1;
    else
        numcomp = 45;
        mu = 0.1;
        coeff = nan(1, numcomp);
    end
    
    rng(3);
    % cluster the weights using hierarchical clustering
    link = linkage(vecwts, 'ward');
    idx = cluster(link, "maxclust", clusters(ctr));
    %idx = kmeans(wts, clusters);

    [idx, meanwts] = orderCluster(vecwts, wts, idx, clusters(ctr), pcaflag, ...
        coeff(:,1:numcomp), mu);
    disp(crosstab(idx))
    
    % plot mean weights per cluster
    debug = 0;
    if debug
        figure;
        for clust = 1:clusters(ctr)
            subplot(1, clusters(ctr), clust)
            imagesc(meanwts{clust});
            set(gca, 'Ydir', 'normal');

            if model==1
                yticks(1:14)
                features = details.features.names;
                features(syllidx) = [];
                yticklabels(features);
%                 caxis([0 3]);
            elseif model==2
                yticks([1 80]);
                yticklabels({'0', '8'});
                ylabel('kHz');
                %set(gca, 'Yscale', 'log')
%                 caxis([-3, 3]);
                colormap(turbo)
            end
        end
    end

    [~, nidx] = sort(idx);

    subplot(2, 2, ctr)
    dendrogram(link);

    wts = squeeze(mean(wts, 3));

    subplot(2, 2, ctr+2)
    imagesc(wts(nidx, :)');
% 
%     % don't show the 0s, find non nan electrodes
    xlim([0 sum(~isnan(wts(nidx, 1)))]);

    % threshold the color axis
    colormap(flipud(pink));
%     caxis([0, 4]);

    % show where clusters ends are
    hold on;
    clustidx = find(diff(nidx(nidx)));
    xline(clustidx+.5, 'k--', 'LineWidth', 2);

    % yticks
    yticks(1:14)
    features = details.features.names;
    features(syllidx) = [];
    yticklabels(features)
    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* fsitable;

%% make a hierarchical clustering of model weights for acoustic features

corpus = {'timit', 'dimex', 'asccd'};
% aggregate all weights into a single matrix
% each row is a model, each column is an electrode
allwts = cell(length(corpus), 1);
for c = 1:length(corpus)
    allwts{c} = nan(height(fsitable), 14);
end

model = 2; % acoustic features
for el = 1:height(fsitable)
    for c = 1:length(corpus)
        if ~isempty(fsitable{el, [corpus{c} '_wts']}{model}) && ismember(c, [1, 2])
            % only use the weights 2-15
            allwts{c}(el, :) = fsitable{el, [corpus{c} '_wts']}{model};
        elseif ~isempty(fsitable{el, [corpus{c} '_wts']}{model}) && ismember(c, 3)
            allwts{c}(el, 1:12) = fsitable{el, [corpus{c} '_wts']}{model}(2:13);
        end
    end
end

%% functions

function [idx] = argmax(arr, dim)
    [~, idx] = max(arr, [], dim);
end

function [idx, meanwts] = orderCluster(vecwts, wts, idx, numclust, pcaflag, coeff, mu)
    
    % get mean of clusters and order by mean weights
    % perform inverse pca and reshape to be electrode x feature x time
    meanwts = cell(numclust, 1);
    for clust = 1:numclust
        if pcaflag
            % perform inverse pca
            tmp = coeff*mean(vecwts(idx==clust, :), 1)' + mu';
    
            % reshape to be electrode x feature x time
            meanwts{clust} = reshape(tmp, size(wts, 2), size(wts, 3));
        else
            meanwts{clust} = reshape(mean(vecwts(idx==clust, :), 1), ...
                size(wts, 2), size(wts, 3));
        end
    end

    % relabel cluster idx by weight order
    [~, sidx] = sort(argmax(squeeze(mean(cat(3, meanwts{:}), 2)), 1));
    nidx = nan(size(idx));
    for clust = 1:length(sidx)
        nidx(idx==sidx(clust)) = clust;
    end
    clear sidx idx
    idx = nidx;

    for clust = 1:numclust
        if pcaflag
            % perform inverse pca
            tmp = coeff*mean(vecwts(idx==clust, :), 1)' + mu';
    
            % reshape to be electrode x feature x time
            meanwts{clust} = reshape(tmp, size(wts, 2), size(wts, 3));
        else
            meanwts{clust} = reshape(mean(vecwts(idx==clust, :), 1), ...
                size(wts, 2), size(wts, 3));
        end
    end
end
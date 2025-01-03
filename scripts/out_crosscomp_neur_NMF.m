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
% This is the list of subjects that have both TIMIT and DIMEx data
SIDs = [sSIDs, eSIDs]; % , {'HS11', 'HS9', 'HS10'}

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
tps = 50:55;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% -------------------- Loading sentence (and stitched) responses  ----- %%

% load sentence responses
[sent_encoding] = loadSentenceResponse(SIDs, timit_details, dimex_details, datapath);

%% -------------------- Raw single sentence responses (mean over repeats) --------------- %%

lang_fields = {'en_rep_resp', 'sp_rep_resp'};
repname_fields = {'en_repname', 'sp_repname'};
titles = {'TIMIT', 'DIMEx'};
lang_exp = {'Spanish mono', 'English mono'};
sp_repname = {'s00104', 's00804', 's01904', 's03004', 's04704', 's05004', 's06104', 's06804', 's07404', 's09004'};
en_repname = {'fcaj0_si1479', 'fcaj0_si1804', 'fdfb0_si1948', 'fdxw0_si2141', 'fisb0_si2209', 'mbbr0_si2315', 'mdlc2_si2244', 'mdls0_si998', 'mjdh0_si1984', 'mjmm0_si625'};
repnames = {en_repname, sp_repname};
min_rep = [10, 8];
X_all = cell(2, 2, 10);

for lang = 1:2
    for ls = 1:2
        idx = sent_encoding.ls==ls;
        tmp_encoding = sent_encoding(idx, :);

        figure;
        sbplt = 1;
        for sent = repnames{lang}

            % elecs x tps
            X = nan(sum(idx), 500);
            
            % first try NMF for single sentence responses
            els = zeros(height(tmp_encoding), 1); % binary for which elecs are included
            ctr = 1;
            for i = 1:height(tmp_encoding)
                
                % find the sentence index from the repname field
                sentid = find(strcmp(tmp_encoding.(repname_fields{lang}){i}, sent));
                rep_resp = tmp_encoding.(lang_fields{lang}){i};
                rep_resp = squeeze(rep_resp(:, :, sentid));

                if ~isempty(rep_resp) && size(rep_resp, 2)>min_rep(lang)
                    % time points x repeats x sentence
                    sent_length = find(isnan(rep_resp(:, 1)), 1, 'first')-1;
                    x = rep_resp(1:sent_length, 1:min_rep(lang));

                    % use all repeats
                    %X(ctr, 1:sent_length) = squeeze(mean(x, 2));

                    % use only the first repeat
                    X(ctr, 1:sent_length) = squeeze(x(:, 1));

                    els(i) = 1;
                    ctr=ctr+1;
                end
            end

            % remove all nan rows
            nanrow = isnan(X(:, 1));
            X(nanrow, :) = [];
            X = X(:, 1:sent_length);
            
            % remove all nan columns (quick fix to error)
            nancol = sum(isnan(X), 1)>0;
            X(:, nancol) = [];

            if ~isempty(X)
                subplot(2, 5, sbplt);
                imagesc(normX(X));
                colormap(flipud(rdylgn));
                title(sent)
            end
            
            % save out X matrix
            X_all(lang, ls, sbplt) = {X};

            sbplt = sbplt+1;
        end

        sgtitle([titles{lang} ' ' lang_exp{ls}]);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% -------------------- NMF for single sentence responses --------------- %%

% Uses X_all from previous section
nd = 12;
H_all = cell(2, 2);

for lang = 1:2
    for ls = 1:2

        H_all{lang, ls} = nan(10, nd, 500);

        figure;
        for sent = 1:10

            X = X_all{lang, ls, sent};

            if isempty(X)
                continue;
            end

            [~, H, X] = NMFwrapper(X, nd);
            
            % sort X by peak latency
            % figure;
            % [~, maxtp] = max(X, [], 2);
            % [~, tpsort] = sort(maxtp);
            % imagesc(X(tpsort,  :));
            % colormap(flipud(spectral))

            [~, newidx] = sort(argmax(H'));
            H = H(newidx, :);

            subplot(2, 5, sent)
            imagesc(H);
            H_all{lang, ls}(sent, 1:size(H, 1), 1:size(H, 2)) = H;
            
            % figure;
            % subplot(nd, 2, 2:2:nd*2)
            % cx= [];
            % xx = [];
            % bx = [];
            % hemis = cellfun(@(x) imgall.(x).hemi, ...
            %         sent_encoding.SID(find(els)), 'UniformOutput', false);
            % langpval = nan(1, nd);

            % for i = 1:nd
            %     idx_sp = sent_encoding.ls(xidx)==1;
            %     idx_en = sent_encoding.ls(xidx)==2;
            % 
            %     cx = [cx; ones(sum(idx_sp), 1); 2*ones(sum(idx_en), 1)];
            %     xx = [xx; i*2*ones(sum(ismember(sent_encoding.ls(xidx), [1, 2])), 1)];
            %     bx = [bx; W(idx_sp, i); W(idx_en, i)];
            %     legend({'Spanish', 'English'})
            % 
            %     % lme
            %     tbl = table();
            %     tbl.effect = W(:, i);
            %     tbl.language = sent_encoding.ls(find(els));
            %     tbl.hemi = hemis;
            %     tbl.subj = sent_encoding.SID(find(els));
            %     tbl.elec = sent_encoding.el(find(els));
            %     lmeHemi = fitlme(tbl,'effect~1+language+hemi+(1|subj:elec)');
            %     disp(lmeHemi)
            %     disp(['pvals: ' num2str(lmeHemi.Coefficients.pValue(2:3)')]);
            % 
            %     langpval(i) = lmeHemi.Coefficients.pValue(2);
            % end
            % boxchart(xx, bx, 'GroupByColor', cx);
            % legend({'Spanish', 'English'})
            % sgtitle([titles{lang} ': sent ' num2str(sent)]);
            % 
            % % plot the NMF components
            % subplot(nd, 2, 1:2:nd*2);
            % for i = 1:nd 
            %     if langpval(i)<0.01
            %         plot(H(i, :), 'LineWidth', 2); hold on;
            %     else
            %         plot(H(i, :), 'LineWidth', 1);
            %     end
            % end
            % xline([50, sent_length-50], 'HandleVisibility', 'off');
            % legend({'cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'})
        end
    end
end

%%

nd = 12;
% normalize time between 0 and 1, relative time 
titles = {'TIMIT', 'DIMEx'};
lang_exp = {'Spanish mono', 'English mono'};

for lang = 1:2
    figure;
    for ls = 1:2
        samples = 200;
        cmap = hsv;
        cmap = cmap(round(linspace(1, 256, nd)), :);
        H_normall = nan(size(H_all{lang, ls}, 1), size(H_all{lang, ls}, 2), samples);
        for sent = 1:10
            % find the first nan
            [~, tp] = find(isnan(H_all{lang, ls}(sent, 1, :)), 1);
            tp = tp-1;
        
            if tp~=0
                % linearly sample 100 time points from the sentence
                sampidx = round(linspace(50, tp, samples));
                H_normall(sent, :, :) = smoothdata(squeeze(H_all{lang, ls}(sent, :, sampidx))', ...
                    'SmoothingFactor',0.1)';
            end
        end
        
        subplot(2, 1, ls)
        for comp = 1:nd
            plot(squeeze(mean(H_normall(:, comp, :), 1, "omitnan"))', ...
                'LineWidth',2, 'Color',cmap(comp, :)); hold on;
        end
        title([titles{lang} ' ' lang_exp{ls}]);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% NMF for average sentence responses

X = cat(1, sent_encoding.sp_mresp{:});

nd = 4;
[W, H, X_mod] = NMFwrapper(X, nd);
%%

% reorder H and W by peak latency
[~, newidx] = sort(argmax(H'));
H = H(newidx, :);
W = W(:, newidx);

figure;
subplot(nd, 2, 2:2:nd*2)
cx= [];
xx = [];
bx = [];
hemis = cellfun(@(x) imgall.(x).hemi, ...
        sent_encoding.SID, 'UniformOutput', false);
langpval = nan(1, nd);
for i = 1:nd
    spidx = sent_encoding.ls==1;
    enidx = sent_encoding.ls==2;

    cx = [cx; ones(sum(spidx), 1); 2*ones(sum(enidx), 1)];
    xx = [xx; i*2*ones(sum(ismember(sent_encoding.ls, [1, 2])), 1)];
    bx = [bx; W(spidx, i); W(enidx, i)];
    legend({'Spanish', 'English'})

    % lme
    tbl = table();
    tbl.effect = W(:, i);
    tbl.language = sent_encoding.ls;
    tbl.hemi = hemis;
    tbl.subj = sent_encoding.SID;
    tbl.elec = sent_encoding.el;
    lmeHemi = fitlme(tbl,'effect~1+language+(1|subj)+(1|elec:subj)');
    disp(lmeHemi)
    disp(['pvals: ' num2str(lmeHemi.Coefficients.pValue(2)')]);

    langpval(i) = lmeHemi.Coefficients.pValue(2);
end
b = boxchart(xx, bx, 'GroupByColor', cx, 'JitterOutliers','on', ...
    'MarkerStyle','.');
xticks(2:2:8);
xticklabels(split(num2str(1:nd)))
ylabel('cluster weight');
xlabel('cluster')
legend({'Spanish mono', 'English mono'});
set(gca, 'FontSize', 13);

% plot the NMF components
subplot(nd, 2, 1:2:nd*2);
cols = hawaii(nd+3);
cols = cols(2:2+nd, :);
for i = 1:nd 
    if langpval(i)<0.001
        plot(H(i, :), 'LineWidth', 2, 'Color', cols(i, :)); hold on;
    else
        plot(H(i, :), 'LineWidth', 2, 'Color', cols(i, :), ...
            'LineStyle','--');
    end
end
xline([50, 150, 200], 'HandleVisibility', 'off');
legend({'cluster 1', 'cluster 2', 'cluster 3', 'cluster 4'});
set(gca, 'FontSize', 13);


%% NMF for multiple sentence with DTW

% aggregate single sentence responses
imgall = load_allimgdata;
lang_fields = {'en_rep_resp', 'sp_rep_resp'};
titles = {'TIMIT', 'DIMEx'};
sents = [10, 8];
min_rep = [10, 8];

% single sentence
lang = 1;
sent = 8;
% elecs x tps
X = nan(height(sent_encoding), 500);
xidx = 1:height(sent_encoding);

% first try NMF for single sentence responses
els = zeros(height(sent_encoding), 1);
for i = 1:height(sent_encoding)
    rep_resp = sent_encoding.(lang_fields{lang}){i};
    if ~isempty(rep_resp) && size(rep_resp, 2)>min_rep(lang)
        % time points x repeats x sentence
        sent_length = find(isnan(rep_resp(:, 1, sent)), 1, 'first')-1;
        x = rep_resp(1:sent_length, 1:min_rep(lang), sent);
        X(i, 1:sent_length) = squeeze(mean(x, 2));
        els(i) = 1;
    end
end

% perform DTW on each sentence response so that they are aligned
s
% perform NMF on the aggregated responses

% plot the NMF results



%% functions

function [W, H, X] = NMFwrapper(X, nd)

    % X should be of form electrode/sentence x time
    X = normX(X);   

    opt = statset('maxiter',500,'display','final','TolX',1e-6,'TolFun',1e-6);
    [W0,H0] = nnmf(X,nd,'rep',8,'opt',opt,'alg','mult');
    opt = statset('maxiter',1000,'display','final','TolX',1e-7,'TolFun',1e-7);
    [W,H] = nnmf(X,nd,'w0',W0,'h0',H0,'opt',opt,'alg','mult','rep',5);
end

function [X] = normX(X)
    % X should be of form electrode/sentence x time

    % Smooth the data:
    smt = round(0.1*100); % smoothing amount (100ms)    
    for k = 1:size(X,1), X(k,:) = smoothdata(X(k,:),smt); end 
    % Baseline substruction:
    X = bsxfun(@minus,X,mean(X(:, [1:50, end-50:end]), 2)); 
    % Min-Max normalization:    
    X = bsxfun(@minus,X,min(X,[],2))./(max(X,[],2)-min(X,[],2));   
end
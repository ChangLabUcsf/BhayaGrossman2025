% Ilina Bhaya-Grossman
% 01.08.2022
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('util'))
zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;
addpath(genpath(datapath))

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
SIDs = [sSIDs, eSIDs, {'HS8', 'HS9', 'HS10'}]; % , {'HS11', 'HS9', 'HS10'}

timit_details = load('out_sentence_details_timit_all_loudness.mat');
dimex_details = load('out_sentence_details_dimex_all_loudness.mat');
tps = 50:55;

% load sentence responses

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% ---------------- STRF Weight and Prediction Comparison -----------------
%% Comparing cross-language tested STRF predictions
% Takes a minute or two to run.

% Repeated sentences in DIMEx and TIMIT
repsentName = {'s00104','s00804', 's01904', 's03004', 's05004', ...
    's06104', 's06804', 's07404', 's09004', ...
    'fcaj0_si1479', 'fcaj0_si1804', 'fdfb0_si1948', 'fdxw0_si2141', ...
    'fisb0_si2209', 'mbbr0_si2315', 'mdlc2_si2244', 'mdls0_si998', ...
    'mjdh0_si1984', 'mjmm0_si625'};

figure;

% Subjects with no repeated sentences for DIMEx so unable to calculate
% cross-language predictions
SIDs(ismember(SIDs, {'EC252', 'EC152', 'HS8', 'HS9', 'HS10'})) = [];
samelang = 'timit';

if strcmp(samelang, 'dimex')
    predcols = {'b', 'r'}; 
    crosslang = 'timit';
    sentdet = dimex_details.sentdet(ismember({dimex_details.sentdet.name}, ...
        repsentName));
elseif strcmp(samelang, 'timit')
    predcols = {'r', 'b'};  
    crosslang = 'dimex';
    sentdet = timit_details.sentdet(ismember({timit_details.sentdet.name}, ...
        repsentName));
end

% Initialize all the correlation structures
samelang_corr = cell(length(SIDs), 1);
crosslang_corr = cell(length(SIDs), 1);
pred_corr = cell(length(SIDs), 1);

testR_same = nan(height(sent_encoding), 1);
testR_cross = nan(height(sent_encoding), 1);

ctr = 1;
for s = SIDs
    SID = s{1};
    idx = strcmp(sent_encoding.SID, SID);
    modelname = 'onset_aud';

    % Calculating same language predictions
    [out_same, testR] = out_addStrfpred(SID, samelang, modelname, 1, ...
        sentdet, samelang);
    testR_same(idx) = testR(sent_encoding.el(idx));
   
    % Calculating cross language predictions
    [out_cross, testR] = out_addStrfpred(SID, samelang, modelname, 1, ...
        sentdet, crosslang);
    testR_cross(idx) = testR(sent_encoding.el(idx));

    % Plot several example sentence predictions
    if strcmp(SID, 'EC100')
        % [~, maxels] = max(mean(cat(3, testR_cross(1:minel, :), 
        % testR_same(1:minel, :)), [2, 3],'omitnan'));

        % Use the same electrodes from above for cross-prediction example
        maxels = [22 150];
        
        % Plot each electrode example sentences as a new figure
        for maxel = maxels
            figure('renderer', 'painters');

            for sent = 1:3          
                subplot(1, 3, sent);                            

                % Smooth neural response 
                data = smoothdata(squeeze(out_same(sent).resp(maxel, :, :)), ...
                    'gaussian', 'SmoothingFactor', 0);
                x = -0.49:0.01:0.01*size(out_same(sent).resp, 2)-0.5;        
                
                % Plot true neural response
                if size(data, 1) == 1
                    plot(x, data, 'LineWidth', 1.75, 'Color', predcols{1});
                    hold on;
                else
                    data = data - mean(data(1:55, :), [1, 2]);
                    addpath(genpath('../../../ecog_scripts'));
                    shadedErrorBar(x, data', {@mean,@nansem}, ...
                        {'color', predcols{1}, 'linewidth', 1.75, ...
                        'DisplayName', 'response'}, 0.3);
                    hold on;
                end
                
                % Plot predicted responses
                yyaxis right
                samepred = squeeze(out_same(sent).predResp(maxel, :, :))';
                crosspred = squeeze(out_cross(sent).predResp(maxel, :, :))';       

                plot(x, samepred, 'LineStyle', '--', ...
                    'LineWidth', 2, 'Color', predcols{1}, ...
                    'DisplayName', 'same-prediction');
                plot(x, crosspred, 'LineStyle', '--', ...
                    'LineWidth', 2, 'Color', predcols{2}, ...
                    'DisplayName', 'cross-prediction');

                xlim([0.2 1.5]);

                % show predicted response correlations
                [r, ~] = corr(samepred, mean(data, 2));
                text(0, 0.8, ['Same corr :' num2str(r)])

                [r, ~] = corr(crosspred, mean(data, 2));
                text(0, 0.6, ['Cross corr :' num2str(r)])

                yticks([]);
                sentidx = strcmp({sentdet(:).name},  out_cross(sent).name);  
                title(join(sentdet(sentidx).wordList, ' '));
            end                                                 
        end
    end
    ctr = ctr + 1;
end
sent_encoding.([samelang '_trained']) = [testR_same, testR_cross]; 

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Scatter of same-language vs. cross-language STRF predictions

ctr=1;
figure;
lss = [1, 2, 4]; % Spanish, English, Bilingual
for ls = lss
    subplot(1, 3, ctr);

    % Aggregate all same language and cross-language predictions
    lsidx = find(sent_encoding.ls==ls);
    testR_Same = mean([sent_encoding.dimex_trained(:, 1), ...
        sent_encoding.timit_trained(:, 1)], 2, 'omitnan'); 
    testR_Cross = mean([sent_encoding.dimex_trained(:, 2), ...
        sent_encoding.timit_trained(:, 2)], 2, 'omitnan'); 

    % Scatter point of R^2 values
    scatter(testR_Same(lsidx).^2, testR_Cross(lsidx).^2, 15, ...
        [0.6 0.6 0.6], 'filled');
    hold on;

    % Calculate correlation coefficient
    [rho, pval] = corr(testR_Same(lsidx).^2, testR_Cross(lsidx).^2, ...
        'rows', 'pairwise');
    text(0.15, 0.7, ['r=' num2str(rho, 2)], 'FontSize', 15);
    text(0.15, 0.65, ['p=' num2str(pval, 2)], 'FontSize', 15);

    % Formatting
    ylim([0.1 0.75]);
    xlim([0.1 0.75]);   
    xlabel({'Mean same-language','prediction R^2'});
    h = refline(1, 0);
    h.LineWidth = 1.5;
    h.Color = 'k';
    set(gca, 'FontSize', 15);

    if ctr==1
        ylabel({'Mean cross-language','prediction R^2'});
    end

    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Compare STRF weights across corpora

corpora = {'timit', 'dimex'};
%modelnames = {'onset_aud'};
modelnames = {'onset_phnfeatConsOnset_maxDtL_formantMedOnset'};
corrstrf = nan(height(sent_encoding), 1);
corrpval = corrstrf;
maxfreq = nan(height(sent_encoding), 2);
maxrsq = nan(height(sent_encoding), 1);
minrsq = nan(height(sent_encoding), 1);
thresh = 0.01;
% window is 5 time points around the max time point for correlation
windsz = 5;

colors = flipud(brewermap(20, 'Spectral'));
cols = [colors(3, :); colors(end-2, :); colors(round(size(colors, 1)/2), :)];

for s =  unique(sent_encoding.SID)'
    
    % load both corpus STRFs
    SID = s{1};
    if strcmp(modelnames{1}, 'onset_aud')
        meanStrf = nan(2, 81, 61, 256);
    else
        meanStrf = nan(2, 12, 61, 256);
    end
    for c = 1:2
        corpus = corpora{c};
        corpusStrf=loadMultModelStrf(SID, modelnames, corpus, datapath, 1, ...
            'v5');
        if ~isempty(corpusStrf{1})
            minel = size(corpusStrf{1}.meanStrf, 3);
            meanStrf(c, :, :, 1:minel) =   corpusStrf{1}.meanStrf;
        end
    end

    % include weight correlation for the STRFs that explain over 0.01 R^2 
    for e = find(strcmp(sent_encoding.SID, SID))'
         
        el = sent_encoding.el(e);
        missing = or(all(meanStrf(1, :, :, el)==0, [1, 2, 3]), ...
            all(meanStrf(2, :, :, el)==0, [1, 2, 3]));

        if (sent_encoding.en_base_rsq(e)>thresh && ...
            sent_encoding.sp_base_rsq(e)>thresh) && ~missing  
    
            % add in TRF weight correlation
            [~, maxtp] = max(squeeze(mean(meanStrf(:, :, :, el), [1 2])));
            wind = max(1, maxtp-windsz):min(maxtp+windsz, size(meanStrf, 3));
        
            % compare frequency selectivity
            x = squeeze(mean(meanStrf(1, :, wind, el), 3));
            y = squeeze(mean(meanStrf(2, :, wind, el), 3));
            
            % scale and smooth the STRF weights to be more comparable
            [~, z] = procrustes(x', y', 'scaling',true);
            if strcmp(modelnames{1}, 'onset_aud')
                z = smoothdata(z, 'SmoothingFactor',0.4); % , 'SmoothingFactor', 0.7
                x = smoothdata(x, 'SmoothingFactor',0.1);
            end     

            % calculate correlation and max frequency
            [corrstrf(e), corrpval(e)] = corr(x', z, 'type', 'Pearson');
            [~, maxfreq(e, :)] = max([x', z]);
            freq(e, :, :) = [x'; z];

            % save out max R^2 and min R^2
            maxrsq(e) = max(sent_encoding.en_base_rsq(e), ...
                sent_encoding.sp_base_rsq(e));
            minrsq(e) = min(sent_encoding.en_base_rsq(e), ...
                sent_encoding.sp_base_rsq(e));
        
            debug = 0;
            % Example electrodes from above (EC100, 22 / 150)
            if debug && ismember(e, [484, 564])
                % (corrstrf(e)>0.9) && minrsq(e)>0.15
                % ismember(e, [184, 120, 133, 136, 137, 473, 864, 1125, 390])
                % (corrstrf(e)>0.3 && corrstrf(e)<0.6 && minrsq>0.15)

                figure;
                subplot(2, 6, 1);
                plot(x, 1:81, 'LineWidth', 2, 'Color',cols(2, :));
                set(gca, 'XDir', 'reverse');
                ylim([1 81]);
                hold on;                

                % Plot the STRF beta weights for the first language
                subplot(2, 6, [2, 3]);
                imagesc(squeeze(meanStrf(1, :, :, el)));
                set(gca, 'YDir', 'normal');
                xline(wind(1));
                xline(wind(end));
                title('timit');
                
                yticks([1 80]);
                yticklabels({'0.01', '8'});
                ylabel('frequency (kHz)')

                subplot(2, 6, 4);
                plot(z, 1:81, 'LineWidth', 2, ...
                    'LineStyle', '-', 'Color', cols(1, :));
                set(gca, 'XDir', 'reverse');
                hold on; 
                ylim([1 81])
        
                % Plot the STRF beta weights for the second language
                subplot(2, 6, [5, 6]);
                imagesc(squeeze(meanStrf(2, :, :, el)));
                set(gca, 'YDir', 'normal');
                xline(wind(1));
                xline(wind(end));
                title('dimex');
                colormap(inferno);
                yticks([]);
                sgtitle([num2str(corrstrf(e)) ', ' SID ', ' ...
                    num2str(sent_encoding.el(e))]);

                % Plot the erp response
                subplot(2, 6, [ 7, 8, 9 10, 11, 12]);
                plotStitchedSentence(e, sent_encoding, 150, 100, 0.5);
                title(['Language exp: ' num2str(sent_encoding.ls(e))]);
            end  
        end
    end
end

% show pie for correlation threshold
idx = minrsq>0.1;
corrthresh = 0.5;

figure;
subplot(2, 3, [1, 2, 3]);
labels = num2str(crosstab(corrstrf(idx)>corrthresh));
p = pie(crosstab(corrstrf(idx)>corrthresh), [1, 1], labels);
p(3).FaceColor = [159, 134, 192]/256;
p(3).EdgeColor = 'none';
p(1).FaceColor = [224, 177, 203]/256;
p(1).EdgeColor = 'none';
p(2).FontSize = 13;
p(4).FontSize = 13;
legend({'r>0.5', 'r<=0.5'});

for ls = [1, 2, 3]
    subplot(2, 3, ls+2)
    labels = num2str(crosstab(corrstrf(idx&sent_encoding.ls==ls)>corrthresh));
    p = pie(crosstab(corrstrf(idx&sent_encoding.ls==ls)>corrthresh), ...
        [1, 1], labels);
    p(3).FaceColor = [159, 134, 192]/256;
    p(3).EdgeColor = 'none';
    p(1).FaceColor = [224, 177, 203]/256;
    p(1).EdgeColor = 'none';
    p(2).FontSize = 13;
    p(4).FontSize = 13;
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq;


%% Show histogram of STRF weight correlation (run cell before!)
figure;
cols = getColorsCrossComp(1);

subplot(1, 2, 1);
ls = sent_encoding.ls;

lss = [1, 2];

sids = sent_encoding.SID;
normsprsq = rescale(sent_encoding.sp_base_rsq);
normenrsq = rescale(sent_encoding.en_base_rsq);
thresh = 0.05;

% Filter data based on thresholds and conditions
idx = (normsprsq > thresh | normenrsq > thresh) ...
    & ~isnan(corrstrf) & corrpval < 0.01; %  ...
%     & corrstrf > 0

% Scatter plot for normalized R^2 values
for l = lss
    scatter(normsprsq(ls == l & idx), ...
        normenrsq(ls == l & idx), 45, cols(l, :), ...
        'filled', 'MarkerFaceAlpha', 0.6);
    hold on;
end

xlabel('Spanish STRF model R^2 (norm)');
ylabel('English STRF model R^2 (norm)');

% Calculate Pearson correlation and display as text
correl = corr(sent_encoding.sp_base_rsq(ls < 5), ...
    sent_encoding.en_base_rsq(ls < 5), ...
    'type', 'Pearson', 'rows', 'complete');
text(0, 0.45, ['Pearson r=' num2str(correl)], 'FontSize', 15);

% Legend and formattingw
lg = legend({'Spanish', 'English', 'Bilingual'});
lg.Title.String = 'Subject Group';
set(gca, 'FontSize', 15);
ylim([0, 1]);
xlim([0, 1]);
h = refline(1, 0);
h.HandleVisibility = 'off';
h.Color = 'k';
h.LineWidth = 2;

ls = sent_encoding.ls;
ctr = 1;
for l = lss
    subplot(length(lss), 2, ctr * 2)

    % Weighted histogram of beta weight correlation
    [histw, intervals] = histwc(corrstrf(ls == l & idx), ...
        maxrsq(ls == l & idx), 13);
    b = bar(intervals, histw, ...
        'EdgeColor', 'none', 'FaceAlpha', 0.5);
    hold on; 
    b.FaceColor = cols(l, :);

    % Plot density of weighted histogram
    x = 0.25:0.01:0.97;
    y = interp1(intervals, smoothdata(histw, 'SmoothingFactor', 0.5), x, 'spline');
    plot(x, y, 'Color', cols(l, :), 'Linewidth', 2);  

    ctr = ctr + 1;
    set(gca, 'FontSize', 13, 'Ytick', []);
    xlim([0.2 1]);
    box off;
end
ylabel('Weighted Probability');
xlabel('STRF weight corr');




clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq;




%% FIX vis: max time point analysis across languages

figure;
subjtitle = {'Spanish mono', 'English mono'};
for ls = 1:2
    subplot(1, 2, ls)
    rng(2)
    lsidx = find(sent_encoding.ls==ls & ...
        (abs(diff(sent_encoding.maxtp'))>5)' & any(sent_encoding.maxresp'>0.5)');

    randjit = rand(1, length(lsidx))*0.3;
    jittered = [ones(1, length(lsidx))-0.15+randjit; ...
        ones(1, length(lsidx))*2-0.15+randjit];
    scatter([jittered(1, :), jittered(2, :)], [sent_encoding.maxtp(lsidx, 1); sent_encoding.maxtp(lsidx, 2)], ...
        10, [0.2, 0.2, 0.2], 'filled', 'MarkerFaceAlpha', 0.5); hold on;
    line(jittered, ...
        [sent_encoding.maxtp(lsidx, 1) sent_encoding.maxtp(lsidx, 2)]', ...
        'Color', [0.8 0.8 0.8]);

    boxchart(ones(length(lsidx), 1), sent_encoding.maxtp(lsidx, 1), ...
        'MarkerStyle','none', 'BoxLineColor','k', 'BoxFaceAlpha', 0)
    boxchart(ones(length(lsidx), 1)*2, sent_encoding.maxtp(lsidx, 2), ...
        'MarkerStyle','none', 'BoxLineColor','k', 'BoxFaceAlpha', 0);
    [~, p] = ttest(sent_encoding.maxtp(lsidx, 1), sent_encoding.maxtp(lsidx, 2));
    disp(p)
    title(subjtitle{ls})
    xticks([1, 2])
    xticklabels({'english', 'spanish'})
end

% imagesc of sentence responclcses
lsidx = sent_encoding.ls==2;
sgttl = 'English monolingual';
lang = {'english', 'spanish'};

figure;
fields = {'en_sent_resp', 'sp_sent_resp'};

val = cell(1, 2);
sorted_idx = cell(1, 2);
for ls = 1:2
    % sorted by maxtp for their native language
    [val{ls}, sorted_idx{ls}] = sort(sent_encoding.maxtp(lsidx, ls));
    % remove cases where maxtp is before onset
    sorted_idx{ls}(isnan(val{ls}) | val{ls}<30) = [];
    val{ls}(isnan(val{ls}) | val{ls}<30) = [];
end

for ls = 1:2
    subplot(1, 2, ls)   

    % aggregate responses
    resp = [cellfun(@(x) mean(x, 3, 'omitnan'), ...
        sent_encoding.(fields{ls}), 'UniformOutput', false)];
    resp = resp(lsidx);

    resp_tp = max(cellfun(@(x) size(x, 2), resp));
    resp_all = nan(length(resp), resp_tp);
    for i = 1:length(resp) 
        resp_all(i, 1:size(resp{i}, 2)) = resp{i}; 
    end

    % normalize by row?   
    imagesc(rescale(resp_all(sorted_idx{ls}, :))); hold on;
    scatter(val{ls}, 1:length(sorted_idx{ls}), 3, 'k', 'filled');

    if ls==1
        opp = double(~logical(ls-1))+1;
        scatter(val{opp}, 1:length(sorted_idx{opp}), 3, [0.5 0.5 0.5], 'filled');
    end

    xl = xlim();
    xlim([40, xl(2)]);
    colormap(flipud(spectral))
    title(lang{ls})
    clear resp resp_all
end
sgtitle(sgttl);

% correlated peak latency with duration of sentence
details = {timit_details, dimex_details};
resp_field = {'en_sent_resp', 'sp_sent_resp'};
titles = {'TIMIT', 'DIMEx'};
for l = 1:2
    r = nan(height(sent_encoding), 1);
    dur = [details{l}.sentdet.duration];
    for el = 1:height(sent_encoding)
    
        el_resp = padcut(sent_encoding.(resp_field{l}){el}, 3, length(dur));
        out_sent = size(el_resp, 3);
    
        nantrl = squeeze(all(isnan(el_resp)));
        peak_tp = arrayfun(@(x) argmax(squeeze(el_resp(:, :, x))), ...
            1:out_sent);
    
        if sum(nantrl)<out_sent
            r(el) = corr(dur(~nantrl)', peak_tp(~nantrl)', 'Type', 'Spearman');
        end
    end
    
    subplot(1, 2, l)
    lsidx = sent_encoding.ls==1;
    randjit = rand(1, sum(lsidx))*0.3;
    scatter(ones(1, sum(lsidx))-0.15+randjit, r(lsidx), 5, 'filled'); hold on;
    boxchart(ones(sum(lsidx), 1), r(lsidx), ...
            'MarkerStyle','none', 'BoxLineColor','k', 'BoxFaceAlpha', 0); hold on;
    
    lsidx = sent_encoding.ls==2;
    randjit = rand(1, sum(lsidx))*0.3;
    scatter(ones(1, sum(lsidx))*2-0.15+randjit, r(lsidx), 5, 'filled')
    boxchart(ones(sum(lsidx), 1)*2, r(lsidx), ...
        'MarkerStyle','none', 'BoxLineColor','k', 'BoxFaceAlpha', 0);
    ylabel('corr( peak latency, sent duration)');
    xlabel('subject group');
    xticks([1, 2]);
    xticklabels({'Spanish', 'English'})
    title(titles{l})
    ylim([-0.3 0.4])

    % probably should do linear mixed effect model here instead
    [p, ~] = ranksum(r(sent_encoding.ls==1), r(sent_encoding.ls==2));
    disp(p)
end


%% FIX vis: max point of difference between average language responses

sent_diff = arrayfun(@(x) sent_encoding.en_mresp{x}-sent_encoding.sp_mresp{x}, ...
    1:height(sent_encoding), 'UniformOutput', false);
fthresh_diff = arrayfun(@(x) sent_encoding.fvals{x}>sent_encoding.fthresh(x), ...
    1:height(sent_encoding), 'UniformOutput', false);

r = nan(height(sent_encoding), 1);
p = r; tp = r; d = r;

f = zeros(height(sent_encoding), 1);
meandiff = nan(height(sent_encoding), 2);
for i = 1:height(sent_encoding)

    % time points that make the threshold cutoff
    numtp = sum(fthresh_diff{i});

    if numtp>5

        % correlate magnitude difference over time
%         [r(i), p(i)]=corr(abs(sent_diff{i}(fthresh_diff{i}))', ...
%             (1:numtp)', 'type', 'Spearman');
        % look at the average difference in the first 150ms
        meandiff(i, 1) = mean(sent_diff{i}(1:150), 'omitnan');

        % look at the average difference in the last 100ms
        meandiff(i, 2) = mean(sent_diff{i}(150:250), 'omitnan');
    
        debug=0;
        if debug
            figure;
    
            subplot(1, 2, 1)
            scatter(1:numtp, sent_diff{i}(fthresh_diff{i}));
            lsline
            title(r(i));
    
            subplot(1, 2, 2);
            plot(sent_encoding.en_mresp{i}, 'Linewidth', 1.5); hold on;
            plot(sent_encoding.sp_mresp{i});
            title(sent_encoding.ls(i))
        end
    
        % find time point of maximal response difference
        [~, tp(i)]=max(abs(sent_diff{i}));
    
        % find magnitude and direction of response difference 
        % if d>0, then english is larger than spanish
        d(i)=sent_diff{i}(tp(i));
    
        % find whether maximal response time difference has a significant fstat
        if ~isempty(sent_encoding.fvals{i})
            f(i)=sent_encoding.fvals{i}(tp(i))>sent_encoding.fthresh(i);
        end
    end
end

figure;   
cols = getColors(1);
corrected_thresh = (0.001)/height(sent_encoding);
idx = p<corrected_thresh & f>0;
    
scatter3(d(idx), tp(idx), find(idx), 15, sent_encoding.ls(idx), 'filled');  
colormap(cols(1:2, :));
xlabel('Difference (+ eng, - span)');
ylabel('time point of greatest difference');
view(2)


%% FIX Look at peak response differences within single subjects

cols = getColorsCrossComp(1);
maxresp = sent_encoding.maxresp';

% find absolute differences for single subjects and find top 10 percentiles to see 
% the magnitude of the biggest differences
ctr = 1;
unsid = unique(sent_encoding.SID);
sidiff = nan(length(unsid), 2);
sidls = nan(length(unsid), 1);
sidiffall = nan(length(unsid), 100);
for s = unique(sent_encoding.SID)'

    SID = s{1};
    els = find(strcmp(sent_encoding.SID, SID));

    if length(els)>5
        % find all max resps for this subject
        mresp = sent_encoding.maxresp(els, :);
        mdiff = diff(mresp')/max(abs(mresp), [], 'all')';
        sidiffall(ctr, 1:length(mdiff)) = mdiff;
    
        debug = 0;
        if debug
            figure; 
            histogram(mdiff); hold on;
            % negative (5%) is english>spanish
            xline(prctile(mdiff, 5), 'linewidth', 1.5);
            % positive (95%) is english<spanish
            xline(prctile(mdiff, 95));
        end    
        
        % find top 10 perctile point   
        sidiff(ctr, 1) = prctile(mdiff, 5);
        sidiff(ctr, 2) = prctile(mdiff, 95);
        sidls(ctr) =  sent_encoding.ls(find(strcmp(sent_encoding.SID, SID), ...
            1, 'first'));
        
    end
    ctr=ctr+1;
end


figure;
binrng = -0.5:0.05:0.5;
counts = nan(size(sidiffall, 1), length(binrng)-1);
for s = 1:size(sidiffall, 1)
    counts(s, :) = histcounts(sidiffall(s, :), binrng);
end
counts(all(counts==0, 2), :) = [];
imagesc(counts, 'Xdata', binrng);
xline(prctile(sidiffall(:), 5), 'LineWidth', 1.5)
xline(prctile(sidiffall(:), 95), 'LineWidth', 1.5)
colormap([1 1 1; crest]);
colorbar;
ylabel('subject');
xlabel('% difference');
set(gca, 'FontSize', 15);

peakdiffs = diff(maxresp)/max(maxresp, [], 'all');
figure;
ctr = 1;
for ls = [1, 2]
    subplot(3, 1, ctr);
    histogram(peakdiffs(sent_encoding.ls==ls), 'FaceColor', ...
        cols(ls , :), 'EdgeColor', 'none', ...
        'Normalization', 'Probability'); hold on;
    yticks([]);
    ylabel('probability');
    yyaxis right

    [f, xi] = ksdensity(peakdiffs(sent_encoding.ls==ls));
    plot(xi, f, 'Color', cols(ls, :), 'Linewidth', 2);
    yticks([])
    xlim([-0.12 0.12]);
    xticks(-0.1:0.1:0.1)
    [~, locs] = findpeaks(f, 'NPeaks',1, 'SortStr','descend');
    xline(xi(locs(1)), 'LineWidth', 2);
    box off;
    set(gca, 'FontSize', 13);
    
    yticks([]);

    ctr = ctr+1;
end


figure;
cols = getColorsCrossComp(1);
ctr = 1;
for ls = [1, 2, 4]
    lidx = sidls==ls;

    subplot(1, 2, 1);
    scatter(rand(sum(lidx), 1)*0.2+ctr-0.1, ...
        abs(sidiff(lidx, 1)), 35, cols(ls, :), ...
        'LineWidth', 1.5); hold on;
    title('En > Sp (%)')

    subplot(1, 2, 2);
    scatter(rand(sum(lidx), 1)*0.2+ctr-0.1, ...
        abs(sidiff(lidx, 2)), 35, cols(ls, :), ...
        'LineWidth', 1.5); hold on;
    title('En < Sp (%)')
    ctr = ctr+1;
end

for i = 1:2
    subplot(1, 2, i);
    ylabel('% diff rel. to max');
    xlabel('subject group');
    xticks(1:3);
    xticklabels({'Spanish', 'English', 'Bilingual'});
    set(gca, 'FontSize', 13);
    xlim([0 4]);
%     ylim([0 0.351])
end


cols = getColorsCrossComp(1);
% [r, c] = gridSize(length(unique(sent_encoding.SID)));
figure;
ctr=1;
[~, i] = sort(sent_encoding.ls); 
for SID = unique(sent_encoding.SID(i), 'stable')'
    
    if sum(strcmp(sent_encoding.SID, SID))>10
        subplot(4, 6, ctr);
        ls = sent_encoding.ls(find(strcmp(sent_encoding.SID, SID), 1));

        scatter(maxresp(2, strcmp(sent_encoding.SID, SID)), ...
            maxresp(1, strcmp(sent_encoding.SID, SID)), 24, ...
            cols(ls, :), 'filled', 'MarkerFaceAlpha', 0.7);    
        
        ylim([0, 4]);
        xlim([0, 4])
        refline(1, 0);
        title(SID);
        yticks([]);
        xticks([]);
        if ctr==1
            xlabel('Spanish HFA');
            ylabel('English HFA')
        end
        ctr=ctr+1;
    end
end

%% FIX: Determine whether peaks differ significantly within versus across language

% p = nan(1, height(sent_encoding));
% p_en = nan(1, height(sent_encoding));
% p_sp = nan(1, height(sent_encoding));
% nsamples = 100;
% for i = 1:height(sent_encoding)
%     if ~strcmp(sent_encoding.SID(i), 'EC219')
%         en = sent_encoding.en_sent_resp{i};
%         sp = sent_encoding.sp_sent_resp{i};
%         maxtp = sent_encoding.maxtp(i);
%     
%         if ~isnan(maxtp)
%             en = squeeze(en(1, maxtp, :));
%             sp = squeeze(sp(1, maxtp, :));
%     
%             for j = 1:20
%                 d(i, j) = mean(randsample(en, nsamples))- mean(randsample(sp, nsamples));
%                 d_en(i, j) = mean(randsample(en, nsamples))- mean(randsample(en, nsamples));
%                 d_sp(i, j) = mean(randsample(sp, nsamples))- mean(randsample(sp, nsamples));
%             end
%             rng(1)
%             [~, p(i)] = ttest2(randsample(en, nsamples), randsample(sp, nsamples));
%             [~, p_en(i)] = ttest2(randsample(en, nsamples), randsample(en, nsamples));
%             [~, p_sp(i)] = ttest2(randsample(sp, nsamples), randsample(sp, nsamples));
%         end
%     end   
% end

typeinfo = struct();
typeinfo.sigField = 'language';
typeinfo.pthresh=0.001;
typeinfo.window = -1:10; % time points that need to be significant
mintp = 300;
data = struct();
mintrl = 100;
data.language = [ones(1, mintrl), ones(1, mintrl)*2];
tmp = nan(height(sent_encoding), mintp, mintrl*2);
for i = 1:height(sent_encoding)
    
    % splitting it up because each subject has a different number of trials
    % electrode by time by trial
    en = cat(1, sent_encoding.en_sent_resp{i});
    sp = cat(1, sent_encoding.sp_sent_resp{i});
        
    %mintp= min([size(en, 2), size(sp, 2)]);
    %disp([num2str(size(en, 3)) ':' num2str(size(sp, 3))] )
    
    nanen = find(~any(isnan(en), [1, 2]));
    nansp = find(~any(isnan(sp), [1, 2]));
    
    if all([length(nanen), length(nansp)]>mintrl)
        tmp(i, :, 1:mintrl) = en(:, 1:mintp, nanen(1:mintrl));
        tmp(i, :, mintrl+1:mintrl*2) = sp(:, 1:mintp, nansp(1:mintrl));
    end
end
data.all.resp = tmp;
[allidx, ~] = getElecs(data, {'all'}, 'bychan', '', 'ftest', typeinfo);
sent_encoding.sig = zeros(height(sent_encoding), 1);
sent_encoding.sig(allidx.all) = 1;
sent_encoding.sig = logical(sent_encoding.sig);

ctr=1;
figure;
for ls = [1, 2, 4]
    subplot(1, 3, ctr);
    lsidx = sent_encoding.ls==ls;
    idx = intersect(find(lsidx), allidx.all);

    sids = cellfun(@(x) str2double(x(3:end)), sent_encoding.SID(idx));
    tsids = unique(sent_encoding.SID(lsidx));
    sidnum = arrayfun(@(x) find(unique(sids)==x), sids);
    
    scatter3(sent_encoding.maxresp(idx, 1), ...
        sent_encoding.maxresp(idx, 2), idx, 25, sidnum, ...
        'filled');
    colormap(brewermap(8, 'Dark2'));
    title([num2str(length(idx)./sum(lsidx), 2) ': ' ...
        num2str(length(unique(sids))) ' - ' num2str(length(tsids))]);
    view(2);
   
    ctr=ctr+1;

    ylim([0, 3]);
    xlim([0, 3]);
    h=refline(1, 0);
    h.LineWidth=1.75;
    h.Color = 'k';
end 

% NNMF on sentence level responses
sent = 200;
X = mean(data.all.resp(:, :, 1:100), 3, 'omitnan');
lidx = sent_encoding.ls;
nanrow = sum(isnan(X), 2)>0 | lidx==2; % nonnati

X(nanrow, :)=[];

lidx(nanrow)=[];

rng(2);
clusters = 4;
[W,H] = nnmf(X, clusters, 'rep', 10);
figure; plot(H', 'LineWidth', 2.5);
legend({'Cluster 1', 'Cluster 2', 'Cluster 3', 'Cluster 4', 'Cluster 5'});
 
figure;
for clust = 1:clusters
    subplot(1, clusters, clust); 
    othclust = 1:clusters;
    othclust(othclust==clust)=[]; 
    [sort_val, sort_idx] = sort(W(:, clust));

    idx = (W(:, clust)-max(W(:, othclust), [], 2))<3 ;
    sort_idx(idx) = [];
    sort_val(idx) = [];

    all_X = mean(data.all.resp, 3, 'omitnan');
    all_X(nanrow, :)=[];
    imagesc(all_X(sort_idx(:), 40:end));
    colormap(flipud(spectral))
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%%

%% vis: comparing speech responsive in TIMIT vs DIMEx

binedges = -0.20:0.02:0.20;
colors = flipud(brewermap(length(binedges)-1, 'Spectral'));
cols = [colors(3, :); colors(end-2, :); colors(round(size(colors, 1)/2), :)];

timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

shared = nan(length(fields(timit_elecs.allidx)), 4);
ctr = 1;
for s = fields(timit_elecs.allidx)'
    SID = s{1};
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs}));
    if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
        dimels = dimex_elecs.allidx.(SID);
        timels = timit_elecs.allidx.(SID);
        shared(ctr, :) = [length(setdiff(dimels, timels)) ...
            length(setdiff(timels, dimels)) ...
            length(intersect(dimels, timels)) ls];        
    end
    ctr = ctr+1;
end

speech_elec = nan(height(sent_encoding), 2);
for i = 1:height(sent_encoding)
    SID = sent_encoding.SID{i};
    if isfield(timit_elecs.allidx, SID)
        if ismember(sent_encoding.el(i), timit_elecs.allidx.(SID)) 
            speech_elec(i, 1) = 1;
        else
            speech_elec(i, 1) = 0;
        end
    end

    if isfield(dimex_elecs.allidx, SID)
        if ismember(sent_encoding.el(i), dimex_elecs.allidx.(SID)) 
            speech_elec(i, 2) = 1;
        else
            speech_elec(i, 2) = 0;
        end
    end
end

% resp 
figure;
labels = {'Spanish mono', 'English mono', 'Bilingual'};
ctr = 1;
for ls = [1, 2, 4]
    subplot(1, 3, ctr);
    p = pie(sum(shared(shared(:, 4)==ls, 1:3)), [1 1 1 ], '%.f%%');

%     colormap([0.19 0.53 0.74; 0.83 0.24 0.30; 0.87 0.77 0.9]);
    colormap(cols);
    p(1).FaceAlpha = 0.4;
    p(3).FaceAlpha = 0.4;
    p(5).FaceAlpha = 0.4;
    title(labels{ctr})
    ctr = ctr+1;
end

figure;
tmp = nan(height(sent_encoding), 1);
tmp(all(speech_elec')) = 3;
nanidx = any(isnan(speech_elec), 2);
tmp(speech_elec(:, 1)==1 & speech_elec(:, 2)==0) = 2;
tmp(speech_elec(:, 1)==0 & speech_elec(:, 2)==1) = 1;
scatter3(rescale(sent_encoding.en_base_rsq(~nanidx)), ...
    rescale(sent_encoding.sp_base_rsq(~nanidx)), find(~nanidx), ...
    45, tmp(~nanidx), 'filled', 'MarkerEdgeColor', 'k', ...
    'MarkerFaceAlpha', 0.6);
colormap(cols);
h = refline(1, 0);
h.LineWidth = 1.5;
h.Color = 'k';
grid off;
view(2)
set(gca, 'FontSize', 15);
ylabel('spanish encoding R^2');
xlabel('english encoding R^2');

% add the rsq values too

%% TOFIX: visualizing result (MNI brain with TRF R^2) 

colors = flipud(brewermap(50, 'Spectral'));
speechrspns = sent_encoding;

bins = 13;
[speechrspns.cond, edges] = discretize(max(speechrspns.en_base_rsq, ...
    speechrspns.sp_base_rsq), bins);

% initialize design electrode structure
desel=struct();
desel.conds = 1:bins;
desel.sz = [30 repmat(50, 1, length(desel.conds)-1)];
desel.cols = colors(9:3:end, :);

% legend
scatter(1:length(desel.cols), 1:length(desel.cols), 75, desel.cols, 'filled');
colorbar('XTick',[0 1], 'XTickLabel', {'0.05', '0.45'});
colormap(desel.cols);

% split up peak rate and phonetic features again for MNI plotting
desel.labels = [];
for s=unique(speechrspns.SID)'
    SID = s{1};
    idx = strcmp(speechrspns.SID, SID);
    desel.(SID).elid = speechrspns.el(idx);
    desel.(SID).condition = speechrspns.cond(idx);
end

% plot histogram of count across y-axis (anterior to posterior)
for ls = 1:3    
    lsid = find(speechrspns.ls == ls);
    mni_rh = [];
    if ls<3 
        [mni_lh] = plotMNIElec(unique(speechrspns.SID(lsid)), desel, 'lh', 0);
        [mni_rh] = plotMNIElec(unique(speechrspns.SID(lsid)), desel, 'rh', 0);
    else
         % no RH HS subjects
        [mni_lh] = plotNativeElec({'HS11'}, desel);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* encodings *D*;

%% frequency selectivity average

modelnames={'onset_phnfeatonset_maxDtL', 'onset_aud'};   % 'onset_phnfeatonset_maxDtL'{'onset_maxDtL_aud'};   % 
corpus = {'timit', 'dimex'};
details = {timit_details, dimex_details};

aud_weights = nan(height(sent_encoding), 80); % max tp frequency selectivity
feat_weights = nan(height(sent_encoding), 14);

word_rsq = [];
rsq_idx = [];
for s = unique(sent_encoding.SID)'
    SID = s{1};
    corpusStrf = loadMultModelStrf(SID, modelnames, corpus{2}, ...
        datapath, 1);
    c = 1;
    if ~isempty(corpusStrf{2})

        els = strcmp(sent_encoding.SID, SID);
        word_rsq = [word_rsq corpusStrf{1}.meanTestR(sent_encoding.el(els)).^2];
        rsq_idx = [rsq_idx find(els)'];

        for row = find(els)'
            el = sent_encoding.el(row);
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

            if corpusStrf{2}.meanTestR(el)^2>0.0001
                aud_weights(row, :)  = y;
            end
        end
    end
end

figure;
[g, ~, ~, labels] = crosstab(max_feat);

x = cellfun(@(x) str2double(x), labels);
b = bar(x, g);
b(1).FaceColor = [0.5 0.5 0.5];

xticks(1:14);
xlim([0.4 14.8]);
xticklabels(featnames);

ylabel('Electrode count');
xlabel('Primary Phonetic Feature');
set(gca, 'FontSize', 15);
box off;

% r-squared comparison
figure;
scatter(sent_encoding.en_base_rsq(rsq_idx, 1), word_rsq, 55, 'k', ...
     'filled', 'MarkerFaceAlpha', 0.5); 

ylim([0 0.55]);
xlim([0 0.55]);

h=refline(1, 0);
h.LineWidth = 1.5;
h.LineStyle = '-';
h.Color = 'k';
ylabel('phonetic feature R^2');
xlabel('mel spectrogram R^2');
set(gca, 'FontSize', 15);

% beta weight comparison (what does the average frequency selectivity look like?)
figure;  
ctr = 1;
subject = {'Spanish', 'English', 'Bilingual'};
for ls = unique(sent_encoding.ls)'
    subplot(1, length( unique(sent_encoding.ls)), ctr)
    
    audls = aud_weights(sent_encoding.ls == ls, :);
    nanidx = any(isnan(audls)');
    audls = audls(~nanidx, :);

    [~, tpidx] = max(audls, [], 2);
    [~, sortidx] = sort(tpidx);
    imagesc(audls(sortidx, :)); %.*rhos(v_minus(sortidx)));
    ylabel('Electrode count');
    xlabel('Frequency Bin');
    set(gca, 'FontSize', 15);
    box off;
    title(subject(ctr))
    ctr = ctr+1;
end

cm = brewermap(50, 'Greys');
colormap(cm(1:end-8, :));

cbh = colorbar();
ylabel(cbh, 'STRF weight (% max)');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *cons*;

%% --------------- Frequency Selectivity: STRFs -----------------------------------------

modelname = {'onset_aud'};
SIDs = [sSIDs bSIDs eSIDs];

figure;
bimola = struct();
for s = SIDs
    SID = s{1};
    corpusStrf{1}=loadMultModelStrf(SID, modelname, 'dimex', datapath, 1, 'v5');   
    corpusStrf{2}=loadMultModelStrf(SID, modelname, 'timit', datapath, 1, 'v5');  
    
    if ~isempty(corpusStrf{1}{1}) && ~isempty(corpusStrf{2}{1})
        els = 1:min(length(corpusStrf{1}{1}.Els), length(corpusStrf{2}{1}.Els));
        bimola.(SID).corr = arrayfun(@(e) corr(reshape(corpusStrf{1}{1}.meanStrf(:, :, e), ...
            [], 1), reshape(corpusStrf{2}{1}.meanStrf(:, :, e), [], 1)), els);       
        bimola.(SID).rsq(1, els) = corpusStrf{1}{1}.meanTestR(els).^2;
        bimola.(SID).rsq(2, els) = corpusStrf{2}{1}.meanTestR(els).^2;
        clear corpusStrf

        idx = all(bimola.(SID).rsq>0.01, 1); % els(idx), repmat(str2num(SID(3:end)), 1, sum(idx))
        scatter3(bimola.(SID).rsq(1, idx), bimola.(SID).rsq(2, idx),...
            repmat(str2double(SID(3:end)), 1, sum(idx)) , 55, ...
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

% comparison of correlation between English, Spanish, Bilingual
% onset_maxDtL_maxDtLOnset_vowelOnset
modelname = {'onset_aud'};

varnames = {'SID', 'el', 'lang', 'corr', 'srsq', 'ersq'};
bimola = array2table(zeros(0,6), 'VariableNames', varnames);
for SIDs = {sSIDs, eSIDs, bSIDs}
    for s = SIDs{1}
        SID = s{1};
        corpusStrf{1}=loadMultModelStrf(SID, modelname, 'dimex', ...
            datapath, 1, 'v5');   
        corpusStrf{2}=loadMultModelStrf(SID, modelname, 'timit', ...
            datapath, 1, 'v5');  

        if ~isempty(corpusStrf{1}{1}) && ~isempty(corpusStrf{2}{1})
            els = 1:min(length(corpusStrf{1}{1}.Els), length(corpusStrf{2}{1}.Els));
            
            rsq = nan(2, length(els));
            correl = nan(1, length(els));
            ctr=1;            
            for e = els
                [~, tp1]=max(mean(corpusStrf{1}{1}.meanStrf(:, :, e)));
                [~, tp2]=max(mean(corpusStrf{1}{1}.meanStrf(:, :, e)));
                correl(ctr) = corr(reshape(corpusStrf{1}{1}.meanStrf(:, tp1, e), ...
                    [], 1), reshape(corpusStrf{2}{1}.meanStrf(:, tp2, e), [], 1)); 
                ctr=ctr+1;
            end
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
% cols = [0.6000 0.4392 0.6706; 0.3529 0.6824 0.3804; 0.8784 0.5098 0.0784];
% violin({bimola.corr(bimola.lang==1), ...
%     bimola.corr(bimola.lang==2), bimola.corr(bimola.lang==3)}, 'facecolor', ...
%     cols, 'medc', []);

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
    betaInfo* *encoding*;





%% --------------- Primary Encoding: TRFs --------------------------------

% using phonetic feature instead of splitting by vowel
% modelname={'onset_phnfeatConsOnset_maxDtL_formant'};  
modelname={'onset_phnfeatConsOnset_maxDtL_formantMedOnset'};  

% 'consfeat', 'formant', 
lang = {'eng', 'sp'};
corpus = {'timit', 'dimex'};

timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

% determine unique variance per feature and primary encoding
varnames = {'SID', 'el', 'ls', ...
    [lang{1} '_base_rsq'], ...
    [lang{1} '_base_beta'], ...
    [lang{2} '_base_rsq'], ...
    [lang{2} '_base_beta'], ...
    };
phn_encoding =  array2table(zeros(0, length(varnames)), 'VariableNames', varnames);

[sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();
for s = [sSIDs eSIDs bSIDs mSIDs]
    SID = s{1}; 
    ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs})); %  mSIDs
    for l = 1:2
        corpusStrf{l} = loadMultModelStrf(SID, modelname, corpus{l}, ...
            datapath, 1, 'v5');
    end

    if ~any(cellfun(@(x) isempty(x), [corpusStrf{1} corpusStrf{2}]))

        % find minimum test R
        minel = min(cellfun(@(x) length(x{1}.meanTestR), corpusStrf));
        if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            els = intersect(timit_elecs.allidx.(SID), dimex_elecs.allidx.(SID));
        elseif isfield(dimex_elecs.allidx, SID)
            els = dimex_elecs.allidx.(SID);
        elseif isfield(timit_elecs.allidx, SID)
            els = timit_elecs.allidx.(SID);
        end

        base = cell(2, 1);
        full = cell(2, 1);
        uvall = cell(2, 1);
        fullBeta = cell(2, 1);

        for l = 1:2
            % full models
            base{l} = (corpusStrf{l}{1}.meanTestR.^2)';
            fullBeta{l} = corpusStrf{l}{1}.strf;

        end

        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);

        for l = 1:2
            betaCell{l} = squeeze(mat2cell(fullBeta{l}{1}(:, :, els), ...
                size(fullBeta{l}{1}, 1), size(fullBeta{l}{1}, 2), ones(length(els),1)));
        end


        tmp = table(sids, els, lss, ...
            base{1}(els), betaCell{1}, ...
            base{2}(els), betaCell{2}, ...
            'VariableNames', varnames);

        phn_encoding = [phn_encoding; tmp];
    else
        warning(['Missing subject ' SID]);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *wrd;

%% get maximum betas and compare distribution across subjects
% id = [52, 258, 671, 800, 896];
feats = {'onset'};
feats = [feats; timit_details.features.names([1:3, 8, 9, 11])];
feats = [feats; {'peakRate', 'F1', 'F2', 'F3', 'F4', 'None'}'];
%[4:7 10]
minthresh = 0.1;
maxthresh = 0.1;
meethresh = all([phn_encoding.eng_base_rsq; ...
    phn_encoding.sp_base_rsq]<maxthresh, 2) & ...
    all([phn_encoding.eng_base_rsq; ...
    phn_encoding.sp_base_rsq]>minthresh, 2);

idx = nan(height(phn_encoding), 2, 3);
fields = {'eng_base_beta', 'sp_base_beta'};

for i = 1:height(phn_encoding)
    if meethresh(i), figure; end
    for l = 1:2
        betas = smoothdata(phn_encoding.(fields{l}){i}, 2);              

        % removing negative betas
        betas(abs(betas)<prctile(abs(betas(:)), 95))=0;

        if meethresh(i)
            ax=subplot(1, 2, l);
            
            imagesc(betas);
            yticks(1:15);
            yticklabels(feats);
            colormap(ax, centered('RdBu', 20));            
        end
        [vals, idx(i, l, :)] = maxk(mean(abs(betas), 2), 3);
        idx(i, l, vals==0)=13;
    end
    
end

corp = {'English', 'Spanish'};
fields = {'eng_base_rsq', 'sp_base_rsq'};
subj = {'Spanish', 'English', 'Mandarin', 'Bilingual'};
for l = 1:2
    figure;
    for ls= [1, 2, 4]
        lsid=phn_encoding.ls==ls & phn_encoding.(fields{l})>0.015;
        lsidx = squeeze(idx(lsid, l, :));
        labels = cell(size(lsidx, 1), 2);
        for c = 1:3
            for i = 1:length(feats)
                if ismember(i, 2:7)
                    labels(ismember(lsidx(:, c), i), c) = {'consonant'};
                elseif ismember(i, 9:12)
                    labels(ismember(lsidx(:, c), i), c) = {'vowel'};
                else
                    labels(ismember(lsidx(:, c), i), c) = feats(i);
                end
                
            end
        end
        
        betaTbl = table();
        betaTbl.first = labels(:, 1);
        betaTbl.second = labels(:, 2);
        betaTbl.third = labels(:, 3);

        subplot(4, 3, find(ls==[1, 2, 4]));
        histogram(categorical([betaTbl.first]), 'FaceColor', [0.5 0.5 0.5]);
        title(subj{ls})

        subplot(4, 3, 3+find(ls==[1, 2, 4]));
        histogram(categorical([betaTbl.second]), 'FaceColor', [0.5 0.5 0.5]);

        ax=subplot(4, 3, [6+find(ls==[1, 2, 4]), 9+find(ls==[1, 2, 4])]);
        heatmap(betaTbl, 'first', 'second');
%         heatmap(betaTbl, 'first', 'third');
    
        colormap(flipud(brewermap(30, 'greys')));
        caxis([0 70])
        title(subj{ls});
    end
    sgtitle(corp{l})
end


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;


%% ------------- Speech Responsive: Feature TRFs -------------------------------
%% ---------------- NEW: TRF surprisal unique variance -------------------

% using phonetic feature instead of splitting by vowel
% to do: add noOns surprisal feature

modelnames_timit={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... %remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... %remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_engSurpNoOnsBin', ... %remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... % remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin'}; % full
%'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... %remove word

modelnames_dimex={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_spSurpNoOnsBin', ... %remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_spSurpNoOnsBin', ... % remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin'}; % full

% modelnames_timit={'phnfeatonset_maxDtL_wordOns_wordL_engSurpBin', ... % no onset feat
%     'onset_maxDtL_wordOns_wordL_engSurpBin', ... % no consonant  feat
%     'onset_phnfeatonset_wordOns_wordL_engSurpBin', ... % no peakrate feat
%     'onset_maxDtL_wordOns_wordL_engSurpBin', ... % no formant feat
%     'onset_phnfeatonset_maxDtL', ... % no word/surp feat (base)
%     'onset_phnfeatonset_maxDtL_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatonset_maxDtL_wordOns_wordL_engSurpBin', ... % full model feat
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_engSurpBin' % larger full model
%     };  
% 
% modelnames_dimex={'maxDtL_phnfeatonset_wordOns_wordL_spSurpBin', ... % no onset feat
%     'onset_maxDtL_wordOns_wordL_spSurpBin', ... % no consonant  feat
%     'onset_phnfeatonset_wordOns_wordL_spSurpBin', ... % no peakrate feat
%     'onset_maxDtL_wordOns_wordL_spSurpBin', ... % no formant feat
%     'onset_phnfeatonset_maxDtL', ... % no word/surp feat (base)
%     'onset_phnfeatonset_maxDtL_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatonset_maxDtL_wordOns_wordL_spSurpBin', ... % full model feat
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_spSurpBin' % larger full model
%     }; 

% modelnames_timit={'phnfeatConsOnset_maxDtL_formant_wordOns_wordL_engSurpBin', ... % no onset feat
%     'onset_maxDtL_formant_wordOns_wordL_engSurpBin', ... % no consonant  feat
%     'onset_phnfeatConsOnset_formant_wordOns_wordL_engSurpBin', ... % no peakrate feat
%     'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_engSurpBin', ... % no formant feat
%     'onset_phnfeatConsOnset_maxDtL_formant', ... % no word/surp feat (base)
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_engSurpBin', ... % full model feat
%     };  
% 
% modelnames_dimex={'phnfeatConsOnset_maxDtL_formant_wordOns_wordL_spSurpBin', ... % no onset feat
%     'onset_maxDtL_formant_wordOns_wordL_spSurpBin', ... % no consonant  feat
%     'onset_phnfeatConsOnset_formant_wordOns_wordL_spSurpBin', ... % no peakrate feat
%     'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_spSurpBin', ... % no formant feat
%     'onset_phnfeatConsOnset_maxDtL_formant', ... % no word/surp feat (base)
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL', ... % no surp feat   
%     'onset_phnfeatConsOnset_maxDtL_formant_wordOns_wordL_spSurpBin', ... % full model feat
%     }; 

SIDs = [sSIDs eSIDs bSIDs mSIDs];
[wordsurp_encoding] = loadUniqueVarTbl(modelnames_timit, modelnames_dimex, SIDs);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *wrd;

%% vis: feature TRF weights + neural response predictions

% id = [52, 258, 671, 800, 896];
feats = {'onset'};
feats = [feats; timit_details.features.names([1:3, 8, 9, 11])];
feats = [feats; {'peakRate', 'F1', 'F2', 'F3', 'F4', 'word onset', 'word length', ...
    'engSurpBin'}'];

% % english
elecs = containers.Map;
% elecs('EC163') = 8; % , 105, 106, 185
% elecs('EC105') = [166, 168, 106, 169];
% elecs('EC183') = [55, 71, 72];
% elecs('EC222') = [89, 90];

elecs('EC183') = [71];
% spanish
% elecs('EC172') = [105, 122];
% elecs('EC100') = 118;

% thresh = 0.15;
% meethresh = all([wordsurp_encoding.eng_full_rsq; ...
%     wordsurp_encoding.sp_full_rsq]>thresh, 2);

% now threshold on unique variance instead of full rsq
thresh = 0.03;
meethresh = any([wordsurp_encoding.eng_uv_all(:, 5), ...
    wordsurp_encoding.sp_uv_all(:, 5)]>thresh, 2);

idx = nan(height(wordsurp_encoding), 5);
fields = {'eng_full_beta', 'sp_full_beta'};
fieldrsq = {'eng_full_rsq', 'sp_full_rsq'};
corpus = {'timit', 'dimex'};

basemdl = 'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_';
model = {[basemdl 'engSurpNoOnsBin'], [basemdl 'spSurpNoOnsBin']};
surpfeat = 15;

for i = 1:height(wordsurp_encoding)

    SID = wordsurp_encoding.SID{i};
    el = wordsurp_encoding.el(i);

    % if plotting
    if (isKey(elecs, SID) && ismember(el, elecs(SID))) % || meethresh(i)
        desel.(SID).selid = el;
        desel.(SID).elid = 1:256;
        desel.(SID).condition = ones(1, 256);
        desel.conds = 1;

        desel.sz = 25;
        desel.cols  = [0 0 0];
        desel.labels = {''};
        plotNativeElec({SID}, desel, 1);

        figure('Renderer', 'Painters'); 
    end

    for l = 1:2
        betas = smoothdata(wordsurp_encoding.(fields{l}){i}, 2);  
        surpbetas = betas(surpfeat:surpfeat+7, :);

        % removing negative betas
        betas(abs(betas)<prctile(abs(betas(:)), 95))=0;

        if  (isKey(elecs, SID) && ismember(el, elecs(SID))) % || meethresh(i)
            % plot the sentence level response for a repeat sentence
            subplot(2, 4, [l+(l-1) l+l]);
            
            plotTRFPred(wordsurp_encoding.SID{i}, wordsurp_encoding.el(i), ...
                corpus{l}, model{l}, dimex_details, timit_details);
            xticks(0:0.5:1.5);
            xlabel('Time (s)');
            yyaxis left
            yticks([]);
            ylabel('HFA (z)');
            set(gca, 'FontSize', 13)
            box off;
            title(['R^2=' num2str(wordsurp_encoding.(fieldrsq{l})(i), 3)]);
 
            % plot sparse-ified the beta weights 
            ax=subplot(2, 4, l+4+(l-1));          
            imagesc(betas);
            yticks(1:15);
            yticklabels(feats);
            colormap(ax, centered);   
            ylim([0.5 14.5]);
            set(gca, 'FontSize', 13)
            xticks([1 60]);
            xticklabels({0 -0.6});
            xlabel('Delay (s)')
            box off;
    
            % plot surprise weights 
            ax=subplot(2, 4, l+5+(l-1));     
            surpcols = flipud(brewermap(8, 'Spectral'));
            % find maximal time point (past first five tps)
            [~, tp] = max(mean(abs(surpbetas(:, 5:55)), 1));
            tp = tp+5;

            y = mean(surpbetas(:, max(1, tp-5):min(tp+5, 61)), 2);
            scatter(1:8, y, 35, 'k', 'filled');hold on;
            %surpcols
            fitln=lsline();
            fitln.LineWidth = 1.8;  
            [r, p] = corr((1:8)', y, 'type', 'Spearman');                   
            text(1, 0.1, {['r=' num2str(r, 3) ],[' p=' num2str(p, 2)]});
            set(gca, 'FontSize', 13);
            xlim([0.5 8.5])
            yticks(-0.2:0.1:0.2);
            xticks([1 8]);
            ylabel('beta weight')
            xticklabels({'min', 'max'});
            xlabel('surprise')

            sgtitle([ SID ' el: ' num2str(el) ...
                ' lang: ' num2str(wordsurp_encoding.ls(i))])

            % plot surprise weight dynamics
            debug = 0;
            if debug   
                figure;
                for j = 1:8
                    x = -0.6:0.01:0;
                    plot(x, flipud(smoothdata(surpbetas(j, :)')), ...
                        'Color', surpcols(j, :), 'LineWidth', 2); hold on;   
                    xline(x(61-tp));
                end
                ylabel('beta weight');
                xlabel('time (s)');
                ylim([-0.25 0.25]);
                set(gca, 'FontSize', 15);
                box off;
            end
        end
        [~, idx(i, :)] = maxk(mean(betas, 2), 5);
    end    
end
labels = cell(size(idx, 1), 2);

for i = 1:5
    labels(ismember(idx(:, i), 1), i) = {'onset'};
    % consonant collapse
    labels(ismember(idx(:, i), 2:7), i) = {'consonant'};
    labels(ismember(idx(:, i), 8), i) = {'peakRate'};
    % vowel collapse
    labels(ismember(idx(:, i), 9:12), i) = {'vowel'};
    labels(ismember(idx(:, i), 13), i) = {'word boundary'};
    labels(ismember(idx(:, i), 14), i) = {'word length'};
    labels(ismember(idx(:, i), 15:22), i) = {'surprise'};
end

betaTbl = table();
betaTbl.first = labels(:, 1);
betaTbl.second = labels(:, 2);
betaTbl.third = labels(:, 3);
figure; heatmap(betaTbl, 'first', 'second');
% figure; heatmap(betaTbl, 'first', 'third');

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% vis: similar phonetic encoding across language contexts

% el = 196; %196
feats = [{'onset'}; timit_details.features.names([1:3, 8, 9, 11])];
feats = [feats; {'peakRate', 'F1', 'F2', 'F3', 'F4'}'];

c = nan(height(wordsurp_encoding), 1);
for el = 1:height(wordsurp_encoding)
    x = wordsurp_encoding.eng_full_beta{el}(1:12, :);
    y = wordsurp_encoding.sp_full_beta{el}(1:12, :);
    SID = wordsurp_encoding.SID(el);
    elid = wordsurp_encoding.el(el);

    x_surp = wordsurp_encoding.eng_full_beta{el}(13:end, :);
    y_surp = wordsurp_encoding.sp_full_beta{el}(13:end, :);

    if prctile(y(:), 90) > 0.0 && prctile(x(:), 90) > 0.0
        % correlation between beta weights
        c(el) = corr(x(:), y(:));

         %if ismember(el, [32, 15, 100, 105, 273, 228, 236, 81, 16]) % [223 81 895 133 300]
        if strcmp(SID{1}, 'EC183') && ...
                wordsurp_encoding.el(el)==22
                %wordsurp_encoding.eng_base_rsq(el)>0.05
            % [223 81 895 133 300]
         
%         if c(el) > 0.95
            figure; 
            subplot(1, 4, 1);
            imagesc(x);
            yticks(1:12);
            yticklabels(feats);
%             caxis([-2 2])
    
            subplot(1, 4, 2);
            imagesc(y);
            yticks(1:12);
            yticklabels(feats);
%             caxis([-2 2]);
            colormap(flipud(rdbu));  
            sgtitle(['Phn Feat: ' SID{1} ', ' num2str(elid)])

            subplot(1, 4, 3);
            imagesc(x_surp);
%             caxis([-2 2])
    
            subplot(1, 4, 4);
            imagesc(y_surp);
%             caxis([-2 2]);
    
            colormap(inferno);       % flipud(rdbu)
        end
    end
end

figure;
ww = min(wordsurp_encoding.eng_base_rsq, wordsurp_encoding.sp_base_rsq);
[histw, vinterval] = histwc(c, ww, 15);
bar(vinterval, histw);
ylabel('R^2 weighted density');
xlabel('phn feat weight corr');

figure;
lsidx = wordsurp_encoding.ls==1;
scatter3(ww(lsidx), c(lsidx), 1:sum(lsidx), 25, ...
    wordsurp_encoding.sp_uv_all(lsidx, 5), 'filled'); hold on;

lsidx = wordsurp_encoding.ls==2;
scatter3(ww(lsidx), c(lsidx), 1:sum(lsidx), 25, ...
    wordsurp_encoding.eng_uv_all(lsidx, 5), 'filled', 'Marker', '^');

colormap(flipud(rdbu));  
xlabel('min model R^2')
ylabel('phn feat weight corr')
set(gca, 'FontSize', 15);
caxis([-0.02 0.02])
view(2);

% Showing r^2 value by phonetic feature correlation
figure('Renderer', 'Painters');
for i = [1:2 4]
    subplot(2, 2, i);
    lsidx = wordsurp_encoding.ls==i & ~isnan(ww) & ~isnan(c);
    scatter3(ww(lsidx), c(lsidx), 1:sum(lsidx), 25, [0.7 0.7 0.7], ...
        'filled', 'MarkerFaceAlpha', 0.7);
    hold on;

    order = 2;  % Or 2 to a quadratic.
    % First the first plot.
    coefficients1 = polyfit(ww(lsidx), c(lsidx), order);
    % Make new x coordinates
    x1 = linspace(min(ww(lsidx)), max(ww(lsidx)), 50);
    y1 = polyval(coefficients1, x1);
    a=plot(x1, y1, 'k-', 'LineWidth', 2);
    uistack(a, 'top')
    
    xlabel('model R^2')
    ylabel({'feature', 'weight correlation'})
    set(gca, 'FontSize', 15);
    caxis([-0.02 0.02]);
    xlim([0 0.35]);
    xticks([0 0.3])
    view(2);
    grid off;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;



%% SCRATCH

%% Example beta weights for aud, TRF, formant ERP

% Dcons = loadDDcons('dimex', 20, 50, {'EC172'}, [], [], ...
%     dimex_details.sentdet, dimex_details.phnnames);
% TDcons = loadDDcons('timit', 20, 50, {'EC172'}, [], [], ...
%     timit_details.sentdet, timit_details.phnnames);

% el = 137; 
% SID = 'EC172';

% example electrodes in the figure
SID = 'EC172';
els = [25, 74,  136, 137]; % 136, 137, 
% SID = 'EC100';
% els = 70;
% Scons = {Dcons, TDcons};
corpusnames = {'dimex', 'timit'};
cols = [0.2 0.2 0.2; 0.6 0.6 0.6];
details = {dimex_details, timit_details};
bef=20;

modelname={'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; 
modelfeatures  = [{'onset'}; timit_details.features.names([1:3, 8, 9, 11]); ... 
    {'peakrate'; 'F1'; 'F2'; 'F3'; 'F4'}];
for el = els
    figure; 
    weights = cell(2, 2);
    for corp = 1:2
        ctr=1;
    %     for feat = {'sonorant', 'obstruent'}
        
            subplot(2, 2, 1 + (corp-1)*2);
    
    %         % previous sonorant vs. obstruent comparison
    %         Scon = Scons{corp};
    %         resp = Scon.(SID).resp(el, :, :);
    % 
    %         phnidx = ismember(Scon.phn, details{corp}.features.(feat{1}));
    %         y = squeeze(resp(1, 1:61, phnidx));
    %         x = (-bef:1:(aft-10))/100;
    %         nanmean = @(x) mean(x, 'omitnan')./max(mean(x, 'omitnan'));
    % 
    %         shadedErrorBar(x,y',{@nanmean,@nansem},{'-', ...
    %             'color', cols(ctr, :),          ...
    %             'markerfacecolor',cols(ctr, :), ...
    %             'Linewidth', 2, 'DisplayName', corpusnames{corp}}); 
    %         hold on;   
    %         yticks([]);
    
            if ctr==1
                % read in strf weights
                [weights{1, corp}] = getTRFweights(SID, el, corpusnames{corp}, ...
                    {'onset_aud'}, datapath);
    
                ax = subplot(2, 2, 1 + (corp-1)*2);
                imagesc(weights{1, corp}(:, 1:40));
                if corp==1
                    colormap(ax, flipud(blues));
                else
                    colormap(ax, flipud(reds));
                end
                %colormap(inferno)
    
                yticks([1 80]);
                yticklabels({'0', '8'});
                xlim([0.5 40])
                xticks([1 40]);
                xticklabels({'0', '-0.4'});
                xlabel('Time (s)');
                set(gca, 'FontSize', 13, 'Ydir', 'normal');
    
                % read in model trf weights   
                [weights{2, corp}] = getTRFweights(SID, el, corpusnames{corp}, ...
                    modelname, datapath);
                
                ax = subplot(2, 2, corp*2);
                imagesc(weights{2, corp}(:, 1:40));
                if corp==1
                    colormap(ax, flipud(blues));
                else
                    colormap(ax, flipud(reds));
                end
    
                yticks(1:length(modelfeatures));
                yticklabels(modelfeatures);
                xlim([0.5 40])
                xticks([1 40]);
                xticklabels({'0', '-0.4'});
                xlabel('Time (s)');
                set(gca, 'FontSize', 13);
                clear strf
            end
            ctr=ctr+1;
        %end
        
        subplot(2, 2, 1 + (corp-1)*2);
        yline(0, 'LineWidth', 1.5);
        xline(0, 'LineWidth', 1.5);
        xlabel('Time (s)');
    %     xticks(-0.2:0.2:0.4)
        set(gca, 'FontSize', 13);
        ylabel('HFA (z)');   
        box off;
        
    end
    sgtitle(corr(weights{2, 1}(:), weights{2, 2}(:)))
end


% addpath(genpath('util'))
% TDvow = loadDD('timit', bef, aft,{SID});
% 
% addpath(genpath('util'))
% Dvow = loadDD('dimex', bef, aft,{SID});
% 
% plotFormantErp(Dvow, SID, el, []);                                   
% plotFormantErp(TDvow, SID, el, []);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *Dcons *wrd*;

%% functions
function [weights] = getTRFweights(SID, el, corpus, modelname, datapath)    
    [strf] = loadMultModelStrf(SID, modelname, corpus, datapath, 1);  
    weights = strf{1}.meanStrf(:, :, el);
end



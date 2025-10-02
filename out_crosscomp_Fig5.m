% Ilina Bhaya-Grossman
% 01.08.2022
out_crosscomp_startup;

[~,~, bSIDs, ~] = getSIDinfo();
SIDs = bSIDs;
tps = 50:55;

% load word structures
if ~exist('Dwrd', 'var')
    load("data/Figure3/Figure3_DIMEXWrd.mat");
end

if ~exist('TDwrd', 'var')
    load("data/Figure3/Figure3_TIMITWrd.mat");
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* *modelnames*;

%% Load in sentence responses

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
        % type is whether selected as speech responsive for timit / dimex
        SID = sent_encoding.SID{el};
        ls = sent_encoding.ls(el);
        if ismember(ls, [1, 2, 4]) && isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            type(:, el) = [ismember(sent_encoding.el(el), dimex_elecs.allidx.(SID)) ...
                    ismember(sent_encoding.el(el), timit_elecs.allidx.(SID))];
            clear SID
        elseif ismember(ls, 3) && isfield(asccd_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            type(:, el) = [ismember(sent_encoding.el(el), asccd_elecs.allidx.(SID)) ...
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
    betaInfo* *encoding* allidx fthresh Dcons *wrd* *modelnames*;

%% Percentage of overlap in electrodes in English vs. Spanish

dimex_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
timit_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');

overlap = nan(length(SIDs), 1);
ctr=1;
for s = SIDs
    SID = s{1};
    if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
        dimex = dimex_elecs.allidx.(SID);
        timit = timit_elecs.allidx.(SID);
        overlap(ctr) = length(intersect(dimex, timit))/length(union(dimex, timit));
        disp([SID ': ' num2str(overlap(ctr)) ' overlap, total: ' num2str(length(union(dimex, timit)))]);
        % number of electrodes in each language
        disp(['Dimex: ' num2str(length(dimex)) ' Timit: ' num2str(length(timit))]);
    end
    ctr=ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* *modelnames*;

%% A - Single scatter of speech responsive electrodes

imgall = load_allimgdata;   
cols = getColorsCrossComp(6);
titles = {'English', 'Spanish'};
anat_counts = cell(2, 1);

types = [0, 1];
for native = types
    ctr = 1;
    figure;
    set(gcf,'Color','w');
    for h = {'lh', 'rh'}
    
        hemi = h{1};
        if strcmp(hemi, 'lh')
            x_add = -20;
        else
            x_add = 10;
        end
    
        numsid = 0;      
    
        if strcmp(hemi, 'lh')
            cortex = imgall.(SIDs{1}).img_mni.cortex;
        else
            cortex = imgall.(SIDs{4}).img_mni.cortex;
        end
    
        subplot(1, 2, ctr)
        hold on
    
        PlotBrainSurface(cortex, hemi,'lateral');
        alpha 0.9
        %light("Style","infinite","Position",[100 100 0]);
        
        % find density map of native speech, lateral side
        native_xyz = [];
        native_hga = [];
        anat_xyz = [];
    
        for si = SIDs
            sid = si{1};
            if strcmpi(imgall.(sid).hemi,hemi)
                subject = imgall.(sid);
                
                if isfield(subject.img_mni, 'elecmatrix')

                    if isfield(subject.img_native, 'anatomy')
                        anat = subject.img_native.anatomy(:, 4);
                    elseif isfield(subject.imgNative, 'anatomy')
                        anat = subject.imgNative.anatomy(:, 4);
                    else
                        disp(['using mni anatomy labels... for ' sid])
                        anat = subject.img_mni.anatomy(:, 4);
                    end

                    elecmatrix = subject.img_mni.elecmatrix;             
                    ch_sid = sent_encoding.el(strcmp(sent_encoding.SID, sid));
                    ls = sent_encoding.ls(find(strcmp(sent_encoding.SID, sid), 1, 'first'));
            
                    if native % spanish for bilinguals
                        if ismember(ls, [1, 2])
                            ch_sel = sent_encoding.el(sent_encoding.type(:, ls)==1);
                        else % bilingual case and mandarin case (where first column is mandarin)
                            ch_sel = sent_encoding.el(sent_encoding.type(:, 1)==1);
                        end
        
                        % for weighting
                        if ismember(ls, [1, 3, 4])
                            native_hga = [native_hga; 
                                sent_encoding.maxresp(intersect(ch_sel, ch_sid), 1)];
                        elseif ls==2
                            native_hga = [native_hga; 
                                sent_encoding.maxresp(intersect(ch_sel, ch_sid), 2)];
                        end
                    else
                        if ismember(ls, [1, 2])
                            ch_sel = sent_encoding.el(sent_encoding.type(:, mod(ls, 2)+1)==1);
                        else % bilingual case and mandarin case
                            ch_sel = sent_encoding.el(sent_encoding.type(:, 2)==1);
                        end
        
                        % for weighting
                        if ismember(ls, [1, 3, 4])
                            native_hga = [native_hga; ...
                                sent_encoding.maxresp(intersect(ch_sel, ch_sid), 2)];
                        elseif ismember(ls, 2)
                            native_hga = [native_hga; ...
                                sent_encoding.maxresp(intersect(ch_sel, ch_sid), 1)];
                        end
                    end
        
                    native_xyz = [native_xyz; elecmatrix(intersect(ch_sel, ch_sid),:)]; 
                    anat_xyz = [anat_xyz; anat(intersect(ch_sel, ch_sid))];
                    numsid=numsid+1;
                end
            end
        end

        % anatomy overview
        [counts, ~, ~, labels] = crosstab(anat_xyz);
        [sorted, idx] = sort(counts, 'descend');
        disp('Anatomy overview');
        % maximum 4 anatomy
        disp(['Max: ' labels{idx(1)} ' ' num2str(sorted(1))]);
        disp(['2nd: ' labels{idx(2)} ' ' num2str(sorted(2))]);
        disp(['3rd: ' labels{idx(3)} ' ' num2str(sorted(3))]);
        disp(['4th: ' labels{idx(4)} ' ' num2str(sorted(4))]);
        disp('--------------------------------------------------------');

        % if anat_counts{native+1} is already a table, add to existing, otherwise make a new table
        if istable(anat_counts{native+1})
            % find equivalent area in table and add count
            for i = 1:length(labels)
                area_idx = find(strcmp(anat_counts{native+1}.area, labels(i)));
                if ~isempty(area_idx)
                    anat_counts{native+1}.count(area_idx) = ...
                        anat_counts{native+1}.count(area_idx) + counts(i);
                else
                    anat_counts{native+1} = [anat_counts{native+1}; ...
                        table(counts(i), labels(i), 'VariableNames', ...
                        {'count', 'area'})];
                end
            end
        else
            anat_counts{native+1} = table(sorted, labels(idx), ...
                'VariableNames', {'count', 'area'});
        end
    
        % number of subjects
        disp(['Is native: ' num2str(native)]);
        disp(['Number of subjects: '  num2str(numsid)]);
        disp(['Number of electrodes: '  num2str(length(native_hga))]);
        disp('--------------------------------------------------------');
        if ~native
            color = cols(1, :);
        else
            color = cols(2, :);
        end
    
        scatter3(native_xyz(:, 1)+x_add, native_xyz(:, 2), native_xyz(:, 3), ...
            round(native_hga*20), color, 'filled', 'MarkerFaceAlpha', 0.7);
        ctr=ctr+1;
    end
    sgtitle(titles{native+1});
end

% Check the maximum size for legend

% Create horizontal bar charts for speech-responsive areas
figure;
subplot(1, 2, 1);
% remove areas with count < 5 and sort by count
anat_counts{1} = anat_counts{1}(anat_counts{1}.count >= 5, :);
anat_counts{1} = sortrows(anat_counts{1}, 'count', 'descend');
barh(anat_counts{1}.area, anat_counts{1}.count, 'FaceColor', ...
    cols(1, :), 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Red color for non-native
title(titles{1});
xlabel('Count');
ylabel('Anatomical Area');
set(gca, 'FontSize', 14);
box off;
xlim([0, 310]);

subplot(1, 2, 2);
% remove areas with count < 5 and sort by count
anat_counts{2} = anat_counts{2}(anat_counts{2}.count >= 5, :);
anat_counts{2} = sortrows(anat_counts{2}, 'count', 'descend');
barh(anat_counts{2}.area, anat_counts{2}.count, 'FaceColor', ...
    cols(2, :), 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Blue color for native
title(titles{2});
xlabel('Count');
ylabel('Anatomical Area');
set(gca, 'FontSize', 14);
sgtitle('Speech-responsive areas');
box off;
xlim([0, 310]);


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* *modelnames*;

%% B/C - Correlating (s)TRF weights across corpora

corpora = {{'timit', 'dimex'}};
%modelnames = {'onset_aud'};
modelnames = {'onset_phnfeatConsOnset_maxDtL_formantMedOnset'};

% consonant feature overlap for mandarin and english
timit_feats = timit_details.features.names([1:3, 8, 9, 11]);
% shared_feats = intersect(timit_feats, asccd_feats);
% timit_fidx = [1; 1+cellfun(@(x) find(strcmp(timit_feats, x)), shared_feats); ...
%     (12-4:12)'];

corrstrf = nan(height(sent_encoding), 1);
corrstrfperm = nan(height(sent_encoding), 10);
corrpval = corrstrf;
maxfreq = nan(height(sent_encoding), 2);
maxrsq = nan(height(sent_encoding), 1);
minrsq = nan(height(sent_encoding), 1);
thresh = 0.01;

% window is 5 time points around the max time point for correlation
windsz = 3;

colors = flipud(brewermap(20, 'Spectral'));
% cols = [colors(3, :); colors(end-2, :); 
%     colors(round(size(colors, 1)/2), :)];

% all windowed TRF weights
alleng = nan(height(sent_encoding), 80);
allsp = nan(height(sent_encoding), 80); % maximum feat length

SIDs =  unique(sent_encoding.SID)';

% remove subjects where information not completely loaded
SIDs(ismember(SIDs, {'EC282', 'EC296'})) = [];

for s = SIDs
    
    % load both corpus STRFs
    SID = s{1};
    if strcmp(modelnames{1}, 'onset_aud')
        meanStrf = nan(2, 81, 61, 256);
    elseif startsWith(SID, 'EC')
        meanStrf = nan(2, 12, 61, 256);
    else
        meanStrf = nan(2, 9, 61, 256);
    end

    if startsWith(SID, 'EC') % all EC patients are English / Spanish
        lang_pair = 1;
        if ~strcmp(modelnames{1}, 'onset_aud')
            fidx = {1:12, 1:12};
        else
            fidx = {1:81, 1:81};
        end
    elseif startsWith(SID, 'HS') % all HS patients are English / Mandarin
        lang_pair = 2;
        % match features for Mandarin and English
        if ~strcmp(modelnames{1}, 'onset_aud')
            fidx = {timit_fidx, asccd_fidx};
        else
            fidx = {1:81, 1:81};
        end
    end

    for c = 1:2
        corpus = corpora{lang_pair}{c};
        corpusStrf=loadMultModelStrf(SID, modelnames, corpus, datapath, 1, ...
            'v5');
        if ~isempty(corpusStrf{1})
            minel = size(corpusStrf{1}.meanStrf, 3);
            meanStrf(c, :, :, 1:minel) =   corpusStrf{1}.meanStrf(fidx{c}, :, :);
        end
    end
    clear fidx

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
            
            % save out weights
            alleng(e, 1:length(x)) = x;
            allsp(e, 1:length(z)) = z;

            % calculate correlation and max frequency
            [corrstrf(e), corrpval(e)] = corr(x', z, 'type', 'Pearson');
            [~, maxfreq(e, :)] = max([x', z]);
            %freq(e, :, :) = [x'; z];

            % calculate the permutation correlation
            for j = 1:10
                rng(1);
                [corrstrfperm(e, j), ~] = corr(x', z(randperm(length(z))), ...
                    'type', 'Pearson');
            end

            % save out max R^2 and min R^2
            maxrsq(e) = max(sent_encoding.en_base_rsq(e), ...
                sent_encoding.sp_base_rsq(e));
            minrsq(e) = min(sent_encoding.en_base_rsq(e), ...
                sent_encoding.sp_base_rsq(e));
        
            debug = 0;
            % Example electrodes from above (EC100, 22 / 150)
            if debug && minrsq(e)>0.1
                % (corrstrf(e)>0.9) && minrsq(e)>0.15
                % ismember(e, [184, 120, 133, 136, 137, 473, 864, 1125, 390])
                % (corrstrf(e)>0.3 && corrstrf(e)<0.6 && minrsq>0.15)

                figure;
                subplot(2, 6, 1);
                plot(x, 1:length(x), 'LineWidth', 2, 'Color',cols(2, :));
                set(gca, 'XDir', 'reverse');
                ylim([1 length(x)+1]);
                hold on;                

                % Plot the STRF beta weights for the first language
                subplot(2, 6, [2, 3]);
                imagesc(squeeze(meanStrf(1, :, :, el)));
                set(gca, 'YDir', 'normal');
                xline(wind(1));
                xline(wind(end));
                title('timit');
                
                if strcmp(modelnames{1}, 'onset_aud')
                    yticks([1 80]);
                    yticklabels({'0.01', '8'});
                    ylabel('frequency (kHz)')
                end

                subplot(2, 6, 4);
                plot(z, 1:length(x), 'LineWidth', 2, ...
                    'LineStyle', '-', 'Color', cols(1, :));
                set(gca, 'XDir', 'reverse');
                hold on; 
                ylim([1 length(x)+1])
        
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
                plotStitchedSentence(e, sent_encoding, 150, 100, 0.5, 1);
                title(['Language exp: ' num2str(sent_encoding.ls(e))]);
            end  
        end
    end
end

% show pie for correlation threshold
idx = minrsq>0.05;
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
legend({'r<=0.5', 'r>0.5'});

% show imagesc of TRF weights comparing the two corpora
% find nan rows for alleng and allsp
if ~strcmp(modelnames{1}, 'onset_aud')
    alleng = alleng(:, 1:12);
    allsp = allsp(:, 1:12);
end
nanrows = any(isnan(alleng), 2) | any(isnan(allsp), 2) | minrsq<0.05;
alleng(nanrows, :) = [];
allsp(nanrows, :) = [];

% absolute value and zscore the weights
alleng = normalize(alleng, 1, 'center');
allsp = normalize(allsp, 1, 'center');

% sort by maximum feature, 
[~, maxidx] = maxk(alleng, 3, 2);
[~, idx] = sortrows(maxidx);
alleng = alleng(idx, :);
allsp = allsp(idx, :);

figure;
subplot(2, 1, 1);
imagesc(alleng');
colormap(flipud(prgn));
brighten(0.2);
clim([-1 1]);
box off;

[corrval, pval] = corr(allsp(:), alleng(:));
disp(['Pearson r(' num2str(length(allsp(:))) '): ' num2str(corrval) ', p = ' num2str(pval)]);

% add xlines where max feature changes (1 to 2, 7 to 8, 8 to 9)
sidx = maxidx(idx, 1);
for i = 1:length(sidx)
    if i>1 && sidx(i)~=sidx(i-1) && ismember(sidx(i-1), [1, 7, 8])
        xline(i-1, 'LineWidth', 2);
    end
end
clear sidx i

% add the feature labels
yticks(1:12);
yticklabels([{'onset'} timit_feats' {'peakrate', 'F1', 'F2', 'F3', 'F4'}]);

subplot(2, 1, 2);
imagesc(allsp');
brighten(0.2);
clim([-1 1]);
box off;

% add xlines where max feature changes (1 to 2, 7 to 8, 8 to 9)
sidx = maxidx(idx, 1);
for i = 1:length(sidx)
    if i>1 && sidx(i)~=sidx(i-1) && ismember(sidx(i-1), [1, 7, 8])
        xline(i-1, 'LineWidth', 2);
    end
end
clear sidx i

% add the feature labels
yticks(1:12);
yticklabels([{'onset'} timit_feats' {'peakrate', 'F1', 'F2', 'F3', 'F4'}]);

% show histogram of correlations
figure;
% lower threshold in Fig. 1, fewer participants
idx = minrsq>0.05;
hbins = -1:0.05:1;
% permutation distribution
corrstrfperm(~idx, :) = NaN;
corrstrfperm(isnan(corrstrfperm)) = [];

histogram(corrstrfperm(:), hbins, 'EdgeColor', 'none', ...
    'FaceColor', [0.7, 0.7, 0.7], 'FaceAlpha', 0.3, ...
    'Normalization','pdf');
hold on;
% make a kdensity plot
[f, xi] = ksdensity(corrstrfperm(:));
plot(xi, f, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 3, ...
    'LineStyle', '-');
yticks(0:1:2);

yyaxis right; 
histogram(corrstrf(idx), hbins, 'EdgeColor', 'none', ...
    'FaceColor', [0.7, 0.1, 0.7], 'FaceAlpha', 0.3, ...
    'Normalization', 'pdf');
hold on;
[f, xi] = ksdensity(corrstrf(idx)); 
plot(xi, f, 'LineWidth', 2, 'Color', [0.7, 0.1, 0.7], 'LineWidth', 3, ...
    'LineStyle', '-');

stathresh = (0.05)*100;
corrstrfperm(isnan(corrstrfperm)) = [];
xline(prctile(corrstrfperm(:), 100-stathresh), 'LineWidth', ...
    2, 'Color', 'k');
xlim([-1 1]);
yticks(0:2:6);
ax = gca();
ax.YAxis(2).TickLabelColor = [0.7, 0.1, 0.7]; 
ax.YAxis(2).Color = [0.7, 0.1, 0.7]; 
% ax.YAxis(2).

% yticks([]);
yyaxis left;
ylabel('Probability Density');
box off;
set(gca, 'FontSize', 13);

disp(['Total number of electrodes: ' num2str(sum(idx))])
disp(['Number of electrodes with higher than 95% for permuted distrib: ' ...
    num2str(sum(corrstrf(idx)>prctile(corrstrfperm(:), 100-stathresh)))]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq *modelnames*;

%% Load in unique variance

[wordsurp_encoding] = loadUniqueVarTbl(modelnames_timit, modelnames_dimex, [bSIDs eSIDs, sSIDs]);
wordsurp_details.featureOrd ...
    = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', ...
    'surp', 'wordO', 'wordL', 'wordF', 'bisurp', 'trisurp', 'pitch', 'env'};
imgall = load_allimgdata;
wordsurp_encoding.hemi = cellfun(@(x) imgall.(x).hemi, ...
    wordsurp_encoding.SID, 'UniformOutput', false);
wordsurp_details.models_dimex = modelnames_dimex;
wordsurp_details.models_timit = modelnames_timit;

% load in p-values from permutation testing
sp_wordsurp_pval = nan(1, height(wordsurp_encoding));
en_wordsurp_pval = nan(1, height(wordsurp_encoding));
prefix = 'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_';
permodel = {[ prefix 'engSurpNoOnsBin_wordFreqLog'], ...
    [prefix 'spSurpNoOnsBin_wordFreqLog']};
for s = unique(wordsurp_encoding.SID)'
    SID = s{1};
    idx = strcmp(wordsurp_encoding.SID, SID);
    elidx = wordsurp_encoding.el(idx);

    for c = 1:2 %{'timit', 'dimex'}
        pvalpath=fullfile(datapath, 'permTest_wordSurp', SID); % c{1}, 
        cmod=dir(fullfile(pvalpath, '*_zX*_*mat')); 
        permidx = find(contains({cmod.name}, permodel{c}));
        permfname=cmod(permidx).name;

        % load pvalues from permutation testing
        pvals = load(fullfile(pvalpath, permfname), 'pval');
        pvals = pvals.pval;
        
        if c==1 % strcmp(c{1}, 'timit')
            en_wordsurp_pval(idx) = pvals(elidx);
        else
            sp_wordsurp_pval(idx) = pvals(elidx);
        end
    end
end
wordsurp_encoding.sp_wordsurp_pval = sp_wordsurp_pval';
wordsurp_encoding.eng_wordsurp_pval = en_wordsurp_pval';

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd *elecs;

%% D - English vs. Spanish: word / sequence surprisal unique variance (box-plots)

% Compare unique variance between languages
rsq_thresh = 0.05;
idx = (wordsurp_encoding.eng_base_rsq>rsq_thresh& wordsurp_encoding.sp_base_rsq>rsq_thresh);
wordsurp_encoding_thresh = wordsurp_encoding(idx, :);

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 

conditionLabels = {'English', 'Spanish'};
thresh = 0.001;

% label = 'word boundary';
% label = 'phonetic feature';
% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% peakrate, phonetic features, word+surp
% 'peakrate', 'formant', 'consonant', 'word+surp',
features = { 'wordO', 'wordF', 'wordL' 'surp'};
titles = { 'onset', 'frequency', 'length' , 'surprisal'};
fields = {'eng_uv_all', 'sp_uv_all'}; 
cols = getColorsCrossComp(6);

[proficiency] = getBilingProf();
% profield = {'eng_prof', 'span_prof'};

% two primary proficiency profiles (high eng, low span and vice versa)
hsle = proficiency.SID(proficiency.span_prof>4 & proficiency.eng_prof<5);
lshe = proficiency.SID(proficiency.span_prof<5 & proficiency.eng_prof>4);

figure('Renderer', 'Painters');
ctr = 1;
for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);
    bidx = wordsurp_encoding_thresh.ls==4; %& ismember(wordsurp_encoding.SID, lshe);

    eng_uv = wordsurp_encoding_thresh.(fields{1})(bidx, index);
    sp_uv = wordsurp_encoding_thresh.(fields{2})(bidx, index);
   
    % Combining native and unfamiliar subjects  
    % subplot(1, length(titles), ctr)
    h=boxchart(ones(sum(eng_uv>thresh), 1)+(ctr-1)*1.5, eng_uv(eng_uv>thresh), ...
        'BoxFaceColor', cols(1, :), 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.3); hold on;
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    h=boxchart(ones(sum(sp_uv>thresh), 1)*1.5+(ctr-1)*1.5, sp_uv(sp_uv>thresh), ...
        'BoxFaceColor', cols(2, :), 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.3);
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([eng_uv(eng_uv>thresh); sp_uv(sp_uv>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    xticks([1 2]);
    xticklabels(conditionLabels);
    xlim([0.5 2.5])
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    tbl.rsq = [eng_uv; sp_uv];
    tbl.sid = [wordsurp_encoding_thresh.SID(bidx); wordsurp_encoding_thresh.SID(bidx)];
    tbl.elec = [wordsurp_encoding_thresh.el(bidx); wordsurp_encoding_thresh.el(bidx)];
    tbl.hemi = [wordsurp_encoding_thresh.hemi(bidx); wordsurp_encoding_thresh.hemi(bidx)];
    tbl.lang = [ones(length(eng_uv), 1); ones(length(sp_uv), 1)*2]; % english and spanish

    lme2 = fitlme(tbl(tbl.rsq>0, :),'rsq~lang+(1|hemi)+(1|sid)+(1|elec:sid)');

    disp(titles{ctr})
    disp(lme2)
    p = lme2.Coefficients.pValue(2);

    disp(['LME L1 language p-value = ' num2str(p)])
    if isempty(getSigStr(p, 2))
        text((ctr-1)*1.5 + 1, maxy(2)-0.05, 'n.s.', 'FontSize', 15);
    else
        text((ctr-1)*1.5 + 1, maxy(2)-0.05, getSigStr(p, 1), 'FontSize', 15);
    end
    %title(titles{ctr});
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2 (log)');

    ctr=ctr+1;
end

xlim([0.5 1.5*length(features)+0.5]);
xticks((1:length(features))*1.5-0.25);
xticklabels(titles);
legend({'English', 'Spanish'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq;
%% E - Monolingual vs. Bilingual: word / sequence surprisal unique variance (box-plots)

rsq_thresh = 0.05;
idx = (wordsurp_encoding.eng_base_rsq>rsq_thresh& wordsurp_encoding.sp_base_rsq>rsq_thresh);
wordsurp_encoding_thresh = wordsurp_encoding(idx, :);

% field = {'eng_rsq_surprisal', 'sp_rsq_surprisal'}; 
% field = {'eng_uv_phnfeat', 'sp_uv_phnfeat'}; 
fieldname = {'monolingual', 'bilingual'};
conditionLabels = {'monolingual', 'bilingual'};
thresh = 0.001 ;

% label = 'word boundary';
% label = 'phonetic feature';

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
% peakrate, phonetic features, word+surp
% features = {'word+surp', 'word', 'surp'};
% titles = {'word+surp', 'word', 'surp'};

features = {'wordO', 'wordF', 'wordL', 'surp'};
titles = {'onset', 'frequency', 'length', 'surp'};
% titles = {'peakrate', 'formant', 'consonant', 'word+surp'};
% titles = {'peakrate', 'phonetic features', 'word+surp', 'word', 'surp'};
fields = {'eng_uv_all', 'sp_uv_all'}; 

% Combine the consonant and vowel features into one
subplts = {1, 2, 3, 4};

ctr = 1;
mono = cell(length(features), 1);
sid_mono = cell(length(features), 1);
bil = cell(length(features), 1);
sid_bil = cell(length(features), 1);
sid = cell(length(features), 1);
el = cell(length(features), 1);
hemi = cell(length(features), 1);

[proficiency] = getBilingProf();
profield = {'eng_prof', 'span_prof'};

for feat = features
    index = ismember(wordsurp_details.featureOrd, feat);

    for f = 1:2 % Aggregate over English and Spanish fields
        field = fields{f};    
        y = wordsurp_encoding_thresh.(field)(:, index);   
        sid_tmp = wordsurp_encoding_thresh.SID; 
        
        % Removing bilinguals so the comparison is more balanced
        if f == 1 
            mono{ctr} = [mono{ctr}; y(ismember(wordsurp_encoding_thresh.ls, 2))];  
            sid{ctr} = [sid{ctr}; sid_tmp(ismember(wordsurp_encoding_thresh.ls, 2))];
            el{ctr} = [el{ctr}; wordsurp_encoding_thresh.el(wordsurp_encoding_thresh.ls==2)];
            hemi{ctr} = [hemi{ctr}; wordsurp_encoding_thresh.hemi(wordsurp_encoding_thresh.ls==2)];
        else
            mono{ctr} = [mono{ctr}; y(ismember(wordsurp_encoding_thresh.ls, 1))];  
            sid{ctr} = [sid{ctr}; sid_tmp(ismember(wordsurp_encoding_thresh.ls, 1))];
            el{ctr} = [el{ctr}; wordsurp_encoding_thresh.el(wordsurp_encoding_thresh.ls==1)];
            hemi{ctr} = [hemi{ctr}; wordsurp_encoding_thresh.hemi(wordsurp_encoding_thresh.ls==1)];
        end

        % subset to bilinguals with full proficiency    
        profSIDs = [proficiency.SID(proficiency.(profield{f})==5)];
        fullprof_biling = wordsurp_encoding_thresh.ls==4 & ...
            ismember(wordsurp_encoding_thresh.SID, profSIDs);
        bil{ctr} = [bil{ctr}; y(fullprof_biling)];
        sid{ctr} = [sid{ctr}; sid_tmp(fullprof_biling)];
        el{ctr} = [el{ctr}; wordsurp_encoding_thresh.el(fullprof_biling)];
        hemi{ctr} = [hemi{ctr}; wordsurp_encoding_thresh.hemi(fullprof_biling)];

    end
    ctr = ctr + 1;
end

% split by sid
ynat_sid = cell(length(subplts), 2);
ynon_sid = cell(length(subplts), 2);

figure('Renderer', 'Painters');
for s = 1:length(subplts)
    
    % Combining native and unfamiliar subjects
    % subplot(1, length(subplts), s)
    rep = length(subplts{s});
    
    % Native boxplot
    y_nat = arrayfun(@(x)  mono{x}, [subplts{s}], 'UniformOutput', false);
    y_nat = cat(1, y_nat{:});

    h=boxchart(ones(sum(y_nat>thresh), 1)+(s-1)*1.5, y_nat(y_nat>thresh), ...
        'BoxFaceColor', [0.5 0.5 0.5], 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.3); hold on;
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Unfamiliar boxplot
    y_non = arrayfun(@(x) bil{x}, [subplts{s}], 'UniformOutput', false);
    y_non = cat(1, y_non{:});

    h=boxchart(ones(sum(y_non>thresh), 1)*1.5 +(s-1)*1.5, y_non(y_non>thresh), ...
        'BoxFaceColor', [0.1 0.1 0.5], 'MarkerColor', 'k', 'Notch','on', 'BoxWidth', 0.3);
    h.JitterOutliers = 'on';
    h.MarkerStyle = '.';
    h.MarkerColor = 'k';

    % Formatting
    ylim([0 prctile([y_nat(y_nat>thresh); y_non(y_non>thresh)], 99)+0.05]);
    yticks(0:0.05:0.1)
    xticks([1 2]);
    xticklabels(conditionLabels);
    xlim([0.5 2.5])
    maxy = ylim();
    set(gca, 'YScale', 'log')

    % Statistical testing with linear mixed effect model
    tbl=table();
    tbl.rsq = [y_nat; y_non];
    
    elecs = el{1};
    sids = sid{1};
    hemis = hemi{1};

    tbl.sid = sids;
    tbl.elec = elecs;
    tbl.hemi = hemis;
    tbl.native = [ones(length(y_nat), 1); ones(length(y_non), 1)*2];

    clear elecs sids

    lme2 = fitlme(tbl(tbl.rsq>0, :),'rsq~native+(1|hemi)+(1|sid)+(1|elec:sid)');

    disp(titles{s})
    disp(lme2.Coefficients(2, :))
    p = lme2.Coefficients.pValue(2);

    disp(['LME mono / bilingual p-value = ' num2str(p)])
%     line([1, 2], [maxy(2)-0.1 maxy(2)-0.1], 'Color', 'k');
    if isempty(getSigStr(p, 2))
        text((s-1)*1.5 + 1, maxy(2)-0.05, 'n.s.', 'FontSize', 15);
    else
        text((s-1)*1.5 + 1, maxy(2)-0.05, getSigStr(p, 1), 'FontSize', 15);
    end
    %title(titles{s});
    set(gca, 'FontSize', 13);
    ylabel('\Delta R^2 (log)');
end
legend({'monolingual', 'bilingual'});
xlim([0.5 1.5*length(subplts)+0.5]);
xticks((1:length(subplts))*1.5-0.25);
xticklabels(titles);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq;

%% F - Word boundary ERPs

f = figure; 

% Bilingual
% EC163 - el55
SIDs = {'EC260'};% 'EC100', 'EC100', 'EC100','EC100'};
els = 221;%,135, 22, 70, 71, 150]; %236, 205 for EC260, 178, 175

Dwrd = loadDwrdMulti('dimex',  20, 50, SIDs, dimex_details);
TDwrd = loadDwrdMulti('timit', 20, 50, SIDs, timit_details);
xdata = -0.2:0.01:0.5;
uv_thresh = 0.001;

plotSingleTrial = 0;
numel = length(els);
Swrds = {Dwrd, TDwrd};
for s = 1:2
    Swrd = Swrds{s};

    for ctr = 1:length(els)
        SID = SIDs{ctr};    
    
        dummy = struct();
        nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
            %| (Swrd.syll<2 & ~isnan(Swrd.syll));
        dummy.(SID).resp = cat(3, Swrd.(SID){~nanidx});
        dummy.wordOns = Swrd.wordOns (~nanidx);
        dummy.syllOns = ones(sum(~nanidx), 1);

        ls = wordsurp_encoding.ls(find(strcmp(wordsurp_encoding.SID, SIDs{1}), 1));
        if s == ls   
            cols = [0 0 0; 0.1, 0.1 0.9]; 
        elseif ls == 4 % bilingual case
            cols = getColorsCrossComp(6);
            cols = [0 0 0; cols(s, :)];
        else
            cols = [0 0 0; 0.9, 0.1 0.1]; 
        end

        subplot(2, numel, ctr+(s-1)*numel)
        addpath(genpath('shadederror'))
        plotWordErp(dummy, SID, els(ctr), ...
            [], f, cols, 1, 0.2, 1); hold on;
        ylabel('HFA (z)');
        set(gca, 'FontSize', 13);

        [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
            dummy.(SID).resp(els(ctr), :, :), dummy.wordOns+1, [1, 2]);
        pval_corrected = 0.05/size(dummy.(SID).resp, 2);
        fthresh = finv(1-pval_corrected, df1, df2);  
    
        scatter(xdata(fvals>fthresh), 0.1*ones(1, sum(fvals>fthresh)), 45, ...
            fvals(fvals>fthresh), 'filled', 'HandleVisibility', 'off');
        cm = colormap("gray");
        colormap(flipud(cm(1:200, :)))
        % ylim([0 1.3]);
        ylim([-0.5 0.8]);
        xlim([-0.2, 0.3]);
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
%
% initialize design electrode structure
fieldnames = {'Spanish', 'English'};
fields = {'sp_uv_all', 'eng_uv_all', '', 'sp_uv_all'}; 
% uv feature order
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant', 'surp', 'word'
%
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
        %binedges = [-1 0.000:0.005:0.01 0.015:0.015:0.1];
        binedges = [-1 uv_thresh:0.01:0.1];
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
        elseif ls == 4
            cls = getColorsCrossComp(6);
            cls = [linspace(cls(lang, 1), 1, 10); ...
                linspace(cls(lang, 2), 1, 10); ...
                linspace(cls(lang, 3), 1, 10)]';
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
        % % add a pie
        % axes('Position',[.6 .15 .3 .3])
        % p = pie([sum(nh.cond>1), sum(nh.cond==1)], [1 1]); 
        % p(1).FaceColor = [desel.cols(5, :)];
        % p(1).EdgeColor = 'none';
        % p(3).FaceColor = [0.6 0.6 0.6];
        % p(3).EdgeColor = 'none';
        % p(2).FontWeight = 'bold';
        % p(2).Color = 'w';
        % p(2).FontSize = 13;
        % p(4).FontWeight = 'bold';
        % p(4).Color = 'w';
        % p(4).FontSize = 13;
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;
%% G - Bilingual uv as a scatter

labels = {'word+surp'};

% Create a figure for the scatter plots
figure('Position', [200, 100, 900, 900]);

% Set a UV threshold for filtering
uv_thresh = 0.001;

% Initialize subplot counter and axis array
ctr = 1;
ax = nan(3, 1);

% Iterate over the labels
for label = labels
    x_all = [];
    y_all = [];
    
    % Find the index of the current label in the featureOrd array
    index = find(ismember(wordsurp_details.featureOrd, label));
    ls = 4;
    
    % Create a subplot for the current label and ls value
    ax(ctr) = subplot(length(labels), 1, ctr);
    
    % Get the corresponding unique variance values for English and Spanish
    x = wordsurp_encoding.sp_uv_all(wordsurp_encoding.ls==ls, index);
    y = wordsurp_encoding.eng_uv_all(wordsurp_encoding.ls==ls, index);
    
    % Get the subject IDs for the current ls value
    sid = cellfun(@(x) str2double(x(3:end)), ...
        wordsurp_encoding.SID(wordsurp_encoding.ls==ls));
    
    % Remove data points that do not meet the UV threshold or contain NaN values
    neg = all([x,y]<uv_thresh, 2) | isnan(x) | isnan(y);
    x(neg) = [];
    y(neg) = [];
    sid(neg) = [];

    % Perform permutation testing to compute correlation coefficient and p-value
    maxlim = prctile([x; y], 100);
    minlim = prctile([x; y], 0);  

    % Plot the scatter plot
    colors = x-y;
    x_all = [x_all; x];
    y_all = [y_all; y];
    scatter(x, y, 25, 'k', 'filled', ... % colors
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', 0.8, ...
                'LineWidth', 0.25); hold on;

    % colorbar should be blue to red going through white (1, 1, 1)
    % first create blue to white
    cols = [linspace(63/256, 1, 50); linspace(211/256, 1, 50); linspace(242/256, 1, 50)]';
   
    % then add white to red
    cols = [cols; [linspace(1, 43/256, 50); linspace(1, 57/256, 50); ...
        linspace(1, 144/256, 50)]'];

    colormap(cols);
    clim([-0.05 0.05]);
    view(2);
    
    % Set the x and y axis limits and add labels
    xlim([minlim maxlim]);
    ylim([minlim maxlim]);
    xlabel(['Spanish ' label]);
    ylabel(['English ' label]);

    [r, p] = corr(x_all, y_all, 'Rows', 'complete');
    title({['r= ' num2str(r) ','], ['p=' num2str(p, 4)]});
     % Increment the subplot counter
    ctr = ctr + 1;
end
set(gca, 'FontSize', 13);

% Link axes for phonetic features and word onset/surp
linkaxes(ax(1:length(labels)));h = refline(1, 0);

% Add reference lines and lines at 0
for i = 1:length(labels)
    subplot(length(labels), 1, i);
    xline(0, 'Color', 'k', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
    xlim([-0.04 0.06]);
    xticks([0 0.06]);
    ylim([-0.04 0.06]);
    yticks([0 0.06]);
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
end

% add inlaid quadrant count on the top right corner
% add a inlaid quadrant plot
ax = axes('Position',[.65 .65 .3 .3]);
quads = rot90(histcounts2(x_all, y_all, [-100 0 100], [-100 0 100]));
imagesc(quads);
% change the colormap of this axis
colormap(ax, flipud(gray));
% colormap([1 1 1; 1 0 0; 1 0 1;0 0 1])
% add in the text overlaid
quads = flipud(rot90(quads));
for x = 1:2
    for y = 1:2
        text(x, y, num2str(quads(x, y)), 'Color', 'w', ...
            'FontSize', 13, 'HorizontalAlignment', 'center');
    end
end
xticks([1 2]);
yticks([1 2]);
xticklabels({'-UV', '+UV'});
yticklabels({'+UV', '-UV'});

% add the number of elecs in each quadrant
%% G - Bilingual uv as scatter v2


% Specify the labels for analysis
labels = {'word+surp'};

% Create a figure for the scatter plots
figure('Position', [200, 100, 900, 900]);

% Set a UV threshold for filtering
uv_thresh = 0.001;
pval_thresh = 0.05;

% Initialize subplot counter and axis array
ctr = 1;
ax = nan(3, 1);

cols = getColorsCrossComp(6);
% cols = [256 256 256; 28 117 188; 8 70 128;]./256;
% make a patch in the lower right quadrant thats light red
patch([0 0 0.08 0.08], [-0.02 0.08 0.08 -0.02], cols(1, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);
% make a patch in the upper left quadrant thats light blue
patch([0.08 0.08 -0.02 -0.02], [0 0.08 0.08 0], cols(2, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35); hold on;
%patch([0 0 -0.01 -0.01], [0 0.08 0.08 0], cols(2, :), 'EdgeColor', 'none', 'FaceAlpha', 0.35);

% Iterate over the labels
for label = labels
    x_all = [];
    y_all = [];
    
    % Find the index of the current label in the featureOrd array
    index = find(ismember(wordsurp_details.featureOrd, label));
    
    % Iterate over the ls values
    ls = 4;

    % Create a subplot for the current label and ls value
    ax(ctr) = subplot(length(labels), 1, ctr);
    lsidx = wordsurp_encoding.ls==ls;
    
    % Get the corresponding unique variance values for English and Spanish
    x = wordsurp_encoding.eng_uv_all(lsidx, index);
    y = wordsurp_encoding.sp_uv_all(lsidx, index);
    pvals_x = wordsurp_encoding.eng_wordsurp_pval(lsidx);
    pvals_y = wordsurp_encoding.sp_wordsurp_pval(lsidx);
    
    % Get the subject IDs for the current ls value
    sid = cellfun(@(x) str2double(x(3:end)), ...
        wordsurp_encoding.SID(wordsurp_encoding.ls==ls));
    
    % Remove data points that do not meet the UV threshold or contain NaN values
    % all([x,y]<uv_thresh, 2)
    neg = any([pvals_x,pvals_y]>pval_thresh, 2) | isnan(x) | isnan(y);
    %neg = all([x,y]<uv_thresh, 2) | isnan(x) | isnan(y);
    x(neg) = [];
    y(neg) = [];
    sid(neg) = [];

    % Perform permutation testing to compute correlation coefficient and p-value
    maxlim = prctile([x; y], 100);
    minlim = prctile([x; y], 3);  

    % Plot the scatter plot
    colors = x-y;
    x_all = [x_all; x];
    y_all = [y_all; y];
    %scatter3(x, y, sid, 20, colors, 'filled', ...
    %            'MarkerEdgeColor', 'k', 'MarkerFaceAlpha', 0.8); hold on;
    
    scatter(x, y, 15, colors, 'filled', ...
                'MarkerEdgeColor', 'k', ...
                'MarkerFaceAlpha', 0.8, ...
                'LineWidth', 0.25); hold on;
    
    % colorbar should be blue to red going through white (1, 1, 1)
    % first create blue to white
    cols = [linspace(43/256, 1, 50); linspace(122/256, 1, 50); linspace(186/256, 1, 50)]';
    % then add white to red
    cols = [cols; [linspace(1, 213/256, 50); linspace(1, 37/256, 50); ...
        linspace(1, 39/256, 50)]'];

    colormap(cols);
    clim([-0.05 0.05]);
    view(2);
    
    % Set the x and y axis limits and add labels
    xlim([minlim maxlim]);
    ylim([minlim maxlim]);
    xlabel(['English ' label{1} ' \Delta R^2']);
    ylabel(['Spanish ' label{1} ' \Delta R^2']);

    [r, p] = corr(x_all(x_all>0&y_all>0), y_all(x_all>0&y_all>0), 'Rows', ...
        'complete', 'type', 'Spearman');
    title({['r= (' num2str(sum(x_all>0&y_all>0)) '),' num2str(r)], ['p=' num2str(p, 4)]});
     % Increment the subplot counter
    ctr = ctr + 1; 
    hold on;
end

% Link axes for phonetic features and word onset/surp
linkaxes(ax(1:length(labels)));
h = refline(1, 0);
set(gca, 'FontSize', 15);

% Add reference lines and lines at 0
for i = 1:length(labels)
    subplot(length(labels), 1, i);
    xline(0, 'Color', 'k', 'LineWidth', 1.5);
    yline(0, 'Color', 'k', 'LineWidth', 1.5);
    xlim([-0.02 0.065]);
    ylim([-0.02 0.065])
    xticks([0 0.08]);
    yticks([0 0.08]);
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    % make colors all black
    colormap(repmat([0 0 0], 256, 1));
end

% add inlaid quadrant count on the top right corner
% add a inlaid quadrant plot
ax = axes('Position',[.65 .65 .3 .3]);
quads = rot90(histcounts2(x_all, y_all, [-100 0 100], [-100 0 100]));

imagesc([1, 3; 0, 2]);
% change the colormap of this axis
cols = getColorsCrossComp(6);
cols = [1 1 1; cols(2, :); cols(1, :); mean(cols)];
cols = brighten(cols, 0.8);
colormap(ax, cols);

% colormap(ax, flipud(gray));
%colormap([1 1 1; 1 0 0; 1 0 1;0 0 1])
% add in the text overlaid
quads = flipud(rot90(quads));
for x = 1:2
    for y = 1:2
        text(x, y, num2str(quads(x, y)), 'Color', 'k', ...
            'FontSize', 16, 'HorizontalAlignment', 'center');
    end
end
box off;
xticks([1 2]);
yticks([1 2]);
xticklabels({'-UV', '+UV'});
yticklabels({'+UV', '-UV'});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;




%% ----------------------- Supplementary Figures --------------------------

%% S2 - Coverage of electrodes across anatomical regions 

SIDs = bSIDs;
dimex_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
timit_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');
asccd_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat');

all_elecs = {dimex_elecs, timit_elecs, asccd_elecs};

label_ord = {'middletemporal', 'superiortemporal', 'supramarginal', 'precentral', ...
    'postcentral', 'inferiorparietal', 'inferiortemporal', 'lateralorbitofrontal', ...
    'parstriangularis', 'fusiform', 'caudalmiddlefrontal', 'parsopercularis', ...
    'rostralmiddlefrontal', 'lateraloccipital', 'temporalpole'};
label_tick = {'middletemp', 'superiortemp', 'supramarg', 'precent', 'postcent', ...
    'infpar', 'inftemp', 'latorbfront', 'parstriang', 'fusif', 'caudmidfront', 'parsoperc', ...
    'rostmidfront', 'latocc', 'temppole'};

for h = {'lh', 'rh'}
    figure;
    set(gcf,'Color','w');
    axes;
    hold on;
    cmax = 2;

    hemi = h{1};
    if strcmp(hemi, 'lh')
        x_add = -20;
    else
        x_add = 10;
    end
    
    numsid = 0;      
    if strcmp(hemi, 'lh')
        cortex = imgall.(SIDs{1}).img_mni.cortex;
    else
        cortex = imgall.(SIDs{3}).img_mni.cortex;
    end
    
    % find density map of native speech, lateral side
    native_xyz = [];
    native_hga = [];
    all_xyz = [];
    anatomy = [];
    speech_resp = [];
    for si = SIDs
        sid = si{1};
        if strcmpi(imgall.(sid).hemi,hemi)
            subject = imgall.(sid);
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix; 
                numel = size(elecmatrix, 1);
                all_xyz = [all_xyz; elecmatrix]; 

                % note speech responsive (union electrodes)
                speech_resp_tmp = zeros(size(elecmatrix, 1), 1);
                for j = 1:length(all_elecs)
                    if isfield(all_elecs{j}.allidx, sid)
                        speech_resp_tmp(all_elecs{j}.allidx.(sid)) = 1;
                    end
                end
                
                if isfield(subject.img_native, 'anatomy')
                    anatomy =[anatomy; subject.img_native.anatomy(1:numel, 4)];
                elseif isfield(subject.imgNative, 'anatomy')
                    anatomy = [anatomy; subject.imgNative.anatomy(1:numel, 4)];
                else
                    disp(['using mni anatomy labels... for ' sid])
                    anatomy = [anatomy; subject.img_mni.anatomy(1:numel, 4)];
                end
                speech_resp = [speech_resp; speech_resp_tmp];
                numsid=numsid+1;
            else
                disp(['Missing...' sid]);
            end
        end
    end

    subplot(2, 1, 1);
    PlotBrainSurface(cortex, hemi,'lateral');
    alpha 0.9
    light("Style","infinite","Position",[100 100 0]);
    yl = ylim;
    zl = zlim;

    % number of subjects
    disp(['Number of subjects: '  num2str(numsid)]);
    disp(['Number of electrodes: '  num2str(length(all_xyz))]);
    disp('--------------------------------------------------------');
    scatter3(all_xyz(:, 1)+x_add, all_xyz(:, 2), all_xyz(:, 3), 3, 'k', 'filled');
    zlim(zl); ylim(yl);

    % Generate coverage contour
    ha = subplot(2, 1, 2);
    PlotBrainSurface(cortex, hemi,'lateral');
    alpha 0.4
    light("Style","infinite","Position",[100 100 0]);
    drawBrainContour(ha, all_xyz, hemi, cmax);
    disp(numsid);
   

    % remove NaN and Unknown, Right-Cerebral-WhiteMatter
    torem = {'NaN', 'Unknown', 'Right-Cerebral-White-Matter', ...
                'Left-Cerebral-White-Matter', 'Left-Hippocampus', ...
                'Right-Hippocampus'};
    % remove all cases where coveage < 20
    mincount = 20;

    [counts, ~, ~, labels] = crosstab(anatomy, speech_resp);
    labels = labels(:, 1);

    % remove all cases where coverage < 20
    labels = labels(sum(counts, 2)>mincount);
    counts = counts(sum(counts, 2)>mincount, :);
    
    % remove NaN and Unknown, Right-Cerebral-WhiteMatter
    counts = counts(~ismember(labels, torem), :);
    labels = labels(~ismember(labels, torem));
    
    idx_cell = cellfun(@(x) find(ismember(labels, x)), ...
        label_ord, 'UniformOutput', false);
    counts_reord = nan(length(label_ord), 2);
    for i = 1:length(label_ord)
        if isempty(idx_cell{i})
            counts_reord(i, :) = 0;
        else
            counts_reord(i, :) = counts(idx_cell{i}, :);
        end
    end
end

% make a dictionary of anatomical labels all set to empty arrays
% this will keep track of the difference between number of speech responsive sites in one particular area
anatomy_diff = containers.Map(label_ord, cell(length(label_ord), 1));

ctr = 1;
figure;
for h = {'lh', 'rh'}
    hemi = h{1};
    numsid = 0;     
    
    % find density map of native speech, lateral side
    all_xyz = [];
    anatomy = [];
    speech_resp = [];
    for si = SIDs
        sid = si{1};
        if strcmpi(imgall.(sid).hemi,hemi)
            subject = imgall.(sid);
            ls = find(cellfun(@(x) ismember(sid, x), {sSIDs, eSIDs, mSIDs, bSIDs}));
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix; 
                numel = size(elecmatrix, 1);
                all_xyz = [all_xyz; elecmatrix]; 

                % note speech responsive (union electrodes)
                speech_resp_tmp = zeros(size(elecmatrix, 1), 2);
                for j = 1:length(all_elecs)
                    native = j;
                    if isfield(all_elecs{j}.allidx, sid)
                        speech_resp_tmp(all_elecs{j}.allidx.(sid), native) = 1;
                    end
                end
                
                if isfield(subject.img_native, 'anatomy')
                    anat_tmp = subject.img_native.anatomy(1:numel, 4);
                elseif isfield(subject.imgNative, 'anatomy')
                    anat_tmp = subject.imgNative.anatomy(1:numel, 4);
                else
                    disp(['using mni anatomy labels... for ' sid])
                    anat_tmp = subject.img_mni.anatomy(1:numel, 4);
                end
                speech_resp = [speech_resp; speech_resp_tmp];
                anatomy = [anatomy; anat_tmp];

                % update the dictionary
                for i = 1:length(label_ord)
                    idx = ismember(anat_tmp, label_ord{i});
                    if sum(idx)>0
                        % append the difference between the two speech responsive sites to array
                        anatomy_diff(label_ord{i}) = [anatomy_diff(label_ord{i}); ...
                            sum(speech_resp_tmp(idx, 1) - speech_resp_tmp(idx, 2))];
                    end
                end
                numsid=numsid+1;
            else
                disp(['Missing...' sid]);
            end
        end
    end

    % number of subjects
    disp(['Number of subjects: '  num2str(numsid)]);
    disp(['Number of electrodes: '  num2str(length(all_xyz))]);
    disp('--------------------------------------------------------');

    % remove NaN and Unknown, Right-Cerebral-WhiteMatter
    torem = {'NaN', 'Unknown', 'Right-Cerebral-White-Matter', ...
                'Left-Cerebral-White-Matter', 'Left-Hippocampus', ...
                'Right-Hippocampus'};
    % remove all cases where coveage < 20
    mincount = 20;

    speech_resp_type = speech_resp(:, 1) + 2*speech_resp(:, 2);
    [counts, ~, ~, labels] = crosstab(anatomy, speech_resp_type);
    labels = labels(:, 1);
    
    % remove all cases where coverage < 20
    labels = labels(sum(counts, 2)>mincount);
    counts = counts(sum(counts, 2)>mincount, :);
    
    % remove NaN and Unknown, Right-Cerebral-WhiteMatter
    counts = counts(~ismember(labels, torem), :);
    labels = labels(~ismember(labels, torem));
    
    idx_cell = cellfun(@(x) find(ismember(labels, x)), ...
        label_ord, 'UniformOutput', false);
    counts_reord = nan(length(label_ord), 4);
    for i = 1:length(label_ord)
        if isempty(idx_cell{i})
            counts_reord(i, :) = 0;
        else
            counts_reord(i, :) = counts(idx_cell{i}, :);
        end
    end
    
    subplot(1, 2, ctr)
    b = barh(counts_reord, 'stacked', 'EdgeColor', 'none'); hold on;
    cols = getColorsCrossComp(6);
    b(1).FaceColor = [0.8 0.8 0.8];
    b(2).FaceColor = cols(2, :);
    b(3).FaceColor = cols(1, :);
    b(4).FaceColor = [0.6 0.1 0.6];
    if ctr==1
        set(gca, 'YTick', 1:length(label_ord), 'YTickLabel', label_tick , 'FontSize', 13);
    else
        set(gca, 'YTick', [], 'FontSize', 13);
    end
    ylabel('anatomy');
    xlabel('count');
    xlim([0 300]);
    ytickangle(30);
    xticks([0 300]);
    box off;
    legend({'non-responsive', 'Spanish', 'English', 'both'});
    ctr = ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% Single example subject response magnitudes, example electrode

%'EC100', 'EC252', 'EC152', 'EC212', 'EC235', 'EC129', 'EC159', 'EC196'
nativeSIDs = {'EC214'};

% Color by language
cm = [0 0 1; 1 0 0];
modelname={'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; 
modelfeatures  = [{'onset'}; timit_details.features.names([1:3, 8, 9, 11]); ... 
    {'peakrate'; 'F1'; 'F2'; 'F3'; 'F4'}];

bins = 15;
% elecs = 203; % for EC183 , 156, 160
%elecs = 137; % for EC172
elecs = [3, 54]; % for EC214
%elecs = [71, 86]; % for HS11
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
%     % code to find electrodes with high response
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
        clim([0 3]);
        colormap(ax, flipud(blues));

        ax = subplot(2, 2, 3);
        imagesc(smoothdata(squeeze(en_resp(1, :, 1:min_sent))', 2)); 
        xline(150, 'LineWidth', 2); 
        xline(50, 'LineWidth', 2, 'LineStyle', '--');
        xline(200, 'LineWidth', 2, 'LineStyle', '--');
        colormap(ax, flipud(reds));
        clim([0 3]);
        sgtitle([SID ': ' num2str(sent_encoding.el(el)) ', ' num2str(el)])

        % imagesc the encoding models
        ax = subplot(2, 2, 2);
        if startsWith(SID, 'EC')
            corp = 'dimex';
        else
            corp = 'asccd';
        end
        [weights] = getTRFweights(SID, sent_encoding.el(el), ...
                corp, modelname, datapath);
        imagesc(weights(:, 1:40));
        colormap(ax, inferno);
        clim([-1.5 1.5]);
        title(corp);
        clear corp

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
        clim([-1.5 1.5]);
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

%% Functions
function [weights] = getTRFweights(SID, el, corpus, modelname, datapath)    
    [strf] = loadMultModelStrf(SID, modelname, corpus, datapath, 1);  
    weights = strf{1}.meanStrf(:, :, el);
end
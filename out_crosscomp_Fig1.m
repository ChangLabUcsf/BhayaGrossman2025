% Ilina Bhaya-Grossman
% 01.08.2022
% SIDs = [sSIDs, eSIDs];
out_crosscomp_startup;

% load in all beta model versions
SIDs = [{'HS8', 'HS9', 'HS10'}, sSIDs, eSIDs]; % , {'HS11', 'HS9', 'HS10'}
tps = 50:55;

% load in sentence responses
if ~exist('sent_encoding', 'var')
    load([datapath '/Figure1/Figure1_SentResp.mat'], 'sent_encoding');
end

if ~exist('imgall', 'var')
    load([datapath '/Figure1/Figure1_ImgData.mat'], 'imgall');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    *encoding* allidx fthresh Dcons *wrd*;

%% Determine how many minutes of speech each participant heard

load('data/Figure1/Figure1_CorpEvnt.mat')

% print the mean and range for each language
disp(['Mean minutes of English speech heard: ' num2str(mean(cell2mat(struct2cell(timin))))]);
disp(['Min of minutes of English speech heard: ' num2str(min(cell2mat(struct2cell(timin))))]);
disp(['Max of minutes of English speech heard: ' num2str(max(cell2mat(struct2cell(timin))))]);
disp('--------------------------------------------------------');
disp(['Mean minutes of Spanish speech heard: ' num2str(mean(cell2mat(struct2cell(dimin))))]);
disp(['Min of minutes of Spanish speech heard: ' num2str(min(cell2mat(struct2cell(dimin))))]);
disp(['Max of minutes of Spanish speech heard: ' num2str(max(cell2mat(struct2cell(dimin))))]);
disp('--------------------------------------------------------');
disp(['Mean minutes of Mandarin speech heard: ' num2str(mean(cell2mat(struct2cell(ascdn))))]);
disp(['Min of minutes of Mandarin speech heard: ' num2str(min(cell2mat(struct2cell(ascdn))))]);
disp(['Max of minutes of Mandarin speech heard: ' num2str(max(cell2mat(struct2cell(ascdn))))]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Display hemisphere and number of participants in each group

% display hemisphere and number of participants in each group
disp('Hemisphere and number of participants in each group');
labels = {'English', 'Spanish', 'Mandarin'};
ctr = 1;
for grp = {eSIDs, sSIDs, {'HS8', 'HS9', 'HS10'}}
    grp_sids = grp{1};
    lh_sids = unique(sent_encoding.SID(strcmpi(sent_encoding.hemi, 'lh') & ...
        ismember(sent_encoding.SID, grp_sids)));
    rh_sids = unique(sent_encoding.SID(strcmpi(sent_encoding.hemi, 'rh') & ...
        ismember(sent_encoding.SID, grp_sids)));

    % number of participants in each group
    disp([labels{ctr} ' monolinguals : ' num2str(length(grp_sids))]);
    % number of participants in each hemisphere
    disp(['LH: ' num2str(length(lh_sids)) ' (' strjoin(lh_sids, ' ') ')']);
    disp(['RH: ' num2str(length(rh_sids)) ' (' strjoin(rh_sids, ' ') ')']);
    disp('--------------------------------------------------------');
    ctr = ctr+1;
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Info about participant neural responses across languages

% how many electrdes are responsive to speech in both languages
blang = all(sent_encoding.type, 2); % when both languages are sig
disp(['Both languages: ' num2str(sum(blang)) ', ' ...
    num2str((sum(blang)/height(sent_encoding))*100) '%']);

% for each subject, identify the percentage of electrodes that are responsive to both languages
perc = nan(length(unique(sent_encoding.SID)), 1);
ctr = 1;
for s = unique(sent_encoding.SID)'
    sid = s{1};
    idx = strcmp(sent_encoding.SID, sid);
    blang = all(sent_encoding.type(idx, :), 2);
    perc(ctr) = sum(blang)/sum(idx);
    disp([sid ': ' num2str(sum(blang)) ', ' ...
        num2str((sum(blang)/sum(idx))*100) '%']);
    ctr = ctr+1;
end

disp(['Median: ' num2str(median(perc)*100) '%']);
disp(['STD: ' num2str(std(perc)*100)]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% B/C - Single example subject response magnitudes, example electrode

nativeSIDs = {'EC100'};
% Color by language
cm = [0 0 1; 1 0 0];

bins = 15;
% elecs = 203; % for EC183 , 156, 160
% elecs = 137; % for EC172
elecs = 22; % for EC100 % [3, 54]; % for EC214
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
    for s=unique(sent_encoding.SID)'
        SID = s{1};
        idx = strcmp(sent_encoding.SID, SID) & sent_encoding.type(:, l)==1;
        desel.(SID).elid = sent_encoding.el(idx);
        desel.(SID).condition = conds(idx);
        desel.(SID).selid = elecs; %sent_encoding.el(idx);
    end
    plotNativeElec(nativeSIDs, desel, 1, imgall);
end

for sid = nativeSIDs
    SID = sid{1};   

    els = find(strcmp(sent_encoding.SID, SID) ...
        & ismember(sent_encoding.el, elecs))';
    [en_mresp, sp_mresp] = plotStitchedSentence(els, sent_encoding, 150, 100, 0.5, 0);

    [r, p] = corr(en_mresp', sp_mresp');
    sgtitle([SID ': ' num2str(sent_encoding.el(elecs(1))) ...
        ', ' num2str(elecs(1))])
    disp(['Pearson r(' num2str(length(en_mresp)) ') = ' num2str(r) ...
        ', p=' num2str(p)]);
end

for i = 1:2
    subplot(1, 2, i);
    xlim([-0.5 2.15]);
    ylim([-0.5 4.5]);
    yticks([0 4]);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% D - Percent speech responsive by language and native / foreign

% electrode selection details
timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");
asccd_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat");

figure;
native_elecs = {timit_elecs, dimex_elecs, asccd_elecs}; % Spanish, English, Mandarin
foreign_elecs = {dimex_elecs, timit_elecs, timit_elecs}; % English, Spanish, English
plts = [2, 1, 3];
for ls = [1, 2, 3] % languages : spanish, english, mandarin
    
    % find number of speech responsive electrodes in table
    idx = sent_encoding.ls == ls;
    sids = unique(sent_encoding.SID(idx));
    % 1st column : native, 2nd column : foreign
    numspeech = zeros(length(sids), 2);
    numel = zeros(length(sids), 2);

    % for each subject, identify the number of total elecs
    ctr = 1;
    for s = sids'
        sid = s{1};
        numspeech(ctr, :) = [length(native_elecs{ls}.allidx.(sid)), ...
            length(foreign_elecs{ls}.allidx.(sid))];
        numel(ctr, :) = [length(native_elecs{ls}.fvals.(sid)), ...
            length(foreign_elecs{ls}.fvals.(sid))];
        ctr = ctr + 1;
    end

    % plot the pie chart for each language group
    subplot(2, 3, plts(ls));
    % explode and show percentage for native and foreign
    pc = pie([sum(numspeech(1, :)), sum(numel(1, :))-sum(numspeech(1, :))], ...
        [0 1]);
    % no edge color, light blue 
    pc(1).FaceColor = [0.6 0.6 1];
    pc(3).FaceColor = [0.6, 0.6, 0.6];
    pc(1).EdgeColor = 'none';
    pc(3).EdgeColor = 'none';
    pc(2).FontSize = 14;
    pc(4).FontSize = 14;

    subplot(2, 3, plts(ls)+3); % foreign
    pc = pie([sum(numspeech(2, :)), sum(numel(2, :))-sum(numspeech(2, :))], ...
        [0 1]);
    % no edge color, light red
    pc(1).FaceColor = [1 0.6 0.6];
    pc(3).FaceColor = [0.6, 0.6, 0.6];
    pc(1).EdgeColor = 'none';
    pc(3).EdgeColor = 'none';
    pc(2).FontSize = 14;
    pc(4).FontSize = 14;

    if ls == 3
        legend({'Speech responsive', 'Non-responsive'}, 'Location', 'southoutside');
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% E - Single scatter plot of speech responsive electrodes in native and foreign

anat_counts = cell(2, 1);

figure;
for h = {'lh', 'rh'}
    hemi = h{1};
    set(gcf,'Color','w');

    if strcmp(hemi, 'lh')
        x_add = -20;
    else
        x_add = 10;
    end

    types = [0, 1];
    for native = types
        numsid = 0;      

        if strcmp(hemi, 'lh')
            cortex = imgall.(SIDs{1}).img_mni.cortex;
            sbplt=1;
        else
            cortex = imgall.(SIDs{6}).img_mni.cortex;
            sbplt=2;
        end

        % lh+native = 1, rh+native = 2
        % lh+foreign = 3, rh+foreign = 4
        subplot_idx = (1 - native)*2 + sbplt;
        subplot(2, 2, subplot_idx);
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
                    elecmatrix = subject.img_mni.elecmatrix; 
                    if isfield(subject.img_native, 'anatomy')
                        anat = subject.img_native.anatomy(:, 4);
                    elseif isfield(subject.imgNative, 'anatomy')
                        anat = subject.imgNative.anatomy(:, 4);
                    else
                        disp(['using mni anatomy labels... for ' sid])
                        anat = subject.img_mni.anatomy(:, 4);
                    end
                    
                    ch_sid = sent_encoding.el(strcmp(sent_encoding.SID, sid));
                    ls = sent_encoding.ls(find(strcmp(sent_encoding.SID, sid), 1, 'first'));
            
                    if native % spanish for bilinguals
                        if ismember(ls, [1, 2])
                            ch_sel = sent_encoding.el(sent_encoding.type(:, ls)==1);
                        else % bilingual case and mandarin case (where first column is mandarin)
                            ch_sel = sent_encoding.el(sent_encoding.type(:, 1)==1);
                        end
        
                        % for weighting
                        if ismember(ls, [1, 3])
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
                        if ismember(ls, [1, 3])
                            native_hga = [native_hga; ...
                                sent_encoding.maxresp(intersect(ch_sel, ch_sid), 2)];
                        elseif ismember(ls, 2)
                            native_hga = [native_hga; ...
                                sent_encoding.maxresp(intersect(ch_sel, ch_sid), 1)];
                        end
                    end
        
                    native_xyz = [native_xyz; elecmatrix(intersect(ch_sel, ch_sid),:)]; 
                    anat_xyz = [anat_xyz; anat(intersect(ch_sel, ch_sid))];
                end
                numsid=numsid+1;
            end
        end

        % number of subjects
        disp(['Is native: ' num2str(native)]);
        disp(['Number of subjects with anatomy labels: '  num2str(numsid)]);
        disp(['Number of electrodes: '  num2str(length(native_hga))]);

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

        if ~native
            color = [1, 0, 0];
        else
            color = [0, 0, 1];
        end
        scatter3(native_xyz(:, 1)+x_add, native_xyz(:, 2), native_xyz(:, 3), ...
            round(native_hga*20), color, 'filled', 'MarkerFaceAlpha', 0.7);
    end
end

% Create horizontal bar charts for speech-responsive areas
figure;
subplot(1, 2, 1);
% remove areas with count < 5 and sort by count
anat_counts{1} = anat_counts{1}(anat_counts{1}.count >= 5, :);
anat_counts{1} = sortrows(anat_counts{1}, 'count', 'descend');
barh(anat_counts{1}.area, anat_counts{1}.count, 'FaceColor', ...
    [1, 0, 0], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Red color for non-native
title('Foreign');
xlabel('Count');
ylabel('Anatomical Area');
set(gca, 'FontSize', 14);
xticks([0 600]);

box off;

subplot(1, 2, 2);
% remove areas with count < 5 and sort by count
anat_counts{2} = anat_counts{2}(anat_counts{2}.count >= 5, :);
anat_counts{2} = sortrows(anat_counts{2}, 'count', 'descend');
barh(anat_counts{2}.area, anat_counts{2}.count, 'FaceColor', ...
    [0, 0, 1], 'FaceAlpha', 0.7, 'EdgeColor', 'none'); % Blue color for native
title('Native');
xlabel('Count');
ylabel('Anatomical Area');
set(gca, 'FontSize', 14);
sgtitle('Speech-responsive areas');
xticks([0 600]);
xlim([0 600]);
box off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% F - Example electrode TRF weights

% example electrodes in the figure
% SID = 'EC172';
% els = [74,  136];
SID = 'EC100';
els = [70 71];

% elecs('EC183') = [72];
% elecs('EC172') = [122];
% elecs('EC100') = 118;

corpusnames = {'dimex', 'timit'};
details = {dimex_details, timit_details};
bef=20;
betacolor = flipud(prgn);

modelname={'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; 
modelfeatures  = [{'onset'}; timit_details.features.names([1:3, 8, 9, 11]); ... 
    {'peakrate'; 'F1'; 'F2'; 'F3'; 'F4'}];

for el = els
    figure; 
    weights = cell(2, 2);
    for corp = 1:2
        subplot(2, 2, 1 + (corp-1)*2);

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
        colormap(betacolor)

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
        colormap(betacolor);
        clim([-1.8 1.8]);

        % formatting
        yticks(1:length(modelfeatures));
        yticklabels(modelfeatures);
        xlim([0.5 40])
        xticks([1 40]);
        xticklabels({'0', '-0.4'});
        xlabel('Time (s)');
        set(gca, 'FontSize', 13);
        clear strf
        box off;
        
        % more formatting
        subplot(2, 2, 1 + (corp-1)*2);
        yline(0, 'LineWidth', 1.5);
        xline(0, 'LineWidth', 1.5);
        xlabel('Time (s)');
        set(gca, 'FontSize', 13);
        ylabel('HFA (z)');   
        box off;
    end
    sgtitle([num2str(el) ' ' num2str(corr(weights{2, 1}(:), ...
        weights{2, 2}(:)))])
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% G/H - Correlating (s)TRF weights across corpora

corpora = {{'timit', 'dimex'}, {'timit', 'asccd'}};
%modelnames = {'onset_aud'};
modelnames = {'onset_phnfeatConsOnset_maxDtL_formantMedOnset'};

% consonant feature overlap for mandarin and english
asccd_feats = asccd_details.features.names(1:5);
timit_feats = timit_details.features.names([1:3, 8, 9, 11]);
shared_feats = intersect(timit_feats, asccd_feats);
asccd_fidx = [1; 1+cellfun(@(x) find(strcmp(asccd_feats, x)), shared_feats); ...
    (11-4:11)'];
timit_fidx = [1; 1+cellfun(@(x) find(strcmp(timit_feats, x)), shared_feats); ...
    (12-4:12)'];

corrstrf = nan(height(sent_encoding), 1);
corrstrfperm = nan(height(sent_encoding), 10);
corrpval = corrstrf;
maxfreq = nan(height(sent_encoding), 2);
maxrsq = nan(height(sent_encoding), 1);
minrsq = nan(height(sent_encoding), 1);
thresh = 0.001;

% window is 5 time points around the max time point for correlation
windsz = 2;

% all windowed TRF weights
alleng = nan(height(sent_encoding), 80);
allsp = nan(height(sent_encoding), 80); % maximum feat length

for s =  unique(sent_encoding.SID)'
    
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
            sent_encoding.sp_base_rsq(e)>thresh) && ~missing % ...
             % &&  ~strcmp(sent_encoding.anatomy(e), 'superiortemporal')
    
            % add in TRF weight correlation
            [~, maxtp] = max(squeeze(mean(meanStrf(:, :, :, el), [1 2])));
            wind = max(1, maxtp-windsz):min(maxtp+windsz, size(meanStrf, 3));
        
            % compare frequency selectivity
            x = squeeze(mean(meanStrf(1, :, wind, el), 3));
            y = squeeze(mean(meanStrf(2, :, wind, el), 3));
            
            % scale and smooth the STRF weights to be more comparable
            [~, z] = procrustes(x', y');
            if strcmp(modelnames{1}, 'onset_aud') % only for STRF case (not in current figure)
                z = smoothdata(z, 'SmoothingFactor',0.4);
                x = smoothdata(x, 'SmoothingFactor',0.1);
            end 
            
            % save out weights
            alleng(e, 1:length(x)) = x;
            allsp(e, 1:length(z)) = z;

            % calculate correlation and max frequency
            % Pearson without procrustes
            [corrstrf(e), corrpval(e)] = corr(x', y', 'type', 'Pearson');

            rng(1);
            for j = 1:10
                [corrstrfperm(e, j), ~] = corr(x', y(randperm(length(y)))', ...
                    'type', 'Pearson');
            end
            [~, maxfreq(e, :)] = max([x', y']);
            %freq(e, :, :) = [x'; z];

            % save out max R^2 and min R^2
            maxrsq(e) = max(sent_encoding.en_base_rsq(e), ...
                sent_encoding.sp_base_rsq(e));
            minrsq(e) = min(sent_encoding.en_base_rsq(e), ...
                sent_encoding.sp_base_rsq(e));
        
            debug = 0;
            % Example electrodes from above (EC100, 22 / 150)
            if debug && startsWith(SID, 'HS') && minrsq(e)>0.1
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

% show imagesc of TRF weights comparing the two corpora
% find nan rows for alleng and allsp
if ~strcmp(modelnames{1}, 'onset_aud')
    alleng = alleng(:, 1:12);
    allsp = allsp(:, 1:12);
end

nanrows = any(isnan(alleng), 2) | any(isnan(allsp), 2) | minrsq<0.1;
alleng(nanrows, :) = [];
allsp(nanrows, :) = [];
% alleng = alleng(idx, :);
% allsp = allsp(idx, :);

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

[r, p] = corr(allsp(:), alleng(:));
disp(['Correlation between beta matrices: ' ...
    'r(' num2str(length(allsp(:))) ') = ' num2str(r) ', p = ' num2str(p)]);

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
hbins = -1:0.05:1;

% permutation distribution
idx = minrsq>0.1 & ~nanrows;
corrstrfperm(~idx, :) = NaN; 
corrstrfperm(isnan(corrstrfperm)) = [];

histogram(corrstrfperm(:), hbins, 'EdgeColor', 'none', 'Normalization', 'pdf', ...
    'FaceColor', [0.7, 0.7, 0.7], 'FaceAlpha', 0.3);
hold on;
% make a kdensity plot
[f, xi] = ksdensity(corrstrfperm(:));
plot(xi, f, 'LineWidth', 3, 'Color', [0.7, 0.7, 0.7], 'LineStyle', '-');
yticks(0:0.5:1.5);

yyaxis right; 
histogram(corrstrf(idx), hbins, 'EdgeColor', 'none', 'Normalization', 'pdf', ...
    'FaceColor', [0.7, 0.1, 0.7], 'FaceAlpha', 0.3);
hold on;
[f, xi] = ksdensity(corrstrf(idx));
plot(xi, f, 'LineWidth', 3, 'Color', [0.7, 0.1, 0.7], 'LineStyle', '-');
stathresh = (0.05)*100;

xline(prctile(corrstrfperm(:), 100-stathresh), 'LineWidth', 2, 'Color', 'k');
xlim([-1 1]);
yticks(0:2:6);
ax = gca();
ax.YAxis(2).TickLabelColor = [0.7, 0.1, 0.7]; 
ax.YAxis(2).Color = [0.7, 0.1, 0.7]; 
yyaxis left;
ylabel('Probability Density');
box off;
set(gca, 'FontSize', 13);

disp(['Total number of electrodes: ' num2str(sum(idx))])
disp(['Number of electrodes with higher than 95% for permuted distrib: ' ...
    num2str(sum(corrstrf(idx)>prctile(corrstrfperm(:), 100-stathresh)))]);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq;

%% ----------------------- Supplementary Figures --------------------------
%% ADDF0: S1 - Spectral and Temporal MTFs for all stimuli
% Note: Will take about two minutes to calculate

corpora = {timit_details, dimex_details, asccd_details};
corpuslabels = {'timit', 'dimex',  'asccd'};

figure;
for c = 1:length(corpora)

    % Load the corpus details
    corpus_details = corpora{c};
    corpus = corpuslabels{c};

    % Create two NaN matrices of size 450x4096
    p = nan(450, 4096);
    d = nan(450, 4096);

    % Create a NaN matrix of size 100x257x512
    mtf = nan(100, 257, 512);
    ctr = 1;

    % Iterate over the first 200 sentences
    for i = 1:58 % minimum number of sentences
        fs = 16000;
        % fs = 44100;

        % Calculate the power spectrum of the sound in the i-th sentence
        [p(i, :),~] = pspectrum(corpus_details.sentdet(i).sound, fs);

        % Define start positions for segmentation
        strts = 50:60:size(corpus_details.sentdet(i).aud, 2)-120;

        % Add F0 (fundamental frequency) for the sentence from corpus_details
        f0 = addF0(corpus_details.sentdet(i), corpus);

        % Remove zero values from f0
        f0(f0==0) = [];

        % Check if the mean of f0 is greater than 170
        if mean(f0) > 170
            % Iterate over the start positions
            for j = 1:length(strts)

                % Extract a segment of audio data and compute the MTF
                segaud = corpus_details.sentdet(i).aud(:, strts(j):strts(j)+60);
                [mtf(ctr, :, :), xax, yax, btm, bsm] = get_MTF_bsm_btm(segaud, ...
                    corpus_details.sentdet(i).dataf, 8);

                % Increment the counter
                ctr = ctr + 1;
            end
        end
    end

    % Display the average MTF as an image
    subplot(1, 3, c);
    imagesc(xax, fliplr(yax), squeeze(mean(mtf, 'omitnan')));
    title(corpus);
    hold on;
    set(gca, 'Ydir', 'normal');
    xlim([-8 8]);
    ylim([0, 0.75]);
    colormap(inferno);
    brighten(0.6)

    % Add contour lines to the image
    % contour(xax, fliplr(yax), squeeze(mean(mtf)), 'Color', [0 0 0], 'LineWidth', 1.75);

    % Set labels and font size
    ylabel('Spectral modulation (cycles/oct)');
    xlabel('Temporal modulation (Hz)');
    set(gca, 'FontSize', 15);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% S2 - Coverage of electrodes across anatomical regions 

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
        cortex = imgall.(SIDs{6}).img_mni.cortex;
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

end


dimex_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
timit_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');
asccd_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat');

all_elecs = {dimex_elecs, timit_elecs, asccd_elecs};

% make a dictionary of anatomical labels all set to empty arrays
% this will keep track of the difference between number of speech responsive sites in one particular area
anatomy_diff = containers.Map(label_ord, cell(length(label_ord), 1));

for h = {'lh', 'rh'}
    hemi = h{1};
    numsid = 0;     
    
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
            ls = find(cellfun(@(x) ismember(sid, x), {sSIDs, eSIDs, mSIDs, bSIDs}));
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix; 
                numel = size(elecmatrix, 1);
                all_xyz = [all_xyz; elecmatrix]; 

                % note speech responsive (union electrodes)
                speech_resp_tmp = zeros(size(elecmatrix, 1), 2);
                for j = 1:length(all_elecs)
                    native = (ls == j)+1;
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
            end
            numsid=numsid+1;
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

    % [counts_foreign, ~, ~, labels_foreign] = crosstab(anatomy, speech_resp(:, 1));
    % [counts_native, ~, ~, labels_native] = crosstab(anatomy, speech_resp(:, 2));
    % counts = [counts_foreign, counts_native];
    % labels = [labels_foreign(:, 1), labels_native(:, 1)];
    % assert(all(strcmp(labels(:, 1), labels(:, 2))));
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
    
    figure;
    b = barh(counts_reord, 'stacked', 'EdgeColor', 'none');
    b(1).FaceColor = [0.8 0.8 0.8];
    b(2).FaceColor = [0.8 0.1 0.2];
    b(3).FaceColor = [0.2 0.1 0.8];
    b(4).FaceColor = [0.6 0.1 0.6];
    set(gca, 'YTick', 1:length(label_ord), 'YTickLabel', label_tick , 'FontSize', 13);
    ylabel('anatomy');
    xlabel('count');
    xlim([0 800]);
    box off;
    legend({'non-responsive', 'foreign', 'native', 'both'});

    % figure;
    % subplot(1, 2, 1)
    % b = barh(counts_reord(:, 1:2), 'stacked', 'EdgeColor', 'none');
    % b(1).FaceColor = [0.8 0.8 0.8];
    % b(2).FaceColor = [0.8 0.1 0.2];
    % set(gca, 'YTick', 1:length(label_ord), 'YTickLabel', label_tick , 'FontSize', 13);
    % ylabel('anatomy');
    % xlabel('count');
    % xlim([0 800]);
    % box off;
    % legend({'non-responsive', 'speech responsive'});

    % subplot(1, 2, 2)
    % b = barh(counts_reord(:, 3:4), 'stacked', 'EdgeColor', 'none');
    % b(1).FaceColor = [0.8 0.8 0.8];
    % b(2).FaceColor = [0.2 0.1 0.8];
    % set(gca, 'YTick', 1:length(label_ord), 'YTickLabel', label_tick , 'FontSize', 13);
    % ylabel('anatomy');
    % xlabel('count');
    % xlim([0 800]);
    % box off;
    % legend({'non-responsive', 'speech responsive'});

end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% S3 - Comparison of HFA peak in Spanish and English

figure;
ctr = 1;
titles = {'Spanish speaker', 'English speaker', 'Mandarin speaker'};

% Iterate over language types
for ls = [1, 2, 3]
    idx = sent_encoding.ls == ls;
    peakresp = sent_encoding.maxresp;

    % Create subplot
    subplot(2, 4, ctr)

    % Marker type for selected from single language or both language
    type = all(sent_encoding.type,2);
    en_type = sent_encoding.type(:, 2) & ~type; % only English
    sp_type = sent_encoding.type(:, 1) & ~type; % only Spanish or Mandarin

    if ismember(ls, [1, 3]) % Spanish or Mandarin, so English is foreign
        cols = [ 1 0 0; 0 0 1;];
    else % English, so Spanish or Mandarin is foreign
        cols = [ 0 0 1; 1 0 0;];
    end

    % Display the proportion of single language vs. both language
    disp(['Language type: ' num2str(ls)]);
    disp(['Both language: ' num2str(sum(idx & type)) ', ' ...
        num2str((sum(idx & type)/sum(idx))*100) '%']);
    disp(['Single languages: ' num2str(sum(idx & ~type)) ', ' ...
        num2str((sum(idx & ~type)/sum(idx))*100) '%']);
    
    % Scatter plot of peak HFA values
    scatter3(peakresp(idx & type, 2), peakresp(idx & type, 1), find(idx&type), ...
        25, [0.6 0.0 0.6], 'MarkerFaceAlpha', 0.8); hold on;

    % plot single language (English)
    scatter3(peakresp(idx & en_type, 2), peakresp(idx & en_type, 1), find(idx&en_type), ...
        25, cols(1, :),  'MarkerFaceAlpha', 0.6); hold on;

    % plot single language (Spanish or Mandarin)
    scatter3(peakresp(idx & sp_type, 2), peakresp(idx & sp_type, 1), find(idx&sp_type), ...
        25, cols(2, :), 'MarkerFaceAlpha', 0.6); hold on;
    view(2);

    % Calculate correlation coefficient
    [rho, pval] = corr(peakresp(idx, 2), peakresp(idx, 1), 'rows', 'pairwise');
    text(0.5, 4, ['r=' num2str(rho, 2) ',' getSigStr(pval, 1)], 'FontSize', 15);
    disp(['r(' num2str(sum(idx)) ')=' num2str(rho, 2) ', ' getSigStr(pval, 2)]);

    ylabel({'English peak', 'HFA (z)'});
    xlabel({'Spanish peak', 'HFA (z)'});
    set(gca, 'FontSize', 13);
    ylim([0 4.5]);
    xlim([0 4.5]);
    yticks(0:2:4);
    xticks(0:2:4);
    title(titles{ctr}, 'FontWeight', 'normal');

    % Reference line
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    h.HandleVisibility = "off";

    % Show count of native, foreign, both counts
    subplot(2, 4, ctr+4);
    if ls == 2
        h=bar([sum(idx & type), sum(idx & en_type), ...
            sum(idx & sp_type)]'./sum(idx), 'FaceColor', 'flat', ...
            'EdgeColor', 'none');
    else
        h=bar([sum(idx & type), sum(idx & sp_type), ...
            sum(idx & en_type)]'./sum(idx), 'FaceColor', 'flat', ...
            'EdgeColor', 'none');
    end
    cols = [0.6 0.0 0.6; 0 0 1; 1 0 0;]';
    h.CData = cols';
    xticks(1:3);
    xticklabels({'Both', 'Native', 'Foreign'});
    ylabel('Proportion');
    ylim([0 1]);
    yticks([0 1]);
    set(gca, 'FontSize', 13);

    %axis off;
    ctr = ctr + 1;
    box off;
end

subplot(2, 4, [4 8]);
% show blox plot of both selected, single selected

both_peak = mean(sent_encoding.maxresp(all(sent_encoding.type,2), :), 2, 'omitnan');
npidx = (~sent_encoding.type(:, 2) & ismember(sent_encoding.ls, [1, 3])) | ...
    (~sent_encoding.type(:, 1) & sent_encoding.ls==2);
native_peak = mean(sent_encoding.maxresp(npidx, :), 2, 'omitnan');
foreign_peak = mean(sent_encoding.maxresp(~npidx & ~all(sent_encoding.type,2), :), 2, 'omitnan');

% prep for boxplot
allpeak = [both_peak; native_peak; foreign_peak];
group = [ones(size(both_peak)); 2*ones(size(native_peak)); 3*ones(size(foreign_peak))];
boxplot(allpeak, group, 'Colors', [0.6 0.0 0.6; 0 0 1; 1 0 0], 'Symbol', 'o');
yticks(0:2:4);
ylabel('HFA peak (z)');
% remove bold
title('Across participants', 'FontWeight', 'normal');
xticklabels({'Both', 'Native', 'Foreign'});
set(gca, 'FontSize', 13);
box off;

% significance tests
disp('Significance tests');
disp('Both vs. Native');
disp(ranksum(both_peak, native_peak));
disp(['n = (' num2str(length(both_peak)) ', ' num2str(length(native_peak)) ')']);
disp('Both vs. Foreign');
disp(ranksum(both_peak, foreign_peak));
disp(['n = (' num2str(length(both_peak)) ', ' num2str(length(foreign_peak)) ')']);
disp('Native vs. Foreign');
disp(ranksum(native_peak, foreign_peak));
disp(['n = (' num2str(length(native_peak)) ', ' num2str(length(foreign_peak)) ')']);

clearvars -except *all subj *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% TRF: S4 - Cross-language transfer of TRF weights
% Takes a minute or two to run.

% Repeated sentences in DIMEx and TIMIT
repsentName = {'s00104','s00804', 's01904', 's03004', 's05004', ...
    's06104', 's06804', 's07404', 's09004', ...
    'fcaj0_si1479', 'fcaj0_si1804', 'fdfb0_si1948', 'fdxw0_si2141', ...
    'fisb0_si2209', 'mbbr0_si2315', 'mdlc2_si2244', 'mdls0_si998', ...
    'mjdh0_si1984', 'mjmm0_si625'};
figure;

% Subjects with no repeated sentences for DIMEx so unable to calculate
% cross-language predictions / Mandarin participants
SIDs(ismember(SIDs, {'EC252', 'EC152', 'HS8', 'HS9', 'HS10'})) = [];
for samelang = {'dimex', 'timit'}
    
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
    dataf = 100;
    for s = SIDs
        SID = s{1};
        idx = strcmp(sent_encoding.SID, SID);
        modelname = 'onset_phnfeatConsOnset_maxDtL_formantMedOnset';
    
        % iterate over languages
        testlangs = {samelang, crosslang};
        for i = 1:2
            testlang = testlangs{i};
            [out, ~] = out_addStrfpred(SID, samelang{1}, modelname, 1, ...
                sentdet, testlang); % test language will either be the same or cross
    
            % calculate the correlation between predResp and pred for each sentence
            % with the padding removed!
            tmp = nan(sum(idx), length(out));
            for j = 1:length(out)
                bef = out(j).befaft(1)*dataf;
                % if the output is 3D, then we need to mean the over the repetitions (3rd dim)
                if length(size(out(j).resp)) == 3
                    x = mean(out(j).resp(sent_encoding.el(idx), bef:end-bef, :), 3);
                else
                    x = out(j).resp(sent_encoding.el(idx), bef:end-bef);
                end
                y = out(j).predResp(sent_encoding.el(idx), bef:end-bef);
                tmp(:, j) = diag(corr(x', y', 'rows', 'pairwise'));
            end
            if i == 1
                testR_same(idx) = mean(tmp, 2);
                out_same = out;
            else
                testR_cross(idx) = mean(tmp, 2);
                out_cross = out;
            end
            clear out tmp x y;
        end
       
        % Plot several example sentence predictions
        if strcmp(SID, 'EC100')
            % [~, maxels] = max(mean(cat(3, testR_cross(1:minel, :), 
            % testR_same(1:minel, :)), [2, 3],'omitnan'));
    
            % Use the same electrodes from above for cross-prediction example
            maxels = 22; %  150
            
            % Plot each electrode example sentences as a new figure
            for maxel = maxels
                figure('renderer', 'painters');
    
                for sent = 2          
                    subplot(1, 1, 1);                            
    
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
                    ylabel('HFA (z)');
                    
                    % Plot predicted responses
                    yyaxis right
                    samepred = squeeze(out_same(sent).predResp(maxel, :));
                    crosspred = squeeze(out_cross(sent).predResp(maxel, :));       
    
                    plot(x, samepred, 'LineStyle', '--', ...
                        'LineWidth', 2, 'Color', predcols{1}, ...
                        'DisplayName', 'same-prediction');
                    plot(x, crosspred, 'LineStyle', '--', ...
                        'LineWidth', 2, 'Color', predcols{2}, ...
                        'DisplayName', 'cross-prediction');
    
                    xlim([0.2 1.5]);
                    xlabel('Time (s)')
    
                    % show predicted response correlations
                    if ~strcmp(samelang, 'dimex')
                        [same_r, ~] = corr(samepred', mean(data, 2));
                        [cross_r, ~] = corr(crosspred', mean(data, 2));
                    else
                        [same_r, ~] = corr(samepred', data');
                        [cross_r, ~] = corr(crosspred', data');
                    end
                    text(0, 0.8, ['Same corr :' num2str(same_r)])
                    text(0, 0.6, ['Cross corr :' num2str(cross_r)]);
    
                    yticks([]);
                    sentidx = strcmp({sentdet(:).name},  out_cross(sent).name);  
                    title(join(sentdet(sentidx).wordList, ' '));


                end                                                 
            end
        end
        ctr = ctr + 1;
    end
    sent_encoding.([samelang{1} '_trained']) = [testR_same, testR_cross]; 
end

% plot the cross- versus same- trained tested model comparison
figure;
titles = {'Spanish speech', 'English speech'};
trained_fields = {'dimex_trained', 'timit_trained'};

% native and foreign (blue and red)
cols = [0 0 1; 1 0 0];
%
figure;
for i = 1:2
    subplot(2, 1, i);
    x = sent_encoding.(trained_fields{i})(:, 1).^2; % same trained
    y = sent_encoding.(trained_fields{i})(:, 2).^2; % cross trained
    % do separate for each ls (native vs foreign)
    for ls = [1, 2]
        idx = sent_encoding.ls == ls;
        if ls == i % native speech condition
            scatter(x(idx), y(idx), 25, cols(1, :), 'filled', 'MarkerFaceAlpha', 0.3); hold on;
            diff_native = y(idx) - x(idx);
        else
            scatter(x(idx), y(idx), 25, cols(2, :), 'filled', 'MarkerFaceAlpha', 0.3); hold on;
            diff_foreign = y(idx) - x(idx);
        end
    end
    ylim([0 1]);
    xlim([0 1]);

    % calculate correlation
    [r, ~] = corr(x, y, 'rows', 'pairwise');
    text(0.5, 4, ['r=' num2str(r, 2)], 'FontSize', 15);

    % plot reference line
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    h.HandleVisibility = "off";

    set(gca, 'FontSize', 15);
    title(titles{i}, 'FontWeight', 'normal');
    ylabel('Cross-trained R^2');
    xlabel('Same-trained R^2');

    % subplot(2, 2, i+2);
    % diff_native = diff_native(~isnan(diff_native));
    % diff_foreign = diff_foreign(~isnan(diff_foreign));
    % boxplot([diff_native; diff_foreign], [ones(size(diff_native)); ...
    %     2*ones(size(diff_foreign))], ...
    %     'Colors', [0 0 1; 1 0 0], 'Symbol', 'o');
    % ylim([-0.2 0.4]);
    % [~, p] = ttest2(diff_native, diff_foreign);
    % text(1.5, 0.3, ['p=' num2str(p, 2)], 'FontSize', 15);
    % ylabel('Cross-trained R^2 - Same-trained R^2');
end
%% UNUSED - Electrode scatter across native only, foreign only, and both

anat_counts = cell(2, 1);
plt = [3, 2, 1];
colors = [0.8 0.1 0.2; 0.2 0.1 0.8; 0.6 0.1 0.6];
for h = {'lh', 'rh'}
    hemi = h{1};
    figure;
    set(gcf,'Color','w');

    if strcmp(hemi, 'lh')
        x_add = -20;
    else
        x_add = 10;
    end

    types = [0, 1, 2];
    
    for native = types
        numsid = 0;      

        if strcmp(hemi, 'lh')
            cortex = imgall.(SIDs{1}).img_mni.cortex;
        else
            cortex = imgall.(SIDs{6}).img_mni.cortex;
        end

        subplot(1, 3, plt(native+1))
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
            fsid = find(strcmp(sent_encoding.SID, sid), 1, 'first');
            if strcmpi(sent_encoding.hemi(fsid), hemi) 
                sidx = strcmp(sent_encoding.SID, sid);
                ls = sent_encoding.ls(fsid);
        
                if native == 1 % spanish for bilinguals
                    if ismember(ls, [1, 2])
                        chidx = sent_encoding.type(:, ls)==1 & ...
                            sent_encoding.type(:, mod(ls, 2)+1)==0;
                    else % bilingual case and mandarin case (where first column is mandarin)
                        chidx = sent_encoding.type(:, 2)==0 & ...
                            sent_encoding.type(:, 1)==1;
                    end                                
                elseif native == 0
                    if ismember(ls, [1, 2])
                        chidx = sent_encoding.type(:, ls)==0 & ...
                            sent_encoding.type(:, mod(ls, 2)+1)==1;
                    else % bilingual case and mandarin case
                        chidx = sent_encoding.type(:, 2)==1 & ...
                            sent_encoding.type(:, 1)==0;
                    end
                else % native == 2
                    chidx = all(sent_encoding.type')';
                end
                chidx = chidx & sidx;
        
                native_xyz = [native_xyz; cell2mat(sent_encoding.mni(chidx))]; 
                anat_xyz = [anat_xyz; sent_encoding.anatomy(chidx)];
                numsid=numsid+1;
            end
        end

        % number of subjects
        disp(['Is native: ' num2str(native)]);
        disp(['Number of subjects with anatomy labels: '  num2str(numsid)]);

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
        
        scatter3(native_xyz(:, 1)+x_add, native_xyz(:, 2), native_xyz(:, 3), ...
            25, colors(native+1, :), 'filled', 'MarkerFaceAlpha', 0.7);
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Functions

function [weights] = getTRFweights(SID, el, corpus, modelname, datapath)    
    [strf] = loadMultModelStrf(SID, modelname, corpus, datapath, 1);  
    weights = strf{1}.meanStrf(:, :, el);
end
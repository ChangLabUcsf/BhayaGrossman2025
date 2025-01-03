% Ilina Bhaya-Grossman
% 01.08.2022
% SIDs = [sSIDs, eSIDs];
out_crosscomp_startup;

% load in all beta model versions
SIDs = [{'HS8', 'HS9', 'HS10'}, sSIDs, eSIDs]; % , {'HS11', 'HS9', 'HS10'}
tps = 50:55;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Load in sentence responses

[sent_encoding] = loadSentenceResponse(SIDs, timit_details, dimex_details, datapath);

% load elecs
dimex_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
timit_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');
asccd_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat');

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
        SID = sent_encoding.SID{el};
        ls = sent_encoding.ls(el);
        if ismember(ls, [1, 2]) && isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
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
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% Determine how many minutes of speech each participant heard

timin = struct();
dimin = struct();
ascdn = struct();

% iterate through subjects in sent_encoding
for s = unique(sent_encoding.SID)'
    SID = s{1};
    % load in the event file for timit
    evntfile = fullfile(datapath, SID, 'timit', [SID '_timitall_evnt.mat']);
    load(evntfile, 'evnt');
    timin.(SID) = sum([evnt.StopTime]-[evnt.StartTime])/60;

    % load in the event file for dimex
    if startsWith(SID, 'EC')
        evntfile = fullfile(datapath, SID, 'dimex', [SID '_dimexall_evnt.mat']);
        load(evntfile, 'evnt');
        dimin.(SID) = sum([evnt.StopTime]-[evnt.StartTime])/60;
    end

    % load in the event file for asccd
    if startsWith(SID, 'HS')
        evntfile = fullfile(datapath, SID, 'asccd', [SID '_asccdall_evnt.mat']);
        load(evntfile, 'evnt');
        ascdn.(SID) = sum([evnt.StopTime]-[evnt.StartTime])/60;
    end
end

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


%% Display hemisphere and number of participants in each group

% load in all img data
imgall = load_allimgdata;

% display hemisphere and number of participants in each group
disp('Hemisphere and number of participants in each group');
labels = {'English', 'Spanish', 'Mandarin'};
ctr = 1;
for grp = {eSIDs, sSIDs, {'HS8', 'HS9', 'HS10'}}
    lh = 0;
    rh = 0;
    lhsids = {};
    rhsids = {};
    for sid = grp{1}
        SID = sid{1};
        if strcmpi(imgall.(SID).hemi, 'lh')
            lhsids = [lhsids, {SID}];
            lh = lh+1;
        else
            rhsids = [rhsids, {SID}]; 
            rh = rh+1;
        end
    end

    % number of participants in each group
    disp([labels{ctr} ' monolinguals : ' num2str(length(grp{1}))]);
    % number of participants in each hemisphere
    lhsids = join(lhsids, ' ');
    rhsids = join(rhsids, ' ');
    disp(['LH: ' num2str(lh) ' (' lhsids{:} ')']);
    disp(['RH: ' num2str(rh) ' (' rhsids{:} ')']);
    disp('--------------------------------------------------------');
    ctr = ctr+1;
end



%% B/C - Single example subject response magnitudes, example electrode

%'EC100', 'EC252', 'EC152', 'EC212', 'EC235', 'EC129', 'EC159', 'EC196'
nativeSIDs = {'EC214'};
nativeSIDs = {'EC100'};

% Color by language
cm = [0 0 1; 1 0 0];
modelname={'onset_phnfeatConsOnset_maxDtL_formantMedOnset'}; 
modelfeatures  = [{'onset'}; timit_details.features.names([1:3, 8, 9, 11]); ... 
    {'peakrate'; 'F1'; 'F2'; 'F3'; 'F4'}];

bins = 15;
% elecs = 203; % for EC183 , 156, 160
%elecs = 137; % for EC172
elecs = [22]; % for EC100 % [3, 54]; % for EC214
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
    plotNativeElec(nativeSIDs, desel, 1);
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

        [r, ~] = corr(mean(en_resp, 3, 'omitnan')', ...
            mean(sp_resp, 3, 'omitnan')');
        sgtitle([SID ': ' num2str(sent_encoding.el(el)) ...
            ', ' num2str(el) ', meancorr = ' num2str(r)])

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

%% E - Same as above but single scatter

imgall = load_allimgdata;
hemi = 'lh';

figure;
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
    else
        cortex = imgall.(SIDs{6}).img_mni.cortex;
    end

    subplot(1, 2, native+1)
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
                else
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
                numsid=numsid+1;
            end
        end
    
    end

    % number of subjects
    disp(['Is native: ' num2str(native)]);
    disp(['Number of subjects: '  num2str(numsid)]);
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
    if ~native
        color = [1, 0, 0];
    else
        color = [0, 0, 1];
    end

    scatter3(native_xyz(:, 1)+x_add, native_xyz(:, 2), native_xyz(:, 3), ...
        round(native_hga*20), color, 'filled', 'MarkerFaceAlpha', 0.7);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% F - Example electrode TRF weights

% example electrodes in the figure
% SID = 'EC172';
% els = [74,  136];
SID = 'EC100';
els = [70 71];

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
    sgtitle([num2str(el) ' ' num2str(corr(weights{2, 1}(:), weights{2, 2}(:)))])
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
            sent_encoding.sp_base_rsq(e)>thresh) && ~missing  
    
            % add in TRF weight correlation
            [~, maxtp] = max(squeeze(mean(meanStrf(:, :, :, el), [1 2])));
            wind = max(1, maxtp-windsz):min(maxtp+windsz, size(meanStrf, 3));
        
            % compare frequency selectivity
            x = squeeze(mean(meanStrf(1, :, wind, el), 3));
            y = squeeze(mean(meanStrf(2, :, wind, el), 3));
            
            % scale and smooth the STRF weights to be more comparable
            [~, z] = procrustes(x', y');
            if strcmp(modelnames{1}, 'onset_aud')
                z = smoothdata(z, 'SmoothingFactor',0.4);
                x = smoothdata(x, 'SmoothingFactor',0.1);
            end     
            
            % save out weights
            alleng(e, 1:length(x)) = x;
            allsp(e, 1:length(z)) = z;

            % calculate correlation and max frequency
            [corrstrf(e), corrpval(e)] = corr(x', y', 'type', 'Pearson');
            % Spearman withou procrustes

            for j = 1:10
                [corrstrfperm(e, j), ~] = corr(x', y(randperm(length(y)))', 'type', 'Pearson');
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
legend({'r<=0.5', 'r>0.5'});

for ls = [1, 2, 3]
    subplot(2, 3, ls+3)
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

disp(['Correlation between beta matrices: ' num2str(corr(allsp(:), alleng(:)))]);

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
%
% show histogram of correlations
figure;
hbins = -1:0.05:1;

% permutation distribution
idx = minrsq>0.1 & ~nanrows;
corrstrfperm(~idx, :) = NaN;
histogram(corrstrfperm(:), hbins, 'EdgeColor', 'none', ...
    'FaceColor', [0.7, 0.7, 0.7], 'FaceAlpha', 0.5, 'Normalization','probability');
hold on;
% make a kdensity plot
yyaxis right; 
[f, xi] = ksdensity(corrstrfperm(:));
plot(xi, f, 'LineWidth', 2, 'Color', [0.7, 0.7, 0.7]);

yyaxis left; 
histogram(corrstrf(idx), hbins, 'EdgeColor', 'none', ...
    'FaceColor', [0.7, 0.1, 0.7], 'FaceAlpha', 0.5, 'Normalization', 'probability');
hold on;
[f, xi] = ksdensity(corrstrf(idx));
yyaxis right; 
plot(xi, f, 'LineWidth', 2, 'Color', [0.7, 0.1, 0.7]);
xline(prctile(corrstrfperm(:), 95), 'LineWidth', 2, 'Color', 'k');
xlim([-1 1]);

yticks([]);
yyaxis left;
ylabel('Normalized Density');
box off;

disp(['Total number of electrodes: ' num2str(sum(idx))])
disp(['Number of electrodes with higher than 95% for permuted distrib: ' ...
    num2str(sum(corrstrf(idx)>prctile(corrstrfperm(:), 95)))]);


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd* corrstrf corrpval maxrsq;

%% Visualizing feature selectivity across subjects

% colors = [0.45 0.67 0.18; ... % onset
%           0.7 0.3 0.7; ... % peakrate
%           0.8 0.5 0.7; ... % phonetic feature
%           0 0.44 0.74; ... % relative pitch
%           0.8 0.35 0.15; ... % absolute pitch
%           0.92 0.69 0.12];    % both
figure;

colors = [0.45 0.67 0.18; ... % onset
          0.8 0.35 0.15; ... % consonant
          0.1 0 0; ... % relative pitch (UNUSED)
          0.1 0.7 0.74; ... % peakrate
          0.92 0.69 0.12];    % vowel

percent = nan(3, 5, ...
    length(unique(sent_encoding.SID(sent_encoding.ls==ls))));

for ls = [1, 2, 3]
    
    ctr = 1;
    for s = unique(sent_encoding.SID(sent_encoding.ls==ls))'
        SID = s{1};
        idx = strcmp(sent_encoding.SID, SID);
        
        % get proportion of electrodes for each primary encoding
        count = histcounts(sent_encoding.primary(idx), 0.5:5.5);
        percent(ls, :, ctr) = (count./sum(count))*100;

        ctr = ctr + 1;
    end
end

% plot bar graph with error bars for this group of subjects where all features (dim 2)
% next to each other
for ls = [1, 2, 3]
    idx = [1:2 4:5];
    ctr=1;
    for i = idx
        x = ctr*8.5 + (ls-1)*2 + 1;
        % only include in legend if it's the first language
        if i==idx(1)
            h = bar(x, mean(percent(ls, i, :), 3, 'omitnan'), 1.2); hold on;
        else
            h = bar(x, mean(percent(ls, i, :), 3, 'omitnan'), 1.2, 'HandleVisibility', 'off'); hold on;
        end
    
        % have each bar be a different color
        h.FaceColor = 'flat';
        h.CData = colors(i, :); % (i, :)

        if ls == 1
            % remove line color
            h.EdgeColor = 'none';        
        elseif ls == 3
            h.EdgeColor = colors(i, :);
            h.LineWidth = 2; 
            h.FaceAlpha = 0;
        else
            h.FaceAlpha = 0.5;
            h.EdgeColor = colors(i, :);
            h.LineWidth = 2; 
        end
        ctr=ctr+1;
    end

    x = (1:length(idx))*8.5 + (ls-1)*2 + 1;
    err = nansem(percent(ls, idx, :), 3);
    er = errorbar(x,mean(percent(ls, idx, :), 3, 'omitnan'),err);    
    er.Color = [0 0 0];                            
    er.LineStyle = 'none';
    er.LineWidth = 1.5;
    er.CapSize = 0;
    er.HandleVisibility = 'off';
end

% % plot bar graph with error bars for this group of subjects
% idx = [1:2 4:5];
% x = (1:length(idx)) + (ls-2)*7;
% h = bar(x, mean(percent(idx, :), 2, 'omitnan')); hold on;
% 
% % have each bar be a different color
% for i = 1:length(idx)
%     h.FaceColor = 'flat';
%     h.CData(i, :) = colors(idx(i), :);
% 
%     % remove line color
%     h.EdgeColor = 'none';
% end
% err = nansem(percent(idx, :), 2);
% er = errorbar(x,mean(percent(idx, :), 2, 'omitnan'),err);    
% er.Color = [0 0 0];                            
% er.LineStyle = 'none';
% er.CapSize = 0;
% er.HandleVisibility = 'off';

box off;
xlabels = {'onset', 'consonant', 'peak rate', 'vowel formant'}; % 'relative/absolute pitch',
llabels = {'Spanish', 'English', 'Mandarin'};
legend(llabels);
ylabel('Proportion of electrodes (%)');
xticks([11.5, 20, 28.5, 37]);
xticklabels(xlabels);
set(gca, 'FontSize', 13);
ylim([0 40]);
yticks(0:20:100);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% ----------------------- Supplementary Figures --------------------------
%% S0 - Proficiency metrics

varnames = {'SID', 'ls', 'spanish', 'english', 'mandarin', 'other'};
proficiency =  array2table(zeros(0, 6), 'VariableNames', varnames);
titles = {'Spanish speakers', 'English speakers', 'Mandarin speakers'};

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
mSIDs={'HS10', 'HS9', 'HS8'}'; % 'HS14',
mprof = {5; 5; 5}; % 5
eprof = {0; 0; 0}; % 2 other language proficiency
sprof = {0; 0; 0}; % 0
other = {0; 0; 0}; % 0
lang = repmat({'m'}, length(mSIDs), 1);
proficiency = [proficiency; [mSIDs, lang, sprof, eprof, mprof, other] ];

% Plotting proficiencies
colors = brewermap(3, 'Dark2');
lss = {'s', 'e', 'm'};
figure('Position',[500, 80, 200, 500])
for i = 1:3
    subplot(2, 2, i);
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
        65, [0.5, 0.5, 0.5], 'filled', 'MarkerFaceAlpha', 0.7); hold on;
    
    % colors(i, :)
    plot(x', y', 'LineWidth', 0.5, ...
        'Color', [0.5, 0.5, 0.5], 'HandleVisibility', 'off');

    % Formatting 
    xlim([-0.5 3.5])
    ylim([-0.5 5.5]);
    ylabel('Proficiency');
    xticks([]);
    
    text(3, 4.5, ['n=' num2str(numid)]);
    set(gca, 'FontSize', 13)
    title(titles{i}, 'Fontweight', 'normal');

    xticks(0:3);
    xticklabels({'Span', 'Eng', 'Mand', 'Other'});
    xtickangle(0);
    xlabel('Language');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% S1 - Spectral and Temporal MTFs for all stimuli
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

imgall = load_allimgdata;

for h = {'lh', 'rh'}
    hemi = h{1};
    figure;
    set(gcf,'Color','w');
    axes;
    hold on
    
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
    
    subplot(1, 2, 1);
    PlotBrainSurface(cortex, hemi,'lateral');
    alpha 0.9
    light("Style","infinite","Position",[100 100 0]);
    
    % find density map of native speech, lateral side
    native_xyz = [];
    native_hga = [];
    all_xyz = [];
    anatomy = [];
    
    for si = SIDs
        sid = si{1};
        if strcmpi(imgall.(sid).hemi,hemi)
            subject = imgall.(sid);
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix; 
                anatomy = [anatomy; subject.img_mni.anatomy];  
                all_xyz = [all_xyz; elecmatrix]; 
                numsid=numsid+1;
            end
        end
    
    end
    
    % number of subjects
    disp(['Number of subjects: '  num2str(numsid)]);
    disp(['Number of electrodes: '  num2str(length(all_xyz))]);
    disp('--------------------------------------------------------');
    scatter3(all_xyz(:, 1)+x_add, all_xyz(:, 2), all_xyz(:, 3), 7, 'k', 'filled');
    
    % show distributions of anatomical coverage
    subplot(1, 2, 2);
    set(gcf,'Color','w');
    hold on
    [counts, ~, ~, labels] = crosstab(anatomy(:, 4));
    % remove all cases where coveage < 20
    mincount = 20;
    labels = labels(counts>mincount);
    counts = counts(counts>mincount);
    
    % remove NaN and Unknown, Right-Cerebral-WhiteMatter
    torem = {'NaN', 'Unknown', 'Right-Cerebral-White-Matter', ...
                'Left-Cerebral-White-Matter', 'Left-Hippocampus', ...
                'Right-Hippocampus'};
    counts = counts(~ismember(labels, torem));
    labels = labels(~ismember(labels, torem));
    
    % order alphabetically
    [labels, idx] = sort(labels);
    counts = counts(idx);
    
    barh(counts, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    set(gca, 'YTick', 1:length(labels), 'YTickLabel', labels, 'FontSize', 13);
    ylabel('anatomy');
    xlabel('count');
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
    text(0.5, 4, ['r=' num2str(rho, 2)], 'FontSize', 15);
    text(0.5, 3.5, ['p=' num2str(pval, 2)], 'FontSize', 15);

    ylabel({'English peak', 'HFA (z)'});
    xlabel({'Spanish peak', 'HFA (z)'});
    set(gca, 'FontSize', 13);
    ylim([0 4.5]);
    xlim([0 4.5]);
    title(titles{ctr}, 'FontWeight', 'normal');

    % Reference line
    h = refline(1, 0);
    h.LineWidth = 2;
    h.Color = 'k';
    h.HandleVisibility = "off";
    %legend({'Both languages', 'Single language'}, 'Location', 'northwest');

    % Show count of native, foreign, both counts
    subplot(2, 4, ctr+4);
    % histogram(peakresp(idx, 2)-peakresp(idx, 1), EdgeColor='none', ...
    %     FaceColor = [0.4, 0.4, 0.4]); 
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
    % for i = 1:3
    %     h(i).FaceColor = cols(i, :);
    % end
    % xline(0);
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
yticks(0:4);
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
disp('Both vs. Foreign');
disp(ranksum(both_peak, foreign_peak));
disp('Native vs. Foreign');
disp(ranksum(native_peak, foreign_peak));

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% S4 - Cross-language transfer of TRF weights
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

%% Info about participant responses

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

%% ----------------------------- UNUSED Panels ----------------------------
%% Discretizing phonetic feature selectivity 

allidx_timit = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat', 'allidx');
allidx_dimex = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat', 'allidx');
allidx_asccd = load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat', 'allidx');

% load in TRF models
% with pitch change feature
modelnames={'onset_phnfeatonset_maxDtL', ...
    'phnfeatonset_maxDtL', ... %2
    'onset_phnfeatonset', ...  %3
    'onset_maxDtL', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0Bin_relPitchBin_F0ChangeBin', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_F0Bin', ...
    'onset_maxDtL_maxDtLOnset_vowelOnset_relPitchBin_F0ChangeBin', ...
    'aud', ...
    'onset', ...
    'onset_maxDtL_vowelOnset_aud', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %11
    'onset_maxDtL_formantMedOnset', ... %12
    'onset_phnfeatConsOnset_maxDtL'}; %13

% determine unique variance per feature and primary encoding
varnames = {'SID', 'el', 'ls', 'base_rsq', 'spect_rsq', 'basepitch_rsq', ...
    'uv_onset', 'uv_pr', 'uv_phn', 'uv_f0', 'uv_relf0', 'uv_formant', 'uv_consonant'};
encodings =  array2table(zeros(0,length(varnames)), 'VariableNames', varnames);

for s = SIDs
    SID = s{1}; 
    corpusStrf=loadMultModelStrf(SID, modelnames, 'timit', datapath, 1);

    % check whether any of the models are empty
    if any(cellfun(@isempty, corpusStrf))
        continue;
    end
    
    % models without pitch
    base = (corpusStrf{1}.meanTestR.^2)';
    spect = (corpusStrf{8}.meanTestR.^2)';
    onset = (corpusStrf{9}.meanTestR.^2)';

    % from Yulia's paper
    uvOnset = onset - spect;
    uvPhnfeat = base - spect; 
    
    % models with pitch
    base_wtpitch = (corpusStrf{5}.meanTestR.^2)';    
    uvRelF0 = base_wtpitch - corpusStrf{6}.meanTestR.^2';
    uvF0 = base_wtpitch - corpusStrf{7}.meanTestR.^2';

    % formant and consonant models
    uvFormant = corpusStrf{11}.meanTestR.^2' - corpusStrf{13}.meanTestR.^2';
    uvConsonant = corpusStrf{11}.meanTestR.^2' - corpusStrf{12}.meanTestR.^2';
    uvPeakrate = corpusStrf{1}.meanTestR.^2' - corpusStrf{3}.meanTestR.^2';
    
    % find which electrodes meet threshold rsq for the basemodel
    tbl_els = sent_encoding.el(strcmp(sent_encoding.SID, SID));
    tidx = find(strcmp(sent_encoding.SID, SID));

    % add to table 
    sent_encoding.base_rsq(tidx) = base(tbl_els);
    sent_encoding.spect_rsq(tidx) = spect(tbl_els);
    sent_encoding.basepitch_rsq(tidx) = base_wtpitch(tbl_els);
    sent_encoding.uv_onset(tidx) = uvOnset(tbl_els);
    sent_encoding.uv_phn(tidx) = uvPhnfeat(tbl_els);
    sent_encoding.uv_f0(tidx) = uvF0(tbl_els);
    sent_encoding.uv_relf0(tidx) = uvRelF0(tbl_els);
    sent_encoding.uv_formant(tidx) = uvFormant(tbl_els);
    sent_encoding.uv_consonant(tidx) = uvConsonant(tbl_els);
    sent_encoding.uv_pr(tidx) = uvPeakrate(tbl_els);
end

% at least one unique variance > 0.01 as a threshold
% initial preferences, primary encoding
pref = zeros(height(sent_encoding), 1);
thresh = 0.005;

% onsets
pref(sent_encoding.uv_onset>0) = 1;
% peakrate and features
%pref(encodings.uv_onset<0 & encodings.uv_phn>0) = 3;
pref(sent_encoding.uv_onset<0 & sent_encoding.uv_pr>0)=2;

% unique variance for pitch features calculated differently
% if uv is greater than 0, primary feature is pitch (both or rel or abs)
% pref(encodings.uv_f0 < 0 & encodings.uv_relf0 > thresh) = 4; % relative
% pref(encodings.uv_f0 > thresh & encodings.uv_relf0 < 0) = 5; % absolute
% pref(encodings.uv_f0 > thresh & encodings.uv_relf0 > thresh) = 6; % both
pref(sent_encoding.uv_f0 > thresh | sent_encoding.uv_relf0 > thresh)=3;

% if consonant feature is greater than the threshold
pref(sent_encoding.uv_formant > thresh & ...
    (sent_encoding.uv_formant>sent_encoding.uv_f0 | ...
    sent_encoding.uv_formant>sent_encoding.uv_relf0)) = 5;
pref(sent_encoding.uv_consonant>thresh) = 4; % consonant

% 1 - onset, 2- peakrate, 3 - pitch, 4 - consonant, 5  - formant
% pref(~any(uvs > 0.02)) = 0;
sent_encoding.primary = pref;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;

%% Proportion of speech responsive elecs for each subject group

% electrode selection details
timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");

figure;
for ls = [1, 2, 3]
    
    % find number of speech responsive electrodes in table
    idx = sent_encoding.ls == ls;
    sids = unique(sent_encoding.SID(idx));
    numspeech = zeros(length(sids), 1);
    numel = zeros(length(sids), 1);

    % for each subject, identify the number of total elecs
    ctr = 1;
    for s = sids'
        sid = s{1};
        numspeech(ctr) = sum(idx & strcmp(sent_encoding.SID, sid));
        numel(ctr) = length(timit_elecs.fvals.(sid));
        ctr = ctr + 1;
    end

    % plot the pie chart for each language group
    subplot(1, 3, ls);
    % explode and show percentage
    pc = pie([sum(numspeech), sum(numel)-sum(numspeech)], ...
        [0 1]);
    % no edge color
    pc(1).FaceColor = [0.7 0.0 0.7];
    pc(3).FaceColor = [0.5 0.5 0.5];
    pc(1).EdgeColor = 'none';
    pc(3).EdgeColor = 'none';
    if ls == 1
        legend({'Speech responsive', 'Non-responsive'}, 'Location', 'southoutside');
    end
end


%% Native vs. unfamiliar contour plot

figure;
set(gcf,'Color','w');

ha = axes;
hemi = 'rh';
hold on

imgall = load_allimgdata;
cmax = 0.3;

types = [0, 1];
anatomy = cell(1, 2);

for native = types
    numsid = 0;      

    if strcmp(hemi, 'lh')
        cortex = imgall.(SIDs{1}).img_mni.cortex;
    else
        cortex = imgall.(SIDs{6}).img_mni.cortex;
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
                        native_hga = [native_hga; sent_encoding.maxresp(intersect(ch_sel, ch_sid), 1)];
                    elseif ls==2
                        native_hga = [native_hga; sent_encoding.maxresp(intersect(ch_sel, ch_sid), 2)];
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
    
        %         ch_use = setdiff(ch_use,ch_temporal);
                native_xyz = [native_xyz; elecmatrix(intersect(ch_sel, ch_sid),:)]; 
                anatomy{native+1} = [anatomy{native+1}; subject.img_mni.anatomy(intersect(ch_sel, ch_sid),4)];
                all_xyz = [all_xyz; elecmatrix]; 
                numsid=numsid+1;
            end
        end
    
    end

    % number of subjects
    disp(['Number of subjects: '  num2str(numsid)]);
    
    yye = min(native_xyz(:,2))-15:1:max(native_xyz(:,2)+15);
    zze = min(native_xyz(:,3))-15:1:max(native_xyz(:,3)+15); 
    ds = histcounts2(native_xyz(:,2),native_xyz(:,3),yye,zze);
    ds_all = histcounts2(all_xyz(:,2),all_xyz(:,3),yye,zze);
    ds_norm = ds./ds_all;
    ds_norm(isnan(ds_norm)) = 0;
    
%     for weighting by high gamma
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
    clim([0 cmax])
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

    figure;
    if ~native
        colormap(flipud(reds));
    else
        colormap(flipud(blues));
    end
    clim([0 cmax]);
    axis off;
    colorbar();
end

% figure show native vs. nonnative anatomy in barh
figure;
set(gcf,'Color','w');
mincount = 10;
axes;
counts_both = [];
for native = types
    hold on
    [counts, ~, ~, labels] = crosstab(anatomy{native+1});
   
    labels = labels(counts>mincount);
    counts = counts(counts>mincount);

    % make labels anatomical
    torem = {'NaN', 'Unknown', 'Right-Cerebral-WhiteMatter', ...
            'Left-Cerebral-White-Matter', 'Left-Hippocampus', ...
            'Right-Hippocampus'};
    counts = counts(~ismember(labels, torem));
    labels = labels(~ismember(labels, torem));
    [labels, idx] = sort(labels);
    counts = counts(idx);

    counts_both = [counts_both counts];
end

% make native and nonnative bars next to each other, set colors
barh(counts_both, 'EdgeColor', 'none');
set(gca, 'YTick', 1:length(labels), 'YTickLabel', labels, 'FontSize', 13);
axis off;

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% All coverage as contour plot

for h = {'lh', 'rh'}
    figure;
    set(gcf,'Color','w');
    
    ha = axes;
    hemi = h{1};
    hold on
    
    imgall = load_allimgdata;
    cmax = 2;
    numsid = 0;      
    
    if strcmp(hemi, 'lh')
        cortex = imgall.(SIDs{1}).img_mni.cortex;
    else
        cortex = imgall.(SIDs{6}).img_mni.cortex;
    end
    
    PlotBrainSurface(cortex, hemi,'lateral');
    alpha 0.9
    light("Style","infinite","Position",[100 100 0]);
    
    ax_pos = get(ha,'Position');
    yy_range = get(ha,'YLim');
    zz_range = get(ha,'ZLim');
    
    % find density map of native speech, lateral side
    % native_xyz = [];
    % native_hga = [];
    all_xyz = [];
    anatomy = [];
    
    for si = SIDs
        sid = si{1};
        if strcmpi(imgall.(sid).hemi,hemi)
            subject = imgall.(sid);
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix;   
                anatomy = [anatomy; subject.img_mni.anatomy];          
    
                all_xyz = [all_xyz; elecmatrix]; 
                numsid=numsid+1;
            end
        end
    
    end
    
    % number of subjects
    disp(['Number of subjects: '  num2str(numsid)]);
    disp(['Number of electrodes: '  num2str(length(all_xyz))]);
    disp('--------------------------------------------------------');
    
    % find density map of native speech, lateral side
    yye = min(all_xyz(:,2))-15:1:max(all_xyz(:,2)+15);
    zze = min(all_xyz(:,3))-15:1:max(all_xyz(:,3)+15); 
    ds_all = histcounts2(all_xyz(:,2),all_xyz(:,3),yye,zze);
    ds_all(isnan(ds_all)) = 0;
    
    % smooth the density map with a gaussian kernel
    gs_kernel = fspecial('gaussian', [15, 15], 3);
    ds_sm = conv2(ds_all,gs_kernel,'same');
    
    yy = yye(1)+diff(yye(1:2))/2:diff(yye(1:2)):yye(end);
    zz = zze(1)+diff(zze(1:2))/2:diff(zze(1:2)):zze(end);
    yyq = yy(1):0.2:yy(end);
    zzq = zz(1):0.2:zz(end); 
    ds_q = interp2(zz,yy,ds_sm,zzq,yyq','cubic');
    
    ha_ct = axes('Position',ax_pos);
    hold on
    [~,hc] = contourf(yyq,zzq,ds_q',10,'LineColor','none');
    % [~,hc] = contourf(yye,zze,ds_sm',15,'LineColor','none');
    
    clim([0 cmax])
    if strcmp(hemi, 'lh')
        set(ha_ct,'Color','None','XDir','reverse');
    else
        set(ha_ct,'Color','None');
    end
    axis equal
    axis off
    xlim(yy_range);
    ylim(zz_range);
    colormap(ha_ct,internet);
    drawnow;
    
    hFills = hc.FacePrims;  % array of TriangleStrip objects
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    for i = 1:length(hFills)
        hFills(i).ColorData(4) = 180;   % default=255
    end
    hFills(1).ColorData(4) = 0;
    disp(numsid);
    
    % show distributions of anatomical coverage
    figure;
    set(gcf,'Color','w');
    ha = axes;
    hold on
    [counts, ~, ~, labels] = crosstab(anatomy(:, 4));
    % remove all cases where coveage < 20
    mincount = 20;
    labels = labels(counts>mincount);
    counts = counts(counts>mincount);
    % move colorbar to the Left
    h = colorbar;
    set(h, 'Location', 'westoutside');
    colormap(internet)
    clim([0 cmax]);
    
    % remove NaN and Unknown, Right-Cerebral-WhiteMatter
    torem = {'NaN', 'Unknown', 'Right-Cerebral-White-Matter', ...
                'Left-Cerebral-White-Matter', 'Left-Hippocampus', ...
                'Right-Hippocampus'};
    counts = counts(~ismember(labels, torem));
    labels = labels(~ismember(labels, torem));
    
    % order alphabetically
    [labels, idx] = sort(labels);
    counts = counts(idx);
    
    barh(counts, 'FaceColor', [0.5 0.5 0.5], 'EdgeColor', 'none');
    set(gca, 'YTick', 1:length(labels), 'YTickLabel', labels, 'FontSize', 13);
    ylabel('anatomy');
    xlabel('count');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;


%% Functions

function [weights] = getTRFweights(SID, el, corpus, modelname, datapath)    
    [strf] = loadMultModelStrf(SID, modelname, corpus, datapath, 1);  
    weights = strf{1}.meanStrf(:, :, el);
end
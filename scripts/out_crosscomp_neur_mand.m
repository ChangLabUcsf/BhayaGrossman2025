%% Mandarin subject validation

% Load data
%% Set up

% Ilina Bhaya-Grossman
% 03.23.2023
addpath(genpath('../../../ecog_scripts'))
addpath(genpath('../../../plotting_scripts'))
addpath(genpath('util'))

zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;

SIDs = {'HS11', 'HS8', 'HS9', 'HS10'}; % 'HS32' 
% english proficiency
prof_all = [4, 0, 0, 0];

% get color map
spec = spectral(8);
spec = flipud(spec([1:3 5:6], :));
cols = arrayfun(@(x) spec(x+1, :), prof_all, 'UniformOutput',false);
cols_all = cat(1, cols{:});

timit_details = load('out_sentence_details_timit_all_loudness.mat');
asccd_details = load('out_sentence_details_asccd_all_loudness.mat');
% tps = 50:55;

bef=50;
aft=50;

% loading in subject data
% TDwrd = loadDwrdMulti('timit', bef, aft, SIDs, timit_details);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% load in all mandarin data (resegmenting mandarin data)

corpus = 'asccd';
carflags = {'noCAR_', 'noCAR_', 'noCAR_', 'noCAR_'};
block = [27, 41, 30];

SIDs = {'EC204', 'EC222', 'EC186'};
for s = 1:length(SIDs)
    SID = SIDs{s};
    disp(['Processing subject: ' SID '...'])

    [out] = loadOutFile(datapath, corpus, SID, block(s), ''); % noCAR_

    % resegment data into shorter sentences
    for i = 1:length(out)
        % find silences (time points to segment at between sentences)
        sound = out(i).sound;
        fs = out(i).soundf;
        [offset_silences, silence_durs] = getSilenceOffsets(sound, fs, 0.05);

        % segment all the fields in out
        [segOut] = resegmentOut(out(i), offset_silences);
        if i == 1
            new_out = segOut;
        else
            new_out = [new_out; segOut];
        end
       
        clear segOut offset_silences offset_silences_dataf silence_dur;
    end

    out = new_out;
    file = fullfile(datapath, SID, corpus, 'block_z', ...
        [SID '_B' num2str(block) '_HilbAA_70to150_8band_noCAR_all_0_mel_' ...
        upper(corpus) '_zflag_global_out_resp_log_SEGMENT.mat']);

    save(file, 'out');
    
    % [outs{2}] = loadOutFile(datapath, 'timit', SID, 10, carflags{s});
end

%% get speech responsive electrodes in mandarin subjects

% english stimuli
[allidx_en, ~] = getSpeechElecs(SIDs, 'bychan', 'timit', 256, 'asccd');

% mandarin stimuli
[allidx_mn, ~] = getSpeechElecs(SIDs, 'bychan', 'asccd', 256, 'timit');
%
% make an empty table with variable names SID, el, resp_en, resp_mn
varnames = {'SID', 'el', 'resp_en', 'resp_mn'};
mand_encoding = array2table(zeros(0, length(varnames)), 'VariableNames', varnames);

for i = 1:length(SIDs)
    SID = SIDs{i};
    els = union(allidx_en.(SID), allidx_mn.(SID));

    % make a binary variable for whether electrode is speech responsive in english
    resp_en = zeros(length(els), 1);
    resp_en(ismember(els, allidx_en.(SID))) = 1;

    % make a binary variable for whether electrode is speech responsive in mandarin
    resp_mn = zeros(length(els), 1);
    resp_mn(ismember(els, allidx_mn.(SID))) = 1;

    % add to table
    tmp = table(repmat({SID}, length(els), 1), els, resp_en, resp_mn, ...
        'VariableNames', varnames); 
    mand_encoding = [mand_encoding; tmp];
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% look at example TRF weights for mandarin subjects

% load in TRF weights
modelnames = {'onset_aud','onset_phnfeatConsOnset_maxDtL_formantMedOnset'};
corpus = {'asccd', 'timit'};

% find phonetic feature similarity between english and mandarin
% consonant features

% change 'is_vowel' to 'syllabic', remove all '_vowel' suffixes
% asccd_feats{strcmp(asccd_feats, 'is_vowel')} = 'syllabic';
asccd_feats = asccd_details.features.names(1:5);
asccd_feats = strrep(asccd_feats, '_vowel', '');
timit_feats = timit_details.features.names([1:3, 8, 9, 11]);

shared_feats = intersect(timit_feats, asccd_feats);
asccd_fidx = cellfun(@(x) find(strcmp(asccd_feats, x)), shared_feats);
timit_fidx = cellfun(@(x) find(strcmp(timit_feats, x)), shared_feats);

for s = 1:length(SIDs)
    SID = SIDs{s};
    
    for corp = 1:2
        
        corpusStrf = loadMultModelStrf(SID, modelnames, corpus{corp}, ...
            datapath, 1, 'v5');
        
        % get TRF weights for speech responsive electrodes
        els = find(strcmp(mand_encoding.SID, SID));
        strf = corpusStrf{1}.meanStrf;
        trf = corpusStrf{2}.meanStrf;
        for e = els'
            if corp == 1
                mand_encoding.strf(e) = {strf(:, :, mand_encoding.el(e))};
                mand_encoding.strf_rsq(e) = corpusStrf{1}.meanTestR(mand_encoding.el(e)).^2;

                % get the shared consonant features and the
                % peakrate/formants
                fidx = [1; 1+asccd_fidx; (size(trf, 1)-4:size(trf, 1))'];
                mand_encoding.trf(e) = {trf(fidx, :, mand_encoding.el(e))};
                mand_encoding.trf_rsq(e) = corpusStrf{2}.meanTestR(mand_encoding.el(e)).^2;
            else
                mand_encoding.strf_en(e) = {strf(:, :, mand_encoding.el(e))};
                mand_encoding.strf_rsq_en(e) = corpusStrf{1}.meanTestR(mand_encoding.el(e)).^2;

                % get the shared consonant features and the
                % peakrate/formants
                fidx = [1; 1+timit_fidx; (size(trf, 1)-4:size(trf, 1))'];

                mand_encoding.trf_en(e) = {trf(fidx, :, mand_encoding.el(e))};
                mand_encoding.trf_rsq_en(e) = corpusStrf{2}.meanTestR(mand_encoding.el(e)).^2;

                % correlate strf weights for maximum window
                window = 10;
                % find maximum
                [r(1), ~] = corr(reshape(mand_encoding.strf_en{e}, [], 1), ...
                    reshape(mand_encoding.strf{e}, [], 1), 'type','Pearson');
                mand_encoding.strf_corr(e) = r(1);
                [r(2), ~] = corr(reshape(mand_encoding.trf_en{e}, [], 1), ...
                    reshape(mand_encoding.trf{e}, [], 1), 'type','Pearson');
                mand_encoding.trf_corr(e) = r(2);
                
                if r(1)>0.7
                    figure
                    ax = subplot(2, 2, 3);
                    imagesc(squeeze(mand_encoding.strf_en{e}));
                    title(sprintf('SID %s, %s, el %d', SID, corpus{corp}, e));
                    colormap(ax, flipud(reds));
                    yticks([1 80]);
                    yticklabels({'0.01', '8'});
                    ylabel('frequency (kHz)');

                    ax = subplot(2, 2, 1);
                    imagesc(squeeze(mand_encoding.strf{e}));
                    title(sprintf('SID %s, %s, el %d', SID, 'asccd', e));
                    colormap(ax, flipud(blues));
                    yticks([1 80]);
                    yticklabels({'0.01', '8'});
                    ylabel('frequency (kHz)');
                    
                    ax = subplot(2, 2, 4);
                    imagesc(squeeze(mand_encoding.trf_en{e}));
                    title(sprintf('SID %s, %s, el %d', SID, corpus{corp}, e));
                    yticks(1:length(shared_feats)+1);
                    yticklabels([{'onset'}, shared_feats' ...
                        {'peakRate', 'F1', 'F2', 'F3', 'F4'}]); % , {'peakRate'}
                    colormap(ax, flipud(reds))
                    
                    ax = subplot(2, 2, 2);
                    imagesc(squeeze(mand_encoding.trf{e}));
                    title(sprintf('SID %s, %s, el %d', SID, 'asccd', e));
                    feats = asccd_details.features.names;
                    yticks(1:length(shared_feats)+1);
                    yticklabels([{'onset'}, shared_feats' ...
                        {'peakRate', 'F1', 'F2', 'F3', 'F4'}]); % , {'peakRate'}
                    colormap(ax, flipud(blues));

                    sgtitle(['trf corr: ' num2str(r(2)) 'strf corr: ' num2str(r(1))])

                    for k = 1:4
                        subplot(2, 2, k)
                        xlim([0, 40]);
                        xticks([0 , 40]);
                        xticklabels({'0', '-0.4'});
                        set(gca, 'YDir', 'normal', 'fontsize', 13);

                    end
                end
            end
        end
    end
end

% show pie for correlation threshold
idx = min(mand_encoding.strf_rsq_en, mand_encoding.strf_rsq) >0.1;
corrthresh = 0.5;
corrstrf = mand_encoding.strf_corr;

figure;
subplot(2, 1, 1);
labels = num2str(crosstab(corrstrf(idx)>corrthresh));
p = pie(crosstab(corrstrf(idx)>corrthresh), [1, 1], labels);
p(3).FaceColor = [159, 134, 192]/256;
p(3).EdgeColor = 'none';
p(1).FaceColor = [224, 177, 203]/256;
p(1).EdgeColor = 'none';
p(2).FontSize = 13;
p(4).FontSize = 13;
legend({'r<=0.5', 'r>0.5'});

idx = min(mand_encoding.trf_rsq_en, mand_encoding.trf_rsq) >0.1;
corrthresh = 0.5;
corrstrf = mand_encoding.trf_corr;

subplot(2, 1, 2);
labels = num2str(crosstab(corrstrf(idx)>corrthresh));
p = pie(crosstab(corrstrf(idx)>corrthresh), [1, 1], labels);
p(3).FaceColor = [159, 134, 192]/256;
p(3).EdgeColor = 'none';
p(1).FaceColor = [224, 177, 203]/256;
p(1).EdgeColor = 'none';
p(2).FontSize = 13;
p(4).FontSize = 13;
legend({'r<=0.5', 'r>0.5'});


clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% show contour of active electrodes in mandarin subjects

% Can't do this because no warped electrodes for these subjects
figure;
set(gcf,'Color','w');

ha = axes;
hemi = 'lh';
hold on

imgall = load_allimgdata;

types = [1];
numsid = 0;      

for native = types

    if strcmp(hemi, 'lh')
        cortex = imgall.(SIDs{1}).img_mni.cortex;
    else
        cortex = imgall.(SIDs{3}).img_mni.cortex;
    end

    PlotBrainSurface(cortex, hemi,'lateral');
    alpha 0.9
    %light("Style","infinite","Position",[100 100 0]);

    ax_pos = get(ha,'Position');
    yy_range = get(ha,'YLim');
    zz_range = get(ha,'ZLim');
    
    % find density map of native speech, lateral side
    native_xyz = [];
    all_xyz = [];

    for si = SIDs
        sid = si{1};
        if strcmpi(imgall.(sid).hemi,hemi)
            subject = imgall.(sid);
            
            if isfield(subject.img_mni, 'elecmatrix')
                elecmatrix = subject.img_mni.elecmatrix;
            
                ch_sid = mand_encoding.el(strcmp(mand_encoding.SID, sid));
                if native
                    ch_sel = mand_encoding.el(find(mand_encoding.resp_mn));
                else
                    ch_sel = mand_encoding.el(find(mand_encoding.resp_en));
                end
    
                native_xyz = [native_xyz; elecmatrix(intersect(ch_sel, ch_sid),:)]; 
                all_xyz = [all_xyz; elecmatrix]; 
                numsid=numsid+1;
            end
        end
    
    end
    
    yye = min(native_xyz(:,2))-15:1:max(native_xyz(:,2)+15);
    zze = min(native_xyz(:,3))-15:1:max(native_xyz(:,3)+15); 
    ds = histcounts2(native_xyz(:,2),native_xyz(:,3),yye,zze);
    ds_all = histcounts2(all_xyz(:,2),all_xyz(:,3),yye,zze);
    ds_norm = ds./ds_all;
    ds_norm(isnan(ds_norm)) = 0;

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
end

figure;
if ~native
    colormap(flipud(reds));
else
    colormap(flipud(blues));
end
caxis([0 0.3]);
axis off;
colorbar();

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh Dcons *wrd*;

%% correlate TRF weights between mandarin and english stimuli  

% for every speech responsive electrode, load in STRF and TRF weights 


%% ------------------------- FUNCTIONS ---------------------------------- %%

function [out] = loadOutFile(datapath, corpus, SID, block, CARflag)

    if nargin<5, CARflag=''; end

    filedir= fullfile([datapath SID '/' corpus '/block_z/']);
    if nargin>3
        filename = [SID ...
        '_B' num2str(block) '_HilbAA_70to150_8band_' CARflag 'all_0_mel_' ...
        upper(corpus) '_zflag_global_out_resp_log'];
    else
        sd = dir(filedir);
        filename = sd(3).name;
    end

    load([filedir filename], 'out');
end

function [offset_silences, silence_durs] = getSilenceOffsets(sound, fs, ...
    threshold, min_silence)

    if nargin<4, min_silence=0.5; end
    if nargin<3, threshold=0.1; end

    % find silences (time points to segment at between sentences)
    silences = abs(sound)<threshold;

    % find contiguous segments of silence
    silences = bwlabel(silences);

    % for the end of every contiguous segment of silence indicate length of silence
    offset_silences = zeros(length(silences), 1);
    for j = 1:max(silences)
        % find the last position of the contiguous segment of silence
        offset = find(silences==j, 1, 'last');
        if sum(silences==j)/fs > min_silence
            offset_silences(offset) = sum(silences==j)/fs;
        end
    end

    % remove last silence 
    last_silence = find(offset_silences>0, 1, 'last');
    offset_silences(last_silence) = 0;
    silence_durs = offset_silences(offset_silences>0);
    offset_silences = find(offset_silences>0);
    
    % for debugging
    debug = 1;
    if debug
        figure;
        plot(sound);
        hold on;
        % plot the silence offset_silences
        for j = 1:length(offset_silences)
            xline(offset_silences(j), 'r', 'LineWidth', silence_durs(j)*3);
            xline(offset_silences(j), 'k');
        end
    end
end

function [segOut] = resegmentOut(out, offset_silences, bef, aft)

    if nargin<3, bef = 0.5; end
    if nargin<4, aft = 0.5; end

    % translate offset_silences into indices based on sampling rate
    dataf = out.dataf;
    soundf = out.soundf;
    offset_silences_dataf = round((offset_silences/soundf)*dataf);
    fields = fieldnames(out);

    % segment all the fields in out
    segOut = struct();
    % make empty fields
    
    % hard code for aud, resp, sound
    for j = 1:length(offset_silences)-1
        tmp = struct();

        for field = fields'
            if ~ismember(field{1}, {'aud', 'resp', 'sound', 'name', 'duration'})
                tmp.(field{1}) = out.(field{1});
            end
        end

        % take clip from onset of sound to the following onset of sound
        sound =  out.sound(offset_silences(j):offset_silences(j+1));
        silences = abs(sound)<0.05;
        silences = bwlabel(silences);
        % find the onset of the last contiguous segment of silence 
        % otherwise known as the offset of sound
        offset =  find(silences==max(silences), 1, 'first')/soundf;
   
        % sound
        segment = max(offset_silences(j)-bef*soundf, 1):...
            min(offset_silences(j)+(offset+aft)*soundf, length(out.sound));
        % ensure this is within range
        tmp.sound = out.sound(round(segment));

        % dur
        tmp.duration = length(segment)/soundf;

        % aud
        segment = max(offset_silences_dataf(j)-bef*dataf, 1):...
            min(offset_silences_dataf(j)+(offset+aft)*dataf, size(out.aud, 2));
        tmp.aud = out.aud(:, segment);

        % resp
        tmp.resp = out.resp(:,segment, :);

        % name
        tmp.name = [out.name '_s' num2str(j)];

        if j==1
            segOut = tmp;
        else
            segOut = [segOut; tmp];
        end
    end
end



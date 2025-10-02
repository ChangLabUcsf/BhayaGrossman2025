%% Set up

out_crosscomp_startup;

% % Note - EC202 has no STG coverage
% dSIDs = {'HS11', 'EC237', 'EC266', 'EC197', 'EC261', 'HS8', 'HS9', 'HS10', };
% prof_all = [4, 4, 3, 2, 0, 0, 0, 0];
dSIDs = [profile_all{5}.SID; profile_all{3}.SID; profile_all{4}.SID; ]; % 
prof_all = [profile_all{5}.EnglishProf; profile_all{3}.EnglishProf; profile_all{4}.EnglishProf];
dLang_all = [profile_all{5}.DominantLanguage; profile_all{3}.L1; profile_all{4}.L1];

[~, idx] = sort(prof_all);
dSIDs = dSIDs(idx);
prof_all = prof_all(idx);
dLang_all = dLang_all(idx);

% get color map
% spec = spectral(8);
% spec = flipud(spec([1:3 5:7], :));
% cols = arrayfun(@(x) spec(x+1, :), prof_all, 'UniformOutput',false);
% cols_all = cat(1, cols{:});

% changed colormap
spec = flipud([35, 100, 170; 61, 165, 217; 115, 191, 184; 254, 198, ...
    1; 234, 115, 23; 234, 115, 23;]./256);
cols = arrayfun(@(x) spec(x+1, :), prof_all, 'UniformOutput',false);
cols_all = cat(1, cols{:});

timit_details = load('out_sentence_details_timit_all_loudness.mat');
% asccd_details = load('stim_info/out_sentence_details_acssd_loudness.mat');
% tps = 50:55;

bef=50;
aft=50;

% loading in subject data
TDwrd = loadDwrdMulti('timit', bef, aft, dSIDs, timit_details);

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;
%% A - Create world map with colored countries
mapfig = figure('Name', 'World Map', 'Position', [100, 100, 1000, 600]);
mapax = axes('Parent', mapfig);

% Create world map
ax = worldmap('World');
setm(ax, 'Frame', 'off', 'Grid', 'off');
mlabel off
plabel off

% Load country borders
load coastlines
geoshow(ax, coastlat, coastlon, 'Color', 'black', 'LineWidth', 0.5);

% generate colormap for each language using muted rainbow colors
cmap = containers.Map();
% Convert HSV to RGB for a muted rainbow palette
colortmp = lines(6);
cmap('Spanish') = colortmp(1,:);    % red
cmap('English') = colortmp(2,:);    % orange
cmap('Mandarin') = colortmp(3,:);   % yellow
cmap('Arabic') = colortmp(4,:);     % light green
cmap('Russian') = colortmp(5,:);    % blue
cmap('Korean') = colortmp(6,:);     % purple
cmap('Catalan') = colortmp(1,:);    % pink

% read in cia_languages_csv.txt as a table from the world factbook
allLanguages = readtable('cia_languages_csv.txt');

languages  = {"English", "Spanish", "Mandarin", ...
    "Arabic", "Russian", "Korean", "Catalan"};

languageCountryInfo = table();
% only keep rows in languageCountryInfo where Language is in languages
for i = 1:length(languages)
    if ismember(languages{i}, allLanguages.Primary_Language)
        languageCountryInfo = [languageCountryInfo; ...
            allLanguages(allLanguages.Primary_Language == languages{i}, :)];
    end
end

% make Primary_Language_Percentage a percentage (remove the % sign and cast to double)
for i = 1:length(languageCountryInfo.Primary_Language_Percentage)
    if strcmp(languageCountryInfo.Primary_Language_Percentage(i), 'N/A')
        languageCountryInfo.PrimaryPercent(i) = NaN;
    else
        prct = languageCountryInfo.Primary_Language_Percentage(i);
        languageCountryInfo.PrimaryPercent(i) = ...
            double(string(prct{1}(1:end-1)));
    end
end

% Read world shapefile
% Note: You'll need to have the world shapefile data
% You can download it from Natural Earth or use MATLAB's built-in data

% Using a shapefile of countries
worldData = shaperead('110m_cultural/ne_110m_admin_0_countries.shp', 'UseGeoCoords', true);
% for each country remove trailing spaces
countryNames = {worldData.SOVEREIGNT};
countryNames = cellfun(@(x) string(deblank(x)), countryNames, 'UniformOutput', false);
% Assign country names to the worldData struct array
for i = 1:length(worldData)
    worldData(i).CountryName = countryNames{i};
end

% remove duplicate countries by keeping only the first occurrence of each country name
% remove the second occurrence of 'United States of America'
idx = find(strcmp([worldData.CountryName], 'United States of America'));
worldData(idx(2:end)) = [];

% remove all fields except CountryName, X, Y, Geometry, BoundingBox, SUBREGION
worldData = rmfield(worldData, setdiff(fieldnames(worldData), ...
    {'CountryName', 'Lon', 'Lat', 'Geometry', 'BoundingBox', 'SUBREGION', ...
    'SOV_A3', 'POP_EST', 'CONTINENT'}));

% Plot each country from the shapefile
% one figure for the pie charts, one for the map
% Create or get the map figure handle

% piefig = figure('Name', 'Pie Charts');
% pieax = axes('Parent', piefig);

% First create all subplots
% figure(piefig);
% for i = 1:36  % 6x6 grid
%     subplot(6, 6, i);
%     axis off;  % Hide axes for cleaner look
% end

% Now plot the pie charts
ctr = 1;
for i = 1:length(worldData)  
    countryName = worldData(i).CountryName;

    if ismember(countryName, languageCountryInfo.Country)
        
        % on top of each country plot a pie chart that shows the percentage of
        % people who speak each language in the country
        prct = languageCountryInfo.PrimaryPercent(languageCountryInfo.Country == countryName, :);
        key = languageCountryInfo.Primary_Language(strcmpi(languageCountryInfo.Country, ...
                countryName));

        % make the lightness of the color inversely proportional to the percentage
        % of people who speak the language
        % desiredColor = cmap(key{1});
        
        % if ~isnan(prct)
        %     desiredColor = desiredColor * (100/prct);
        % else
        %     desiredColor = desiredColor * (100/50);
        % end
        % desiredColor = desiredColor + (1-desiredColor) * 0.5;
        % % make sure all values are between 0 and 1
        % desiredColor = max(min(desiredColor, 1), 0);

        desiredColor = [0 0 0];
        if ~isnan(prct)
            desiredColor = desiredColor + 1-(prct/100);
        else
            desiredColor = desiredColor + 0.5;
        end

        % if ~isnan(prct) & prct > 50
        %     % Get the population of the country 
        %     pop = worldData(i).POP_EST;
        %     disp(['Population of ' countryName{1} ' is ' num2str(pop), ' million']);
        %     % height = log(pop)/100;
        %     % width = height;
        %     % x = ax.Position(1);
        %     % y = ax.Position(2);
        %     %ax.Position = [x, y, width, height];
        %     % make the pie chart proportional to the population
        %     figure(piefig);
        %     nexttile;
        %     p = pie([prct/100 1-(prct/100)], [1, 1]);
        %     % remove the labels from the pie chart
        %     p(1).FaceColor = cmap(key{1});
        %     p(3).FaceColor = [1 1 1];
        %     p(3).EdgeColor = 'none';
        %     title(countryName, 'FontSize', 12, 'FontWeight', 'normal');
        %     ctr = ctr + 1;
        % end
    else
        desiredColor = [1 1 1];  % Default gray
    end
    
    % Plot the country on the map
    figure(mapfig);
    if ~isnan(prct) & prct>30 & prct<80
        geoshow(worldData(i).Lat, worldData(i).Lon, ...
            'DisplayType', 'polygon', ...
            'FaceColor', desiredColor, ...
            'EdgeColor', 'black', ...
            'LineWidth', 0.5, 'DisplayName', key{1})
    else
        geoshow(worldData(i).Lat, worldData(i).Lon, ...
            'DisplayType', 'polygon', ...
            'FaceColor', desiredColor, ...
            'EdgeColor', 'black', ...
            'LineWidth', 0.5,'HandleVisibility', 'off')
    end
end

% Add title
title('World Map with Colored Countries', 'FontSize', 16);

% When you need to switch back to the map figure later
figure(mapfig);  % This makes the map figure current
%% Load in unique variances

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

modelnames = {modelnames_timit};

% 'consfeat', 'formant', 
lang = {'eng'};
corpus = {'timit'};

timit_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');
dimex_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
asccd_elecs = load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat');

% determine unique variance per feature and primary encoding
varnames = {'SID', 'el', 'ls', 'prof', ...
    [lang{1} '_base_rsq'], ...
    [lang{1} '_full_rsq'], ...
    [lang{1} '_uv_all'], ...
    [lang{1} '_full_beta'], ...
    };
wordsurp_encoding =  array2table(zeros(0, length(varnames)), 'VariableNames', varnames);

for s = dSIDs'
    SID = s{1}; 
    ls = dLang_all(ismember(dSIDs, SID)); 
    prof = prof_all(ismember(dSIDs, SID));
    corpusStrf{1} = loadMultModelStrf(SID, modelnames{1}, corpus{1}, ...
        datapath, 1, 'v5');

    if ~any(cellfun(@(x) isempty(x), [corpusStrf{1}]))

        % find minimum test R
        minel = min(cellfun(@(x) length(x{3}.meanTestR), corpusStrf));
        els = timit_elecs.allidx.(SID);


        base = cell(2, 1);
        full = cell(2, 1);
        uvall = cell(2, 1);
        fullBeta = cell(2, 1);

        % preallocating cell space
        % full models
        noOns = cell(2, 1);
        noCons = cell(2, 1);
        noPeakr = cell(2, 1);
        noForm = cell(2, 1);
        noWordSurp = cell(2, 1);
        noWord = cell(2, 1);
        noSurp = cell(2, 1);

        % unique variance
        uvOns = cell(2, 1);
        uvCons = cell(2, 1);
        uvPeakr = cell(2, 1);
        uvForm = cell(2, 1);
        uvWordsurp = cell(2, 1);
        uvWord = cell(2, 1);
        uvSurp = cell(2, 1);

        for l = 1:1
            % full models
            base{l} = (corpusStrf{l}{5}.meanTestR.^2)';
            full{l} = (corpusStrf{l}{8}.meanTestR.^2)';
            fullBeta{l} = corpusStrf{l}{8}.strf;

            % single feature excluded models
            noOns{l} = (corpusStrf{l}{1}.meanTestR.^2)';
            noCons{l} = (corpusStrf{l}{2}.meanTestR.^2)';
            noPeakr{l} = (corpusStrf{l}{3}.meanTestR.^2)';
            noForm{l} = (corpusStrf{l}{4}.meanTestR.^2)';
            
            % single feature excluded models (surprisal)
            noWordSurp{l} = (corpusStrf{l}{5}.meanTestR.^2)';        
            noSurp{l} = (corpusStrf{l}{6}.meanTestR.^2)';
            noWord{l} = (corpusStrf{l}{7}.meanTestR.^2)';
            
            % unique variances
            uvOns{l} = full{l}-noOns{l};
            uvPeakr{l} = full{l}-noPeakr{l};
            uvForm{l} = full{l}-noForm{l};
            uvCons{l} = full{l}-noCons{l};
            uvWordsurp{l} = full{l}-noWordSurp{l};

            % this one is calculated without the surprisal feature
            % uvWord{l} = noSurp{l}-noWordSurp{l};

            % now calculated with the surprisal feature
            uvWord{l} = noSurp{l}-noWordSurp{l};
            uvSurp{l} = noWord{l}-noWordSurp{l};
        end

        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        profs = repmat(prof, length(els), 1);

        for l = 1:1
            betaCell{l} = squeeze(mat2cell(fullBeta{l}{1}(:, :, els), ...
                size(fullBeta{l}{1}, 1), size(fullBeta{l}{1}, 2), ones(length(els),1)));
        end


        tmp = table(sids, els, lss, profs, ...
            base{1}(els), full{1}(els), ...
            [uvOns{1}(els), uvPeakr{1}(els), uvForm{1}(els), uvCons{1}(els), ...
                uvWordsurp{1}(els), uvWord{1}(els), uvSurp{1}(els)], betaCell{1}, ...
            'VariableNames', varnames);

        wordsurp_encoding = [wordsurp_encoding; tmp];
    else
        warning(['Missing subject ' SID]);
    end
end

% % load in p-values from permutation testing
% en_wordsurp_pval = nan(1, height(wordsurp_encoding));
% prefix = 'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_';
% permodel = [ prefix 'engSurpNoOnsBin_wordFreqLog'];
% for s = unique(wordsurp_encoding.SID)'
%     SID = s{1};
%     idx = strcmp(wordsurp_encoding.SID, SID);
%     elidx = wordsurp_encoding.el(idx);
% 
%     pvalpath=fullfile(datapath, 'permTest_wordSurp', SID); % c{1}, 
%     cmod=dir(fullfile(pvalpath, '*_zX*_*mat')); 
%     permidx = find(contains({cmod.name}, permodel));
%     permfname=cmod(permidx).name;
% 
%     % load pvalues from permutation testing
%     pvals = load(fullfile(pvalpath, permfname), 'pval');
%     pvals = pvals.pval;
%     en_wordsurp_pval(idx) = pvals(elidx);
% end
% wordsurp_encoding.eng_wordsurp_pval = en_wordsurp_pval';

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *_details *wrd;

%% Run model uv comparison


%% Vis: uv comparison across subject

thresh = 0.001;

featureOrd = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word'};
features = {'word+surp', 'formant', 'consonant', 'peakrate'}; %'peakrate', , 
% titles = {'peakrate', 'phonetic features', 'word+surp'};

figure;
ctr=1;
fields = {'eng_uv_all'}; 
for feat = features
    index = ismember(featureOrd, feat);

    % show individual electrodes per subject
    field = fields{1};    
    y = wordsurp_encoding.(field)(:, index);   
    meethresh = y>thresh;

    uvmed=nan(length(dSIDs), 1);
    subplot(2, length(features), ctr)
    for i = 1:length(dSIDs)
        SID = dSIDs{i};
        
%         boxchart(repmat(i, sum(sidx), 1), y(sidx), ...
%             'MarkerStyle','none', ...
%             'BoxFaceColor', cols_all(i, :)); hold on;
        sidx = strcmp(wordsurp_encoding.SID, SID);
        h = bar(i, sum(sidx & meethresh)/sum(sidx)); 
        hold on; box off;
        h.FaceColor = cols_all(i, :);
        h.EdgeColor = "none";
        set(gca, 'FontSize', 15, 'XTickLabelRotation', 60)

        uvmed(i) = sum(sidx & meethresh)/sum(sidx);

%         % plot those electrodes that do not meet threshold
%         scatter(repmat(i, sum(sidx & ~meethresh), 1) ...
%             - 0.1 + 0.2 *rand(sum(sidx & ~meethresh), 1), ...
%             y(sidx & ~meethresh), 5, [0.7 0.7 0.7], 'filled', ...
%             'MarkerFaceAlpha', 0.6);
%         hold on;
% 
%          % plot those electrodes that do meet threshold
%         scatter(repmat(i, sum(sidx & meethresh), 1) ...
%             - 0.1 + 0.2*rand(sum(sidx & meethresh), 1), ...
%             y(sidx & meethresh), 10, cols_all(i, :), 'filled', ...
%             'MarkerFaceAlpha', 0.6);
% 
%         % find median for uvs greater than zero
%         uvmed(i) = median(y(sidx & meethresh));
%         line([i-0.25 i+0.25], [uvmed(i) uvmed(i)], 'LineWidth', ...
%             3, 'Color', cols_all(i, :));   
%         ylim([-0.01 0.025]);
    end

    % Plot a best fit line through subject medians
%     coefs = polyfit(prof_all, uvmed, 1);
    
    % Use polyval to calcualte function values of that fit
%     plot(1:0.1:length(dSIDs), fliplr(polyval(coefs, 1:0.1:length(dSIDs))),'-', ...
%         'LineWidth',2.5, 'Color',[0.6 0.6 0.6]);
    clear coefs

    title(feat);
    xticks(1:length(dSIDs));
    xticklabels(dLang_all);
    ctr=ctr+1;

    % show proportion of electrodes per subject
    ax = subplot(2, length(features), ctr+length(features)-1);

    perc = nan(1, length(dSIDs));
    count = nan(1, length(dSIDs));
    total = perc;
    for i = 1:length(dSIDs)
        SID = dSIDs{i};
        sidx = strcmp(wordsurp_encoding.SID, SID);

        % percentage of electrodes with positive unique values
        perc(i) = sum(meethresh & sidx) / sum(sidx);
        count(i) = sum(meethresh & sidx);
        total(i) = sum(sidx);
    end
    scatter(1:length(dSIDs), perc, total*5, cols_all, ...
            'filled', 'MarkerFaceAlpha', 0.5); hold on;
    lsline;
    title(feat);
    xticks(1:length(dSIDs));
    xlim([0 9])
    xticklabels(dLang_all);

    tbl = wordsurp_encoding(:, {'SID', 'prof', 'el'});
    tbl.uv = y;
    lme = fitlme(tbl,'uv~1+prof+(1|SID:el)');
    disp(lme)
    disp([feat{1} ': ' num2str(lme.Coefficients.pValue(2))]);
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: uv comparison within subject
thresh = 0.0001;

featureOrd = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word', 'surp'};
features = {'peakrate', 'formant', 'consonant', 'word+surp'}; %'peakrate', , 
% titles = {'peakrate', 'phonetic features', 'word+surp'};

fields = {'eng_uv_all'}; 
cols = lines(length(features));
for i = 1:length(dSIDs)
    SID = dSIDs{i};
    sidx = strcmp(wordsurp_encoding.SID, SID);

    figure;
    j = 1;
    for feat = features
        index = ismember(featureOrd, feat);

        % show individual electrodes per subject
        field = fields{1};    
        y = wordsurp_encoding.(field)(:, index);   
        meethresh = y>thresh;

        boxchart(repmat(j, sum(sidx & meethresh), 1), y(sidx & meethresh), ...
            'MarkerStyle','none', 'BoxFaceColor',cols(j, :)); hold on;

        scatter(repmat(j, sum(sidx & meethresh), 1) - 0.1 + ...
            0.2*rand(sum(sidx&meethresh), 1), ...
            y(sidx & meethresh), 35, cols(j, :), 'filled');

        j = j+1;
    end
    set(gca, 'YScale', 'log');
    title([upper(dLang_all{i}) ', english proficiency: ' ...
        num2str(prof_all(i)) '/5']);
    xticks(1:length(features));
    xticklabels(features);  
    ylabel('feature \Delta R^2');
    xlabel('speech feature');
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;

%% vis: native brains

% initialize design electrode structure
fieldnames = {'English'};
fields = {'eng_uv_all'}; 
thresh = 0.001;

% label = 'word+surprisal';
% label = 'peakrate';
% label = 'word';
feat = 'formant';
% uv feature order
featureOrd = {'onset', 'peakrate', 'formant', 'consonant', 'word+surp', 'word'};
feats = { 'word+surp'}; % 'peakrate', 'formant', 'consonant',

fig = figure();
for lang = 1:1
    for f = 1:length(feats)
        feat = feats{f};
        index = find(ismember(featureOrd, feat));
    
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
    
        % make desel structure
        desel=struct();
        desel.conds = 1:7;
        ls = [1, 2, 4]; % can only do 1-2
        
        % size and color
        desel.sz = [35; 35*ones(length(desel.conds), 1)]; %1:20:200; %ones(1, 10)*0.00000001; %1
        
        % split up peak rate and phonetic features again for MNI plotting wordsurp_encoding.ls(x)
        desel.labels = [];
        desel.yval = arrayfun(@(x) wordsurp_encoding.(fields{lang})(x, index), ...
            1:height(wordsurp_encoding));
        
        % discretize values
        % yvals = sort(desel.yval(desel.yval>0 & ismember(wordsurp_encoding.ls,ls)'));
        % binedges = yvals(1:ceil(length(yvals)/8):length(yvals));
        % [~, binedges] = discretize(desel.yval(desel.yval>0 & ismember(wordsurp_encoding.ls,ls)'), ...
        %     length(desel.conds)-1);
        % manual non-linear edges
        binedges = [0:0.005:0.01 0.015:0.015:0.045];
        binedges = [-1, binedges];
        
        for s=unique(wordsurp_encoding.SID)'
            SID = s{1};
            idx = strcmp(wordsurp_encoding.SID, SID);
            desel.(SID).elid = wordsurp_encoding.el(idx);
            desel.(SID).condition = discretize(desel.yval(idx), ...
                binedges);
        end
        
        desel.cols = [1 1 1; [linspace(0, featcol(1), length(binedges)); ...
            linspace(1, featcol(2), length(binedges)); ...
            linspace(1, featcol(3), length(binedges))]'];
        if strcmp(feat, 'word+surp')
                
            gns = flipud(fpurple(length(binedges)-2));
            desel.cols = [0 0 0 ; gns];
    
            % colorbar
            figure;
            colormap(desel.cols);
            colorbar;            
        end
        
%         lsid = find(ismember(wordsurp_encoding.ls, ls));
        SIDs = dSIDs;
    
        [native_plot] = plotNativeElec(SIDs, desel, 1);
        sgtitle(lang)
    end
end

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding*;


%% run word boundary decoding

Swrd = TDwrd;
corpus = 'timit';

% single window before onset
startp = 31;
timing = 57;
timelabel = '600ms';

tps = startp(1):startp(end)+max(timing);
nreps = 15;

tic
varnames = {'SID', 'elecs', 'ls', 'acc', 'weights', 'auc'}; 

wordBoundary_logisticwrapper(Swrd, dSIDs, wordsurp_encoding, corpus, ...
    startp, timing, timelabel, 'diverse', 1900, nreps);

% nreps 15, mintrials 1900

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* *wrd*;


%% vis: word boundary decoding

corpus = 'TIMIT';
timelabel = '600ms';

filename = [corpus '_word_decode_diverse_' timelabel '_bysubj.mat']; % dimex filename
load([datapath 'Figure3/decode/' filename], 'decode_details');

[~, idx] = sort(prof_all);

ctr = 1;
for i = idx'
    nreps = length(decode_details.tbl.auc(i, :));
%     boxchart(i*ones(nreps,1 ), decode_details.tbl.auc(i, :), ...
%         'BoxFaceColor', cols_all(i, :), 'MarkerColor','k', ...
%         'LineWidth',2, 'JitterOutliers','on', 'MarkerStyle','.'); hold on;
    scatter(ctr*ones(nreps,1)-0.1+rand(nreps, 1)*0.2, decode_details.tbl.auc(i, :), ...
        15,cols_all(i, :), "filled", 'HandleVisibility','off', ...
        'MarkerFaceAlpha', 0.3); hold on;
    aucmean = mean(decode_details.tbl.auc(i, :)); % auc(i, :)
    line([ctr-0.25 ctr+0.25], [aucmean aucmean], 'LineWidth', 3, 'Color', cols_all(i, :));
    ctr=ctr+1;
end

xticks(1:length(dSIDs));
yline(0.5);
xticklabels(decode_details.tbl.ls(idx));
l = legend();

title(l, 'English Proficiency');

set(gca, 'FontSize', 13);
ylabel('AUC');
ylim([0.4 0.7]);
yticks([0.3 0.9]);
xlim([0 18]);
xlabel('Native Language');

% run stats
auc = [decode_details.tbl.auc];
prof = repmat(prof_all, 1, nreps);
elecs = repmat(cellfun(@(x) length(x), decode_details.tbl.elecs), 1, nreps);
% trials = 

lme_tbl = table();
lme_tbl.auc = auc(:);
lme_tbl.prof = prof(:);
lme_tbl.elecs = elecs(:);
lme = fitlme(lme_tbl,'auc~1+prof+(1|elecs)');
disp(lme)

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;
%%
corpus = 'TIMIT';
timelabel = '600ms';

filename = [corpus '_word_decode_diverse_' timelabel '_bysubj.mat']; % dimex filename
load([datapath 'ecog_decode/wordOnset/' filename], 'decode_details');

[~, idx] = sort(prof_all);

ctr = 1;
for i = idx'
    nreps = length(decode_details.tbl.auc(i, :));
%     boxchart(i*ones(nreps,1 ), decode_details.tbl.auc(i, :), ...
%         'BoxFaceColor', cols_all(i, :), 'MarkerColor','k', ...
%         'LineWidth',2, 'JitterOutliers','on', 'MarkerStyle','.'); hold on;
    % scatter(ctr*ones(nreps,1)-0.1+rand(nreps, 1)*0.2, decode_details.tbl.auc(i, :), ...
    %     15,cols_all(i, :), "filled", 'HandleVisibility','off', ...
        % 'MarkerFaceAlpha', 0.3); hold on;
    aucmean = mean(decode_details.tbl.auc(i, :)); % auc(i, :)
    scatter(ctr, aucmean, length(decode_details.tbl.elecs{i})*5, ...
        cols_all(i, :), 'filled', 'MarkerFaceAlpha', 0.8); hold on;
    % line([ctr-0.25 ctr+0.25], [aucmean aucmean], 'LineWidth', 3, 'Color', cols_all(i, :));
    ctr=ctr+1;
end

xticks(1:length(dSIDs));
yline(0.5);
xticklabels(decode_details.tbl.ls(idx));
l = legend();

title(l, 'English Proficiency');

set(gca, 'FontSize', 13);
ylabel('AUC');
ylim([0.45 0.65]);
yticks([0.3 0.9]);
xlim([0 18]);
xlabel('Native Language');

% run stats
auc = [decode_details.tbl.auc(:, :)];
prof = repmat(prof_all, nreps, 1)';

lme_tbl = table();
lme_tbl.auc = auc(:);
lme_tbl.prof = prof(:);
lme = fitlme(lme_tbl,'auc~1+prof');
disp(lme)

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%%


% find top 3 weighted electrodes and look at word erps
f = figure; 
pthresh = 0.05;

% English
SIDs = {'HS11', 'HS11'};% 'EC100', 'EC100', 'EC100','EC100'};
els = [104:105];%, 56, 135, 151, 22, 70, 71, 150]; %236, 205 for EC260, 178, 175

plotSingleTrial = 0;
numel = length(els);
Swrds = TDwrd;

for ctr = 1:length(els)
    SID = SIDs{ctr};    

    dummy = struct();
    nanidx = cellfun(@(x) isempty(x), Swrd.(SID));% ...
        %| (Swrd.syll<2 & ~isnan(Swrd.syll));
    dummy.(SID).resp = cat(3, Swrd.(SID){~nanidx});

    baseline = mean(dummy.(SID).resp(els(ctr), 30:50, :), 'all');
    dummy.(SID).resp  = dummy.(SID).resp - baseline;

    dummy.wordOns = Swrd.wordOns (~nanidx);
    dummy.syllOns = ones(sum(~nanidx), 1);

    ls = wordsurp_encoding.ls(find(strcmp(wordsurp_encoding.SID, SIDs{1}), 1));
    cols = [0 0 0; 0.1, 0.1 0.9]; 

    subplot(1, numel, ctr)
    addpath(genpath('shadederror'))
    plotWordErp(dummy, SID, els(ctr), ...
        [], f, cols, 1, 0.5, 1); hold on;
    ylabel('HFA (z)');
    set(gca, 'FontSize', 13);

    [fvals, betweenVar, withinVar, df1, df2] = Fstat_TIMIT(...
        dummy.(SID).resp(els(ctr), :, :), dummy.wordOns+1, [1, 2]);
    corrected_pval = pthresh / size(dummy.(SID).resp, 2);
    fthresh = finv(1-corrected_pval, df1, df2);  

    x = -0.5:0.01:0.5;
    scatter(x(fvals>fthresh), 0.1*ones(1, sum(fvals>fthresh)), 45, ...
        fvals(fvals>fthresh), 'filled', 'HandleVisibility', 'off');
    cm = colormap("gray");
    colormap(flipud(cm(1:200, :)))
    % ylim([0 1.3]);
    ylim([-0.3 0.8]);
    xlim([-0.2, 0.21]);
    h=xline(0);
    h.Color = 'k';
    legend('off');
end
clear dummy

%% ------------------------------ functions -------------------------------


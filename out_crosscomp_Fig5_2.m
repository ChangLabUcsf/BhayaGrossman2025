% Ilina Bhaya-Grossman
% 11.19.2025
% WARNING: Running this script will delete all variables in your current
% workspace. Proceed with caution.
out_crosscomp_startup;

dSIDs = [profile_all{5}.SID; profile_all{3}.SID; profile_all{4}.SID;]; 
prof_all = [profile_all{5}.EnglishProf; profile_all{3}.EnglishProf; profile_all{4}.EnglishProf];
dLang_all = [profile_all{5}.DominantLanguage; profile_all{3}.L1; profile_all{4}.L1];

[~, idx] = sort(prof_all);
dSIDs = dSIDs(idx);
prof_all = prof_all(idx);
dLang_all = dLang_all(idx);

% changed colormap
spec = flipud([35, 100, 170; 61, 165, 217; 115, 191, 184; 254, 198, ...
    1; 234, 115, 23; 234, 115, 23;]./256);
cols = arrayfun(@(x) spec(x+1, :), prof_all, 'UniformOutput',false);
cols_all = cat(1, cols{:});

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% H - Create world map with colored countries

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
worldData = shaperead('diverse_lang/110m_cultural/ne_110m_admin_0_countries.shp', 'UseGeoCoords', true);
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

% Now plot the pie charts
for i = 1:length(worldData)  
    countryName = worldData(i).CountryName;

    if ismember(countryName, languageCountryInfo.Country)
        
        % on top of each country plot a pie chart that shows the percentage of
        % people who speak each language in the country
        prct = languageCountryInfo.PrimaryPercent(languageCountryInfo.Country == countryName, :);
        key = languageCountryInfo.Primary_Language(strcmpi(languageCountryInfo.Country, ...
                countryName));

        desiredColor = [0 0 0];
        if ~isnan(prct)
            desiredColor = desiredColor + 1-(prct/100);
        else
            desiredColor = desiredColor + 0.5;
        end

    else
        desiredColor = [1 1 1];  % Default gray
    end
    
    % Plot the country on the map
    figure(mapfig);
    if ~all(isnan(prct)) && prct>30 && prct<80
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

clearvars -except *all subj *vow* *details *SIDs datapath bef aft tps ...
    betaInfo* *encoding* allidx fthresh *cons* *wrd;

%% I - Diverse word boundary decoding

corpus = 'TIMIT';
timelabel = '600ms';

filename = [corpus '_word_decode_diverse_' timelabel '_bysubj.mat']; % dimex filename
load([datapath 'Figure3/decode/' filename], 'decode_details');

[~, idx] = sort(prof_all);

ctr = 1;
legendidx = [];
for i = idx'
    nreps = length(decode_details.tbl.auc(i, :));
    scatter(ctr*ones(nreps,1)-0.1+rand(nreps, 1)*0.2, decode_details.tbl.auc(i, :), ...
        15,cols_all(i, :), "filled", 'HandleVisibility','off', ...
        'MarkerFaceAlpha', 0.3); hold on;
    aucmean = mean(decode_details.tbl.auc(i, :)); % auc(i, :)
    if ~ismember(prof_all(i), legendidx)
        vis = 'on';
        legendidx = [legendidx, prof_all(i)];
    else
        vis = 'off';
    end
    line([ctr-0.25 ctr+0.25], [aucmean aucmean], 'LineWidth', 3, 'Color', ...
        cols_all(i, :), 'HandleVisibility', vis);
    ctr=ctr+1;
end

xticks(1:length(dSIDs));
yline(0.5, 'HandleVisibility', 'off');
xticklabels(decode_details.tbl.ls(idx));
l = legend(num2str(unique(prof_all)));
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

%% UNUSED - Diverse word boundary decoding by number of speech responsive
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


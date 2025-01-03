
% load in all subjects
[sSIDs, eSIDs, bSIDs, mSIDs, dSIDs] = getSIDinfo();

% if select_elec folder doesn't exist, create it
if ~exist('select_elec', 'dir')
    mkdir('select_elec');
end

%% for timit electrode selection

% for each set of subjects, run getSpeechElecs
[allidx_es, fvals_es] = getSpeechElecs([eSIDs, sSIDs, bSIDs], 'bychan', 'timit', [], 'dimex');

% add in mSIDs and dSIDs
[allidx_mo, fvals_mo] = getSpeechElecs(mSIDs, 'bychan', 'timit', [], 'asccd');

% add in dSIDs
[allidx_d, fvals_d] = getSpeechElecs(dSIDs, 'bychan', 'timit', [], 'timit');

% copy all fields from each set of subjects into one struct
allidx = struct();
fvals = struct();
subject_sets = {allidx_es, allidx_mo, allidx_d};
fval_sets = {fvals_es, fvals_mo, fvals_d};

for i = 1:length(subject_sets)
    for f = fieldnames(subject_sets{i})'
        allidx.(f{1}) = subject_sets{i}.(f{1});
        fvals.(f{1}) = fval_sets{i}.(f{1});
    end
end

% save output as 
outfile = "select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat";
save(outfile, 'allidx', 'fvals');

clear fvals_es fvals_mo allidx_es allidx_mo allidx fvals

%% for dimex electrode selection

% for each set of subjects, run getSpeechElecs
[allidx, fvals] = getSpeechElecs([eSIDs, sSIDs, bSIDs], 'bychan', 'dimex', [], 'timit');

% save output as
outfile = "select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat";
save(outfile, 'allidx', 'fvals');

clear fvals allidx

%% for Mandarin electrode selection

% for each set of subjects, run getSpeechElecs
[allidx, fvals] = getSpeechElecs([mSIDs, {'EC186', 'EC204', 'EC222'}], 'bychan', 'asccd', [], 'timit');

% save output as
outfile = "select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat";
save(outfile, 'allidx', 'fvals');

clear fvals allidx


%% For Spanish / English 

% get the speech electrodes for each language
load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');
allidx_timit = allidx;
fvals_timit = fvals;

load('select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat');
allidx_dimex = allidx;
fvals_dimex = fvals;

% make a venn diagram of the overlap between the two languages for each subject
SIDs = [sSIDs, eSIDs, bSIDs];
figure;
for i = 1:length(SIDs)
    SID = SIDs{i};
    subplot(5, 6, i);
    if isfield(allidx_timit, SID) && isfield(allidx_dimex, SID)
        timit_elecs = allidx_timit.(SID);
        dimex_elecs = allidx_dimex.(SID);
        overlap = intersect(timit_elecs, dimex_elecs);
        venn([length(timit_elecs), length(dimex_elecs)], length(overlap));
        percent = length(overlap) / length(union(timit_elecs, dimex_elecs));
        title([SID ': ' num2str(percent*100, 2) '%']);
        axis off;
    end
end

%% For Mandarin / English

% get the speech electrodes for each language
load('select_elec/out_elecs_speechtypeftest_bychan_asccd_all.mat');
allidx_asccd = allidx;
fvals_asccd = fvals;

load('select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat');
allidx_timit = allidx;
fvals_timit = fvals;

% make a venn diagram of the overlap between the two languages for each subject
SIDs = [mSIDs, {'EC186', 'EC204', 'EC222'}];
figure;
for i = 1:length(SIDs)
    SID = SIDs{i};
    subplot(2, 4, i);
    if isfield(allidx_asccd, SID) && isfield(allidx_timit, SID)
        asccd_elecs = allidx_asccd.(SID);
        timit_elecs = allidx_timit.(SID);
        overlap = intersect(asccd_elecs, timit_elecs);
        venn([length(asccd_elecs), length(timit_elecs)], length(overlap));
        percent = length(overlap) / length(union(asccd_elecs, timit_elecs));
        title([SID ': ' num2str(percent*100, 2) '%']);
        axis off;
    end
end


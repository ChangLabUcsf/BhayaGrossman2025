addpath(genpath('../../../../ecog_scripts'))
addpath(genpath('../../../../plotting_scripts'))
addpath(genpath('../util'))
zFolder = 'block_z'; % 'block_z'
[datapath, dpath] = setDatapath;
addpath(genpath(datapath))
bef=20;
aft=50;

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs, mSIDs, dSIDs] = getSIDinfo();

% load in all corpus details
timit_details = load([datapath 'stim_info/out_sentence_details_timit_all_loudness.mat']);
dimex_details = load([datapath 'stim_info/out_sentence_details_dimex_all_loudness.mat']);
asccd_details = load([datapath 'stim_info/out_sentence_details_asccd_all_loudness.mat']);
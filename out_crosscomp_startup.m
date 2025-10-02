% addpath(genpath('../../../../ecog_scripts'))
% addpath(genpath('../../../../plotting_scripts'))
addpath(genpath('utils'))
% zFolder = 'block_z'; % 'block_z'

bef=20;
aft=50;
datapath = 'data/';

% Note - EC202 has no STG coverage
[sSIDs, eSIDs, bSIDs, mSIDs, dSIDs] = getSIDinfo();

% load in all corpus details
if ~exist('timit_details', 'var')
    load([datapath 'stim_info/out_sentence_details_timit_all_loudness.mat'], 'timit_details');
end

if ~exist('dimex_details', 'var')
    load([datapath 'stim_info/out_sentence_details_dimex_all_loudness.mat'], 'dimex_details');
end

if ~exist('asccd_details', 'var')
    asccd_details = load([datapath 'stim_info/out_sentence_details_asccd_all_loudness.mat']);
end

%% load in patient language profile data (xlsx)

english_profile = readtable('LanguageProfiles.xlsx', 'UseExcel',true, ...
    'Sheet','EnglishMono');
% remove all rows where L1 is not English and monolingual
english_profile = english_profile(strcmp(english_profile.L1, 'English'), :);
english_profile = english_profile(strcmp(english_profile.Group, 'monolingual'), :);

% load in spanish profile
spanish_profile = readtable('LanguageProfiles.xlsx', 'UseExcel',true, ...
    'Sheet','SpanishMono');
% remove all rows where L1 is not Spanish and monolingual
spanish_profile = spanish_profile(strcmp(spanish_profile.L1, 'Spanish'), :);
spanish_profile = spanish_profile(strcmp(spanish_profile.Group, 'monolingual'), :);

% load in mandarin profile
mandarin_profile = readtable('LanguageProfiles.xlsx', 'UseExcel',true, ...
    'Sheet','MandarinMono');
% remove all rows where L1 is not Mandarin and monolingual
mandarin_profile = mandarin_profile(strcmp(mandarin_profile.L1, 'Mandarin'), :);
mandarin_profile = mandarin_profile(strcmp(mandarin_profile.Group, 'monolingual'), :);

% load in Spanish-English bilingual profile
bilingual_profile = readtable('LanguageProfiles.xlsx', 'UseExcel',true, ...
    'Sheet','SEBilingual_Truncated');
% remove all rows where L1 is not Spanish or English and bilingual
bilingual_profile = bilingual_profile(strcmp(bilingual_profile.Group, 'bilingual'), :);

% load in bilingual-other profile
other_profile = readtable('LanguageProfiles.xlsx', 'UseExcel',true, ...
    'Sheet','Other_Fixed');
% remove all rows where L1 is not Spanish or English and bilingual
other_profile = other_profile(cellfun(@(x) ~isempty(x), other_profile.DominantLanguage), :);

profile_all = {spanish_profile, english_profile, mandarin_profile, ...
    bilingual_profile, other_profile};

%% load in model names

% could add _wordFreqLog to the full model
% load in surprisal values
% using phonetic feature instead of splitting by vowel
modelnames_timit={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % 2 remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % 3 remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_engSurpNoOnsBin', ... % 4 remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... % 5 remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % 6 remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... % 7 remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin', ... % 8 full
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns', ... % 9 remove surprise and length
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordFreqLog', ... 
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_engSurpNoOnsBin_wordFreqLog', ... % add in frequency
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engBiSurpBin', ... % biphone surprisal (review)
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engTriSurpBin', ... % biphone surprisal (review)
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_F0_relPitchBin_F0ChangeBin', ... % pitch (review)
    'onset_phnfeatConsOnset_L_formantMedOnset', ... % continuous envelope (review) 
    }; 
%'onset_phnfeatConsOnset_maxDtL_formantMedOnset_engSurpNoOnsBin', ... %remove word

modelnames_dimex={'phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ...%remove onset
    'onset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove consonant features        
    'onset_phnfeatConsOnset_formantMedOnset_wordOns_wordL_spSurpNoOnsBin', ... %remove peakrate
    'onset_phnfeatConsOnset_maxDtL_wordOns_wordL_spSurpNoOnsBin', ... %remove formant
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset', ... %remove word feat/base
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL', ... % remove surprise
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_spSurpNoOnsBin', ... % remove word only
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin' ... % full
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordFreqLog', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_wordOns_wordL_spSurpNoOnsBin_wordFreqLog', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_spBiSurpBin', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_spTriSurpBin', ...
    'onset_phnfeatConsOnset_maxDtL_formantMedOnset_F0_relPitchBin_F0ChangeBin', ... % pitch (review)
    'onset_phnfeatConsOnset_L_formantMedOnset', ... % continuous envelope (review) }; % full v2
    }; 




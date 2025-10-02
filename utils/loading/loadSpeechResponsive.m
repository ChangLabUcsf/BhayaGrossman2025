function [elec_tbl] = loadSpeechResponsive(SIDs)

    varnames = {'SID', 'el', 'ls'}; % , 'uv_peakRate'
    elec_tbl = array2table(zeros(0, length(varnames)), 'VariableNames', varnames);
    [sSIDs, eSIDs, bSIDs, mSIDs] = getSIDinfo();

    timit_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_timit_all.mat");
    dimex_elecs = load("select_elec/out_elecs_speechtypeftest_bychan_dimex_all.mat");

    % determine unique variance per feature and primary encoding
    for s = SIDs
        SID = s{1}; 
        ls = find(cellfun(@(x) ismember(SID, x), {sSIDs, eSIDs, mSIDs, bSIDs}));

        % Include all speech responsive electrodes instead
        if isfield(dimex_elecs.allidx, SID) && isfield(timit_elecs.allidx, SID)
            els = union(timit_elecs.allidx.(SID), dimex_elecs.allidx.(SID));
        elseif isfield(dimex_elecs.allidx, SID)
            els = dimex_elecs.allidx.(SID);
        elseif isfield(timit_elecs.allidx, SID)
            els = timit_elecs.allidx.(SID);
        end

        sids = repmat({SID}, length(els), 1);
        lss = repmat(ls, length(els), 1);
        
        tmp = table(sids, els, lss,'VariableNames', varnames);
        elec_tbl = [elec_tbl; tmp];
    end
end
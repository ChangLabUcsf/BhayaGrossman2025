function [sSIDs, eSIDs, bSIDs, mSIDs, dSIDs] = getSIDinfo()

    % spanish speaking subjects
    % EC100, EC105, EC116, EC128, EC152, EC163, EC172, EC203, EC214, EC252, EC225
    sSIDs={'EC100', 'EC105', 'EC163', 'EC172','EC214', 'EC252', 'EC225', 'EC116', 'EC320'}; 
    % bad subj: EC128, EC203

    % english speaking subjects
    % EC183, EC186, EC195, EC196, EC204, EC208, EC212, EC219, EC221, EC222, EC230, EC235, EC242
    eSIDs={'EC183', 'EC186', 'EC195', 'EC196', 'EC212',  'EC222', 'EC235', 'EC204'}; 
    % bad subj: EC208, EC242, EC221, EC230
    % 'EC219' --> not clear about timit!

    % spanish-english bilingual speaking subjects
    % EC122, EC129, EC139, EC159, EC161, EC200, EC202, EC260, EC266, EC273,
    % EC274, EC282
    bSIDs={'EC122', 'EC129', 'EC139', 'EC159', 'EC161', 'EC200', 'EC266', ...
        'EC260', 'EC282', 'EC296', 'EC322'}; % EC273, EC274, 
    % updated to exclude porto, noisy & catalan trilingual (161, 139, 266)
    bSIDs={'EC122', 'EC129',  'EC159',  'EC200',  ...
        'EC260', 'EC282', 'EC296', 'EC322'};
    % bad subj: EC202
    
    % mandarin speaking subjects
    mSIDs={'HS8', 'HS9', 'HS10', 'HS11'}; 

    % other languages (Russian, Arabic, Korean)
    dSIDs = {'EC237', 'EC266', 'EC197', 'EC261'};
end
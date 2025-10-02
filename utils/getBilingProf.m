function [proficiency] = getBilingProf()
    varnames = {'SID', 'L1', 'span_prof', 'eng_prof', 'mand_prof', 'other_prof', ...
    'span_freq', 'eng_freq', 'span_age', 'eng_age'};
    proficiency =  array2table(zeros(0, length(varnames)), 'VariableNames', varnames);
    
    % double checked
    % missing language info for: EC100, EC105, EC152
    bSIDs = {'EC122', 'EC129', 'EC139', 'EC159', 'EC161', 'EC202', 'EC200', ...
        'EC260', 'EC266', 'EC282', 'EC292'}';
    
    % Proficiency in each language
    sprof = {5, 5, 5, 3, 3, 2, 5, 5, 5, 3, 5}';
    eprof = {5, 5, -1, 5, 5, 5, 4, 4, 3, 5, 5}'; 
    mprof = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0}';
    other = {3, 0, -1, 0, 4, 0, 0, 0, 5, 0, 0}'; 
    
    % L1, S - Spanish, E - English, O - Other, B - both
    L1 = {'S', 'S', 'S', 'E', 'O', 'S', 'S', 'S', 'S', 'B', 'S'}';
    
    % frequency of use (1- daily, 2-weekly, 3-monthly, 
    % 4-irregular, -1-don't know)
    efreq = {-1, -1, -1, 1, 1, 1, 4, 1, 1, 1, 1}'; 
    sfreq = {-1, -1, -1, 1, 1, 2, 1, 1, 1, 1, 1}'; 
    
    % Age of acquisition
    eaa = {-1, -1, -1, 1, 5, 5, 5, 3, 8, 3, 5}'; 
    saa = {-1, 1, -1, 1, 3, 1, 1, 1, 1, 3, 1}'; 
    
    % Compbining all proficiency, frequency, and AoA information
    proficiency = [proficiency; 
        [bSIDs, L1, sprof, eprof, mprof, other, ...
        sfreq, efreq, saa, eaa]];
end
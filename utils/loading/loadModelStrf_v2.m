function [corpusStrf]=loadModelStrf_v2(modelname, corpus, SID, datapath, scaleflag, ...
    zscoreflag, version, pitchSel)

    if nargin<5, scaleflag=0; end
    if nargin<6, zscoreflag=0; end
    if nargin<7, version='v5'; end
    
    % pitch selection criteria for sentences included in training + test
    % set
    if nargin<8, pitchSel=''; end % 'lte170', 'gt170', ''
    
    if ~strcmp(pitchSel, '')
        pitchSel = ['_' pitchSel];
    end
    
    filename = [SID '_strf_' modelname '_' ...
        'zX' num2str(zscoreflag) '_zY0_scX1_scY0_hp0_SentOns1_sentScale' ...
        num2str(scaleflag) '_El1to*_600ms_edge1_boot1' pitchSel '.mat'];
    fulldir = dir(fullfile(datapath, 'trf', SID, [corpus '/'], filename));
    
    corpusStrf = [];
    try
        corpusStrf = load(fullfile(fulldir.folder, fulldir.name));
    catch
        warning(['No ' corpus ' ' modelname ' strf data exists for ' SID])
    end
end
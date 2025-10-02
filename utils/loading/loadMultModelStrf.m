%% get corpusStrfs
function corpusStrf=loadMultModelStrf(SID, modelnames, corpus, datapath, zscoreflag, version, pitchSel)

    if nargin<5, zscoreflag=0; end
    if nargin<6, version='v5'; end
    if nargin<7, pitchSel=''; end
    
    corpusStrf = cell(1, length(modelnames));
    for i=1:length(modelnames)
        if length(zscoreflag)>1 % can use different zscore flags for the multiple models
            z = zscoreflag(i);
        else
            z = zscoreflag;
        end
        corpusStrf{i} = loadModelStrf_v2(modelnames{i}, corpus, SID, ...
                datapath, 0, z, version, pitchSel);
    end 
end


%% plot erp
function plotResp(resp, desmat, singleSet, SID, f, cols, gridLayout, bef, plotmean)
    if nargin < 5, f = []; end 
    if nargin < 6, cols = ametrine(length(desmat.names)); end    
    if nargin < 9, plotmean = 1; end
    if nargin < 7 || isempty(gridLayout)
        if strcmp(SID, 'EC161') 
            elect = 1:128;
            gridLayout = rot90(reshape(1:length(elect),[8,16])',2);
        elseif strcmp(SID, 'EC225') || strcmp(SID, 'EC249')
            elect = 1:128;
            gridLayout = flipud(reshape(1:length(elect),[16,8])')';  
        elseif strcmp(SID, 'EC252')
            elect = 1:128;
            gridLayout = fliplr(rot90(reshape(elect,[16,8]), 2));
        elseif strcmp(SID, 'EC200')
            elect = 1:64;
            gridLayout = rot90(reshape(1:length(elect),[8,8])',2);
        elseif strcmp(SID, 'EC212') || strcmp(SID, 'EC219')
            elect = 1:432;
            gridLayout = rot90(reshape(1:length(elect),[16, 27])',2);
        elseif strcmp(SID, 'EC237') || strcmp(SID, 'EC260') || strcmp(SID, 'EC186')
            elect = 1:256;
            gr1 = fliplr(rot90(reshape(1:length(elect)/2,[16,8])',2));
            gr2 = fliplr(rot90(reshape(length(elect)/2+1:length(elect),[16,8])',2));
            gridLayout = [gr1; gr2]';
        elseif strcmp(SID, 'EC159') || strcmp(SID, 'EC129')  
            elect = 1:256;
            gridLayout = rot90(reshape(1:length(elect),[sqrt(length(1:length(elect))),...
                sqrt(length(1:length(elect)))])');
        else
            elect = 1:256;
            gridLayout = rot90(reshape(1:length(elect),[sqrt(length(1:length(elect))),...
                sqrt(length(1:length(elect)))])',2);
        end
    end
    if nargin < 8 || isempty(bef), bef = 0.2; end
    
    if isempty(singleSet)
        ecog_erpPlot(resp, desmat, elect, gridLayout,[], ...
            f, bef, [], plotmean, 100, [], cols);   
    elseif length(singleSet) > 20
        ecog_erpPlot(resp, desmat, singleSet, gridLayout,[], ...
            f, bef, [], plotmean, 100, [], cols);       
    else
        [rows, colus] = gridSize(numel(singleSet));
        gridLayout = reshape(singleSet, colus, rows);
        
        if isempty(f)                       
            ecog_erpPlot(resp, desmat,singleSet, ...
                gridLayout,[], [], bef, [], plotmean, 100, [], cols);
        elseif isscalar(singleSet)
            ecog_erpPlotSingle(resp, desmat,[singleSet(1)], ...
                    gridLayout,[], f, bef, [], plotmean, 100, [], cols);
        else
            for i = 1:numel(singleSet)
                % one half of the grid            
                subplot(rows, colus, i)         
                ecog_erpPlotSingle(resp, desmat,[singleSet(i)], ...
                    gridLayout,[], f, bef, [], plotmean, 100, [], cols);
                xticks([]); yticks([]);
            end
        end        
    end
    
end
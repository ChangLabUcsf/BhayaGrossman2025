%% plot electrodes on MNI brain
% desel is like desmat where electrodes are specified under condition
function [mni_plot] = plotMNIElec(SIDs, desel, hemi, plotsil, plotbrain, imgall)
    if nargin<5, plotbrain = 1; end

        
    SIDs = SIDs(cellfun(@(x) strcmp(imgall.(x).hemi, hemi), SIDs));  
    if isempty(SIDs)
        warning('No SIDs included with correct hemisphere type.');
        mni_plot = [];
        return
    end 

    if plotbrain
        fig = figure();
        axh = axes('Parent', fig);
        hold(axh, 'on');

        % plot glass brain    
        ctmr_gauss_plot(imgall.(SIDs{1}).img_mni.cortex,[0 0 0], 0, hemi)
        alpha 0.6;
    end

    % last four columns correspond to coordinates and whether electrode is
    % selected (highlighted)
    varnames = {'SID', 'el', 'cond', 'hemi', 'x', 'y', 'z', 'sel', 'anatomy'};
    mni_plot = array2table(zeros(0,9), 'VariableNames', varnames);
    for s=1:length(SIDs)
        SID=SIDs{s};
        
        % check if montage information exists
        if isfield(imgall, SID)
            img_mni=imgall.(SID).img_mni;
            
            x_add=70;
            % z dimension for specific subjects
            if ismember(SID, {'EC172'}), x_add=140; end
            if ismember(SID, {'EC100'}), x_add=90; end
            if ismember(SID, {'EC195'}), x_add=40; end
            if ismember(SID, {'EC159'}), x_add=150; end
            if strcmp(hemi,'lh'), x_add=x_add-90; end
            
            try 
                % plot electrodes corresponding to each condition
                for i=desel.conds
                    
                    % find electrodes corresponding to current condition
                    condel = find(desel.(SID).condition == i);  
                    if ~isempty(condel)
                        
                        % electrodes in condition, make sure not indexing out of mni elecs
                        elidx = desel.(SID).elid(condel);                        
                        elidx = elidx(elidx<size(img_mni.elecmatrix, 1));
                        n = length(elidx);
                        
                        % plotting in three dimensional space
                        if strcmp(hemi,'lh') 
                            cond_add = -i*10;
                        else
                            cond_add = i*10;
                        end
                        x = img_mni.elecmatrix(elidx, 1)+x_add-10+cond_add;
                        y = img_mni.elecmatrix(elidx, 2);
                        z = img_mni.elecmatrix(elidx, 3);
                        selid = zeros(n, 1);
                        if isfield(desel.(SID), 'selid')
                            selid = ismember(elidx, desel.(SID).selid);
                        end
                        %x = cat(1, x, [b, c]);                        
                        
                        % make sure dimensions match
                        elidx = reshape(elidx, [n, 1]);
                        selid = reshape(selid, [n, 1]);
                        anatomy = img_mni.anatomy(elidx, 4);
                        t2 = table(repmat({SID}, n, 1), elidx, repmat(i, n, 1), ...
                            repmat({hemi}, n, 1), x, y, z, selid, anatomy, 'VariableNames', varnames);
                        mni_plot = [mni_plot; t2];   
                        clear x y z n
                    end
                end 
            catch
                warning(['Missing electrode information for ' SID])
            end
        end
    end
    
    % convert the conditions to be indices for the size field
    % unique sorts
    %condidx = arrayfun(@(x) find(unique(mni_plot.cond)==x), mni_plot.cond);
    condidx = arrayfun(@(x) find(unique(desel.conds)==x), mni_plot.cond);
    
    if plotbrain
        scatter3(mni_plot.x, mni_plot.y, mni_plot.z, desel.sz(condidx), ...
            desel.cols(condidx, :), 'o', 'filled', ...
            'MarkerFaceAlpha', 0.7); % 'MarkerEdgeColor', [0 0 0]
    end

    % can look SID by SID if electrode coloring looks incorrect
    debug = 0;
    if debug
        for s = unique(mni_plot.SID)'
            SID = s{1};
            idx = strcmp(SID, mni_plot.SID);
            disp(sum(idx));
            scatter3(mni_plot.x(idx), mni_plot.y(idx), mni_plot.z(idx), ...
                desel.sz(condidx(idx)), ...
                 desel.cols(condidx(idx), :), 'o', 'filled', ...
                'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0 0 0]); hold on; % 'MarkerEdgeColor', [0 0 0]
        end
    end
    
    if plotbrain
        hold on;
        grid(axh, 'on');    
        
        sel = mni_plot.sel>0;
        scatter3(mni_plot.x(sel), mni_plot.y(sel), mni_plot.z(sel),...
            75, 'r', 'o',  'LineWidth', 1.7); % 'MarkerEdgeColor', 'k', 
    end

    if ~isempty(mni_plot)
        disp('------------------- Stats ---------------')
        disp(['Number of subjects included: ' strjoin(unique(mni_plot.SID))]);
        disp(['Total number of electrodes: ' num2str(height(mni_plot))]);
        disp(['Number of electrodes over cond 1: ' num2str(sum(mni_plot.cond>1))]);   
    end 
end
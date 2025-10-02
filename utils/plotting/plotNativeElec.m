%% plot electrodes on NATIVE brains
% desel is like desmat where electrodes are specified under condition
function [native_plot] = plotNativeElec(SIDs, desel, plot, imgall)
    
    % last four columns correspond to coordinates and whether electrode is
    % selected (highlighted)
    varnames = {'SID', 'el', 'cond', 'hemi', 'x', 'y', 'z', 'sel', 'anatomy'};
    native_plot = array2table(zeros(0,9), 'VariableNames', varnames);
    for s=1:length(SIDs)
        SID=SIDs{s};  
        
        % check if montage information exists and is on correct hemisphere
        if isfield(imgall, SID)
            img_native=imgall.(SID).img_native;
            hemi=imgall.(SID).hemi;
            
            x_add=15;
            % z dimension for specific subjects
%             if strcmp(SID, 'EC172'), x_add=90; end
            if strcmp(hemi,'lh'), x_add=-10; end        
            if strcmp(SID, 'EC221'), x_add=-90; end
            % try 
                % plot electrodes corresponding to each condition
                for i=desel.conds
                    
                    % find electrodes corresponding to current condition
                    condel = find(desel.(SID).condition == i);  
                    if ~isempty(condel) && isfield(img_native, 'elecmatrix')
                        
                        elidx = desel.(SID).elid(condel);
                        % make sure not indexing out of mni elecs
                        elidx = elidx(elidx<size(img_native.elecmatrix, 1));
                        n = length(elidx);
                        
                        % plotting in three dimensional space
                        x = img_native.elecmatrix(elidx, 1)+x_add-10;
                        y = img_native.elecmatrix(elidx, 2);
                        z = img_native.elecmatrix(elidx, 3);
                        selid = zeros(n, 1);
                        if isfield(desel.(SID), 'selid')
                            selid = ismember(elidx, desel.(SID).selid);
                        end
                        % x = cat(1, x, [b, c]);                        
                        
                        % make sure dimensions match
                        elidx = reshape(elidx, [n, 1]);
                        selid = reshape(selid, [n, 1]);
                        if  isfield(img_native, 'anatomy')
                            anatomy = img_native.anatomy(elidx, 4);
                        else
                            anatomy = nan(n, 1);
                        end

                        t2 = table(repmat({SID}, n, 1), elidx, repmat(i, n, 1), ...
                            repmat({imgall.(SID).hemi}, n, 1), x, y, z, selid, ...
                            anatomy, 'VariableNames', varnames);
                        native_plot = [native_plot; t2];   
                        clear x y z n
                    end
                end 
                
                idx = strcmp(SID, native_plot.SID);
                % if there are >0 electrodes to plot for this subject
                if ~isempty(idx) && plot
                    % plot subject specific glass brain
                    fig = figure();
                    axh = axes('Parent', fig);
                    hold(axh, 'on');   
                    ctmr_gauss_plot(img_native.cortex,[0 0 0], 0, imgall.(SID).hemi)

                    % plot all other electrodes on native brain
                    elidx = find(~ismember(1:size(img_native.elecmatrix, 1), ...
                        desel.(SID).elid));
                    scatter3(img_native.elecmatrix(elidx, 1)+x_add, ...
                        img_native.elecmatrix(elidx, 2), ...
                        img_native.elecmatrix(elidx, 3), 5, [0 0 0], 'o', ...
                        'filled', 'MarkerFaceAlpha', 0.8);  hold on;


                    condidx = arrayfun(@(x) find(desel.conds==x), native_plot.cond(idx));               
                    scatter3(native_plot.x(idx), native_plot.y(idx), native_plot.z(idx), ...
                        desel.sz(condidx), desel.cols(condidx, :), 'o', 'filled', ...
                         'MarkerEdgeColor', 'none', ...
                        'LineWidth',1.5); % , 'MarkerEdgeColor', [0 0 0] 'MarkerFaceAlpha', 1,
                    legend(desel.labels);
                    hold on;

                    grid(axh, 'on');  
                    legend('off');

                    sel = native_plot.sel & idx;
                    scatter3(native_plot.x(sel), native_plot.y(sel), native_plot.z(sel),...
                        75, 'k', 'o', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5); 
                    title(SID)
                end
            % catch
            %     warning(['Missing electrode information for ' SID])
            % end
        end
        
    end

    % only if native plot is not empty        
    if ~isempty(native_plot)
        disp('------------------- Stats ---------------')
        disp(['Number of subjects included: ' strjoin(unique(native_plot.SID))]);
        disp(['Total number of electrodes: ' num2str(height(native_plot))]);
        disp(['Number of electrodes over cond 1: ' num2str(sum(native_plot.cond>1))]);  
    end
end
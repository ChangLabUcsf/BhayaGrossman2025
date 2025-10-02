function [chl] = sentenceERP(meanErp, semErp, color, bef, winds, stitchgap)

    lineWidth = 2;
    
    addpath(genpath('shadederror/'))
    dataf = 100;

    % calculating erp and sem
    % meanErp=squeeze(mean(resp, 3, 'omitnan'));
    % semErp=squeeze(nansem(resp, 3));

    if nargin<5, stitchgap = 20; end

    if nargin<4
        time_axis= (0:(size(resp, 2)-1))/dataf - bef;   
    
        chl = shadedErrorBar(time_axis,meanErp,semErp,...
            'lineprops', {'color', color, 'LineWidth', lineWidth}); hold on;
    else % response is stitched

        swind= winds(1);
        ewind= winds(2);

        time_axis= (0:(swind-1))/dataf - bef;   
        chl(1) = shadedErrorBar(time_axis,meanErp(1:swind),semErp(1:swind),...
            'lineprops', {'color', color, 'LineWidth', lineWidth}); hold on;

        time_axis= ((swind-1)+stitchgap:((swind-1)+stitchgap)+ewind-1)/dataf-bef;
        chl(2) = shadedErrorBar(time_axis,meanErp(end-ewind+1:end),semErp(end-ewind+1:end),...
            'lineprops', {'color', color, 'LineWidth', lineWidth}); hold on;

        % indicating a break in the ERP response
        xline(swind/dataf-bef, 'LineWidth',2, 'LineStyle',':');
        xline((swind+stitchgap)/dataf-bef, 'LineWidth',2, 'LineStyle',':');
        xline((sum(winds)+stitchgap)/dataf-2*bef, '-k', 'LineWidth', 2, ...
        'HandleVisibility','off');
        
    end

    set(gca, 'FontSize', 15);
    % Mark onset and offset
    % Assumes bef = aft
    xline(0, '-k', 'LineWidth',2, 'HandleVisibility','off');

    yline(0, '-k', 'LineWidth',1, 'HandleVisibility','off');
    xlabel('Time (s)');
    ylabel('HFA (z)')
    set(gca, 'FontSize', 15);
    box off;     
    hold on;
end


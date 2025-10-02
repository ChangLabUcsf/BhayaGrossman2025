function [] = drawBrainContour(ha, all_xyz, hemi, cmax, cols, plotbar)
    ax_pos = get(ha,'Position');
    yy_range = get(ha,'YLim');
    zz_range = get(ha,'ZLim');

    if nargin<5
        cols = internet;
    end

    if nargin<5
        plotbar=1;
    end

    % find density map of native speech, lateral side
    yye = min(all_xyz(:,2))-15:1:max(all_xyz(:,2)+15);
    zze = min(all_xyz(:,3))-15:1:max(all_xyz(:,3)+15); 
    ds_all = histcounts2(all_xyz(:,2),all_xyz(:,3),yye,zze);
    ds_all(isnan(ds_all)) = 0;
    
    % smooth the density map with a gaussian kernel
    gs_kernel = fspecial('gaussian', [15, 15], 3);
    ds_sm = conv2(ds_all,gs_kernel,'same');
    
    yy = yye(1)+diff(yye(1:2))/2:diff(yye(1:2)):yye(end);
    zz = zze(1)+diff(zze(1:2))/2:diff(zze(1:2)):zze(end);
    yyq = yy(1):0.2:yy(end);
    zzq = zz(1):0.2:zz(end); 
    ds_q = interp2(zz,yy,ds_sm,zzq,yyq','cubic');
    
    ha_ct = axes('Position',ax_pos);

    hold on
    [~,hc] = contourf(yyq,zzq,ds_q',10,'LineColor','none');
    %[~,hc] = contourf(yye(1:end-1),zze(1:end-1),ds_sm',15,'LineColor','none');
    
    clim([0 cmax])
    if strcmp(hemi, 'lh')
        set(ha_ct,'Color','None','XDir','reverse');
    else
        set(ha_ct,'Color','None');
    end
    axis equal
    axis off
    xlim(yy_range);
    ylim(zz_range);
    colormap(ha_ct,cols);
    drawnow;
    
    hFills = hc.FacePrims;  % array of TriangleStrip objects
    [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
    for i = 1:length(hFills) % changes alpha
        hFills(i).ColorData(4) = 230;   % default=255
    end
    hFills(1).ColorData(4) = 0;

    if plotbar
        figure; 
        colorbar;
        clim([0 cmax])
        colormap(cols);
        cbh = colorbar;
        cbh.Label.String = 'Density';
        cbh.Label.FontSize = 14;
    end
end

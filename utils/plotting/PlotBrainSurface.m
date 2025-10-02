function x_offset = PlotBrainSurface(cortex, hemi, side)
    % wrapper function for ctmr_gauss_plot
    % side = 'lateral', or 'medial'​
    ctmr_gauss_plot(cortex,[0 0 0], 0, hemi)

    if (strcmpi(side,'lateral') && strcmpi(hemi,'lh')) || ...
        (strcmpi(side,'medial') && strcmpi(hemi,'rh'))
        view_angle = [-90,0];
        x_offset = -5;
    elseif (strcmpi(side,'medial') && strcmpi(hemi,'lh')) || ...
            (strcmpi(side,'lateral') && strcmpi(hemi,'rh'))
        view_angle = [90,0];
        x_offset = 5;
    end
    view(view_angle)
    % add light
    if strcmpi(side,'lateral') && strcmpi(hemi,'lh')
        hl = light('Color',[1 1 1]*0.15);
        lightangle(hl, -90, -20)
        hl = light('Color',[1 1 1]*0.2);
        lightangle(hl, -30, 70)
        hl = light('Color',[1 1 1]*0.2);
        lightangle(hl, -150, 70)
        
    elseif strcmpi(side,'lateral') && strcmpi(hemi,'rh')
        hl = light('Color',[1 1 1]*0.2);
        lightangle(hl, 90, -20)
        hl = light('Color',[1 1 1]*0.2);
        lightangle(hl, 30, 70)
        hl = light('Color',[1 1 1]*0.2);
        lightangle(hl, 150, 70)
        
    elseif strcmpi(side,'medial') && strcmpi(hemi,'lh')
        hl = light('Color',[1 1 1]*0.8);
        lightangle(hl, 50, 30)
        hl2 = light('Color',[1 1 1]*0.6);
    %     lightangle(hl2, 130, -30)
        lightangle(hl2, 120, -30)
        
    elseif strcmpi(side,'medial') && strcmpi(hemi,'rh')
        hl = light('Color',[1 1 1]*0.8);
        lightangle(hl, -50, 30)
        hl2 = light('Color',[1 1 1]*0.6);
    %     lightangle(hl2, -130, -30)
        lightangle(hl2, -120, -30)
    end
end
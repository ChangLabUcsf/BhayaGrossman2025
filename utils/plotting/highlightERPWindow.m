function highlightERPWindow(wind, col)
    % reshape window for rectangle function

    % plot a rectangle that is the full height of current figure
    yl = ylim();

    h = rectangle('Position', ...
        [wind(1) yl(1) diff(wind([1 find(~isnan(wind), 1, 'last')])) diff(yl)], ...
        'FaceColor', col, 'EdgeColor',[1 1 1]); hold on;

    % make rectangle opaque
    h.FaceAlpha = 0.3;
end

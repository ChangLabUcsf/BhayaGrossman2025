function [en_mresp, sp_mresp] = plotStitchedSentence(elec, sent_encoding, swind, ewind, bef, overlaid)

    colors = flipud(brewermap(20-1, 'Spectral'));
    % spanish is 1, english is 2

    cols = [colors(3, :); colors(end-2, :); colors(round(size(colors, 1)/2), :)];
    %cols = [0.4 0.4 0.4; 0.8, 0.8, 0.8];
    
    range = 5;
    for el = elec
        
        % aggregate all sentences, find mean 
        en_mresp = sent_encoding.en_mresp{el};
        en_sem = sent_encoding.en_sem{el};
        sp_mresp = sent_encoding.sp_mresp{el};
        sp_sem = sent_encoding.sp_sem{el};

        figure('Renderer','painters');    
        gap = 20;
        xlim([-0.5 (ewind+swind+gap)./100-bef]);
        ylim([-0.5 4.5]);

        maxtp = sent_encoding.maxtp(el, :);
        wind=cell(1, 2);
        for lang = 1:2 
            if ~isnan(maxtp(lang))
                tmpwind = maxtp(lang)-range:maxtp(lang)+range; 
                if maxtp(lang)>swind, tmpwind=tmpwind+gap; end
                wind{lang} = tmpwind; 
            end 
        end
    
        if overlaid
            % Plot the response traces
            sentenceERP(sp_mresp, sp_sem, cols(1, :), bef, ...
                [swind, ewind], gap); hold on; %[0.19 0.53 0.74]
            if ~isempty(wind{2})
                highlightERPWindow(wind{2}./100-bef, cols(1, :));
            end
        
            sentenceERP(en_mresp, en_sem, cols(2, :), bef, ...
                [swind, ewind], gap); %[0.83 0.24 0.30]
            if ~isempty(wind{1})
                highlightERPWindow(wind{1}./100-bef, cols(2, :));
            end

            % Calculate significant differences per time point
            fvals = sent_encoding.fvals{el};
            fthresh = sent_encoding.fthresh(el);
            dataf = 100;
            aftergap = (swind+gap)/100+(1/dataf);
            x = [-1*bef+(1/dataf):(1/dataf):(swind/dataf)-bef ...
                (aftergap:(1/dataf):aftergap+ewind/dataf)-bef];
    
            scatter(x(fvals>fthresh), -0.1*ones(1, sum(fvals>fthresh)), 15, ...
            [fvals(fvals>fthresh)]', 'filled');
            colormap(flipud(gray));
        else
            subplot(2, 1, 1);
            % Plot the response traces
            sentenceERP(sp_mresp, sp_sem, cols(1, :), bef, ...
                [swind, ewind], gap); hold on; %[0.19 0.53 0.74]
            % if ~isempty(wind{2})
            %     highlightERPWindow(wind{2}./100-bef, cols(1, :));
            % end
        
            subplot(2, 1, 2);
            sentenceERP(en_mresp, en_sem, cols(2, :), bef, ...
                [swind, ewind], gap); %[0.83 0.24 0.30]
            % if ~isempty(wind{1})
            %     highlightERPWindow(wind{1}./100-bef, cols(2, :));
            % end
        end

        sgtitle(['Language exp: ' num2str(sent_encoding.ls(el)) ...
            ', el: '  num2str(sent_encoding.el(el)) ' SID: ' ...
            sent_encoding.SID{el}]);     
    end

end


function [colors] = getColorsCrossComp(type)
    switch type
        case 1 % language type
            % spanish, english, bilingual mandarin
            colors = brewermap(4, 'Dark2');
            colors(4, :) = [35 170 225]./256;
        case 2 % primary feature encoding map
            colors = brewermap(3, 'Dark2');
        case 3 % word syllable distinction
%             colors = [0.6 0 0.6; 0.1 0.7 0.2];
            colors = [0.3 0.3 0.3; 0.6 0 0.6];
        case 4 % VOT
            colors = nuclear(8);
            colors = colors([2, 7], :);
        case 5
            colors = nuclear(9);
            colors = colors(3:8, :);
        case 6 % bilingual colors
            colors = [34 174 209; 43, 57, 144]./256; % Spanish / English
    end
end

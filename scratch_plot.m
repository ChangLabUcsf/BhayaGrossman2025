el = 88; % 70
SID = 'EC183';
addpath brewer
ctr = 1;
figure;
label = {'Spanish (DIMEX)', 'English (TIMIT)'};
for c = {DDcons, TDcons}

    consStruct = c{1};
    resp = consStruct.(SID).resp;
    onsIdx = consStruct.wordOns==1;
    syllOns = consStruct.syllOns==1;
    consStruct.sentOns = [1 diff(consStruct.sent)>0];
    xs = -0.2:0.01:0.6;

    subplot(2, 2, ctr)    
    nantrl = squeeze(sum(isnan(resp(el, :, :))))>0;
    r = squeeze(resp(el, :, onsIdx'&~nantrl))';
    if ctr == 3
        consStruct.wrdLen = cellfun(@(x) length(x),consStruct.wrd);
        [~, idx] = sort(consStruct.wrdLen(onsIdx'&~nantrl), 'ascend');
        r = r(idx, :);
    end    
    imagesc(xs, 1:sum(onsIdx'&~nantrl), r); 
    colormap(flipud(brewermap(20, 'RdBu')));
    ylabel('Trials');
    yticks([]);
    set(gca, 'FontSize', 15);
    xline(0, 'LineWidth',2);
    title(label{floor(ctr/2)+1});

    subplot(2, 2, ctr+1)   
    % word onsets
    shadedErrorBar(xs, squeeze(mean(resp(el, :, onsIdx & ~consStruct.sentOns), ...
        3, 'omitnan')), squeeze(nansem(resp(el, :, onsIdx & ~consStruct.sentOns), 3)), ...
        {'-', 'LineWidth', 1.25, 'Color', [0.5 0.2, 0.6], 'HandleVisibility', 'off'}); hold on;
    % phoneme onsets
    shadedErrorBar(xs, squeeze(mean(resp(el, :,~onsIdx), 3, 'omitnan')),...
        squeeze(nansem(resp(el, :,~onsIdx), 3)), ...
        {'-', 'LineWidth', 1.25, 'Color', [0.3 0.7, 0.3], 'HandleVisibility', 'off'}); hold on;
    % syllable onsets
    if ctr == 3
        shadedErrorBar(xs, squeeze(mean(resp(el, :,~onsIdx & syllOns), 3, 'omitnan')),...
            squeeze(nansem(resp(el, :,~onsIdx & syllOns), 3)), ...
            {'-', 'LineWidth', 1.25, 'Color', [0.7 0.7, 0.3], 'HandleVisibility', 'off'});
    end
    xline(0, 'LineWidth',2);
    ylabel('HGA')
    set(gca, 'FontSize', 15);
    ylim([0.2 1])
   
    ctr = ctr+2;
end
sgtitle(['word onset ERP: ' SID ' , el = ' num2str(el)]);
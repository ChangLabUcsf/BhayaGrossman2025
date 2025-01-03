figure;
sent = 20;
subplot(2, 1, 1);
imagesc(flipud(timit_details.sentdet(sent).aud));
cols = colormap('bone');
colormap(flipud(cols));
xline(find(timit_details.sentdet(sent).wordOns), 'LineWidth', 2);
xticks([]);
subplot(2, 1, 2)
plot(timit_details.sentdet(sent).sound, 'Color',[0.1 0.1 0.1]); box off
dataf = timit_details.sentdet(sent).dataf;
soundf = timit_details.sentdet(sent).soundf;
xline((find(timit_details.sentdet(sent).wordOns)./dataf)*soundf, 'LineWidth', 2);
xlim([0 length(timit_details.sentdet(sent).sound)]);
xticks([]);
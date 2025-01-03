%% read in text grid
% Phoneme Times
phonemeTimes = [
    0.028167497985877272,
    0.1181488106824311,
    0.2157185624319476,
    0.2993317169204928,
    0.35987476589322226,
    0.43076743364757886,
    0.5410038417443638,
    0.6403349834761254,
    0.7142897125038924,
    0.8371131884978156,
    0.9569087126517216,
    1.0433702566835625,
    1.2423024620347953,
    1.3860300472333118,
    1.427737385146842,
    1.507824236204427,
    1.5815395311213643,
    1.7075796176082412,
    1.7929341696173264,
    1.9245616565967,
    1.9914873848765513,
    2.084601441613735,
    2.1410312583517106,
    2.2549217861218027,
    2.3025586895641013,
    2.5092036637779347,
    2.6376550448495237,
    2.695537457307709
];

% Syllable Times
syllableTimes = [
    0.028167497985877272,
    0.2157185624319476,
    0.35987476589322226,
    0.5410038417443638,
    0.7142897125038924,
    0.9569087126517216,
    1.3860300472333118,
    1.507824236204427,
    1.7075796176082412,
    1.9245616565967,
    2.084601441613735,
    2.2549217861218027,
    2.4496869200518376,
    2.6376550448495237
];

% Word Times
wordTimes = [
    0.028167497985877272,
    0.2157185624319476,
    0.7142897125038924,
    1.507824236204427,
    1.7075796176082412,
    2.084601441613735,
    2.2549217861218027,
    2.4496869200518376
];

% Phonemes
phonemes = {
    'l', 'a', 'm', 'o', 'n', 'a', 'x', 'o', 's', 'i', 'd', 'i', 's', 'ch', 'e', 'l', 'a', 'b', 'a', 'z', 'o', 'd', 'e', 'l', 'a', 'r', 'b', 'o'
};

words = {'la' 'monaxo', 'sidis', 'che', 'la', 'bazo', 'de' 'la' 'arbo'};

%%  audio file for esperanto

% read in text grid
phones = {''};
phntimes = [];

syltimes = [];
wordtimes = [];

% read in audio file
[esperanto,fs] = audioread('esperanto_audio.mp3'); % read in audio file

% generate the spectrogram
t1 = resample(esperanto,16000,fs);
pspectrum = powspec(t1, 16000);
aspectrum = audspec(pspectrum, 16000, 80, 'mel');
fs = 16000;
aud = aspectrum.^( 0.077);

figure; 
subplot(4, 1, 1);
% get timepoints for the sounds using fs
t = 0:1/fs:(length(t1)-1)/fs; % get timepoints for the sounds using fs
plot(t, t1, 'Color', [0.5 0.5 0.5]); % plot the audio file
% generate and plot the speech envelope
env = envelope(t1, 500, 'peak');
hold on;
plot(t, env, 'k', 'LineWidth', 1.5); % plot the speech envelope
set(gca, 'Ydir', 'normal', 'Fontsize', 13); % set the y-axis to be normal and set the font size

subplot(4, 1, 2);
imagesc(aud, 'Xdata', [0 max(t)], 'Ydata', [0 8]); % plot the spectrogram
ylabel('Frequency (Hz)'); % y-axis label
yticks([1, 80]);
xticks([]);
colormap(flipud(gray));
yticklabels({'0', '8'});
set(gca, 'Ydir', 'normal', ...
    'Fontsize', 13); % set the y-axis to be normal and set the font size

% plot text grid
subplot(4, 1, 3);
xline(phonemeTimes);
xticks([]);
% write out the phonemes
for i = 1:length(phonemeTimes)
    text(phonemeTimes(i), 0.5, phonemes{i}, 'HorizontalAlignment', ...
        'left', 'Fontsize', 13);
end
set(gca, 'Ydir', 'normal', ...
    'Fontsize', 13);

subplot(4, 1, 4);

% syllables times as dashes and words as vertical lines
xline(wordTimes);
xticks(0:3);
xlabel('Time (s)')
xline(syllableTimes, '--', 'Color', [0.5 0.5 0.5]);
% write out the words
for i = 1:length(wordTimes)
    text(wordTimes(i), 0.5, words{i}, 'HorizontalAlignment', 'left', ...
        'Fontsize', 13);
end
set(gca, 'Ydir', 'normal', ...
    'Fontsize', 13);

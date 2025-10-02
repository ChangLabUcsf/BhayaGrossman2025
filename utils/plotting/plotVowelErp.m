%% plot erps over vowels
function plotVowelErp(Dvow, SID, singleSet, corpus, f, cols, ...
    unique_vowels)
    % dimex
    addpath brewer
    desmat.names = [];
    desmat.names = unique(Dvow.vowel);
    desmat.condition = Dvow.vowelType;
    %unique_vowels = {'e'}; % five in spanish
    if nargin<6
        cols = brewermap(5, 'Dark2'); 
    end
    
    if nargin<7
        switch corpus
            case 'dimex'
                unique_vowels = {'a', 'e', 'i', 'o', 'u'};
            case 'timit'
                unique_vowels = {'aa', 'ey', 'iy', 'ow', 'uw'};
        end
    end
      
    vowelNameIdx=ismember(unique(Dvow.vowel), unique_vowels);
    names = unique(Dvow.vowel);
    desmat.names = names(vowelNameIdx);

    % all indices with vowels to keep
    vowelCorpusIdx = ismember(Dvow.vowelType,find(vowelNameIdx));

    % relabel vowels in order of how they appear in vowel names
    desmat.condition = nan(1, sum(vowelCorpusIdx, 'omitnan'));
    for i=1:length(find(vowelNameIdx))
        vowelInds = find(vowelNameIdx);
        desmat.condition(1, Dvow.vowelType(vowelCorpusIdx)==vowelInds(i))=i;
    end
    desmat.vowelNames = Dvow.vowel(vowelCorpusIdx);

    %plot response
    % remove util shadedErrorBar folder from path!
    plotResp(Dvow.(SID).resp(:, :, vowelCorpusIdx), desmat, singleSet, SID, f, cols);
end
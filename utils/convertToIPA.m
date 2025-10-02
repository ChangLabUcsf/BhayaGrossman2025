function ipaLabels = convertToIPA(mapping, arpabetLabels)
   
    % Initialize the output cell array
    ipaLabels = cell(size(arpabetLabels));

    % Iterate over each label in the input cell array
    for i = 1:numel(arpabetLabels)
        arpabetLabel = arpabetLabels{i};

        % Find the corresponding ipa label in the mapping
        mappingIndex = find(strcmp(mapping(:, 1), arpabetLabel), 1);

        % If a mapping is found, assign the ipa label to the output cell array
        if ~isempty(mappingIndex)
            ipaLabel = mapping{mappingIndex, 2};
            ipaLabels{i} = ipaLabel;
        else
            % If no mapping is found, keep the original arpabet label
            ipaLabels{i} = arpabetLabel;
        end
    end
end
function [phenotypeVector,ia] = MatchMeCellDensity(structInfo,theCellType)

switch theCellType
case {'VIP','SST','PV'}
    % Interneuron cell density measures from our friends, Kim et al.:
    dataOutput = ImportInterneuronCellDensities();
    regionNames = dataOutput.name;
    theCellType = sprintf('%s_mean',theCellType);
otherwise
    % Cell density measures from Ero et al.:
    dataOutput = ImportCellAtlas('density');
    regionNames = dataOutput.regionName;
end

% Match names to acronyms:
structInfoNames = regexprep(structInfo.name,',','');
[~,ia,ib] = intersect(lower(structInfoNames),lower(regionNames));
dataOutput = dataOutput(ib,:);
fprintf(1,'%u names match to set of %u acronyms\n',...
                height(dataOutput),height(structInfo));
phenotypeVector = dataOutput.(theCellType);

end

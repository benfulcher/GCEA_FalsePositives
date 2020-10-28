function [phenotypeVector,ia,regionNamesOrdered] = MatchMeCellDensity(structInfo,theCellType)

switch theCellType
case {'VIP','SST','PV'}
    % Interneuron cell density measures from our friends, Kim et al.:
    dataOutput = ImportInterneuronCellDensities();
    regionNames_cell = dataOutput.acronym;
    regionNames_ref = regexprep(structInfo.acronym,',','');
    theCellType = sprintf('%s_mean',theCellType);
otherwise
    % Cell density measures from Ero et al.:
    dataOutput = ImportCellAtlas('density');
    regionNames_cell = dataOutput.regionName;
    regionNames_ref = regexprep(structInfo.name,',','');
end

% Match names to acronyms:
[regionNamesOrdered,ia,ib] = intersect(lower(regionNames_ref),lower(regionNames_cell));
dataOutput = dataOutput(ib,:);
phenotypeVector = dataOutput.(theCellType);

% Output to commandline:
fprintf(1,'%u names match to set of %u acronyms\n',...
                height(dataOutput),height(structInfo));

end

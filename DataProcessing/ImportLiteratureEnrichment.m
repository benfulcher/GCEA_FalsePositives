function [resultsTables,mouseOrHuman] = ImportLiteratureEnrichment(filterOnOurGenes,doSave)
% Loads in a table of enrichment results for all prior studies using GO
% enrichment
%-------------------------------------------------------------------------------

% Whether to only look at categories with annotations for genes in our set:
if nargin < 1
    filterOnOurGenes = false;
end
if nargin < 2
    doSave = true;
end

%-------------------------------------------------------------------------------
% Load in data from the file in which I've manually annotated enriched GO terms
% From a range of published studies:
manualDataFile = 'TableGOBPs.csv';
TableGOBPs = readtable(manualDataFile);
% Split off species identifier:
% speciesID = TableGOBPs(1,3:end);
% TableGOBPs = TableGOBPs(3:end,:);
% Give user feedback:
numManual = height(TableGOBPs);
numStudies = length(TableGOBPs.Properties.VariableNames)-2;
fprintf(1,'Loaded manual annotations for %u GO Terms across %u studies from %s\n',...
            numManual,numStudies,manualDataFile);

%-------------------------------------------------------------------------------
% Load in data on GO category annotations:
params = GiveMeDefaultParams();
params.e.sizeFilter = [1,1e6];
if filterOnOurGenes
    % (Only look at categories with annotations for genes in our set)
    [~,geneInfo] = LoadMeG(params.g);
    restrictEntrez = geneInfo.entrez_id;
else
    restrictEntrez = [];
end
GOTerms = GiveMeGOData(params,restrictEntrez);
numGOCategories = height(GOTerms);

%-------------------------------------------------------------------------------
% Try to match GO categories by name:
fprintf(1,'Matching %u manually-compiled GO categories to the full list',...
                        ' of GO BPs by name\n',numManual);
matchMe = nan(numManual,1);
GOIDsMatches = nan(numManual,1);
for i = 1:numManual
    itsHere = find(strcmp(TableGOBPs.GOCategory{i},GOTerms.GOName));
    if isempty(itsHere)
        warning('There is no GO category called: ''%s''',TableGOBPs.GOCategory{i})
    else
        matchMe(i) = itsHere;
        GOIDsMatches(i) = GOTerms.GOID(matchMe(i));
    end
end

%-------------------------------------------------------------------------------
% Convert each entry to an element in resultsTables (of the form GOID [numeric], p-value)
resultsTables = struct();

theManualResultNames = TableGOBPs.Properties.VariableNames;
% exclude GOCategory, ID:
theManualResultNames = setxor(theManualResultNames,{'GOCategory','ID'});
numManualResults = length(theManualResultNames);
fprintf(1,'%u manual result tables\n',numManualResults);

for i = 1:numManualResults
    theData = TableGOBPs.(theManualResultNames{i});

    if isnumeric(theData)
        hasHits = ~isnan(theData);
    else
        hasHits = ~cellfun(@isempty,theData);
    end

    theGOIDs = GOIDsMatches(hasHits);
    isValid = ~isnan(theGOIDs);
    GOID = theGOIDs(isValid);

    theData = theData(hasHits);

    if isnumeric(theData)
        hasPVals = any(theData < 1);
    else
        try
            theData = cellfun(@(x)str2num(x),theData);
            hasPVals = any(theData < 1);
        catch
            fprintf(1,'Non-numeric data for %s---treating all marked categories as significant\n',...
                                theManualResultNames{i});
            hasPVals = false;
        end
    end

    if hasPVals
        pValCorr = theData(isValid);
        resultsTables.(theManualResultNames{i}) = table(GOID,pValCorr);
    else
        pValCorr = zeros(sum(isValid),1);
        resultsTables.(theManualResultNames{i}) = table(GOID,pValCorr);
    end
end

%-------------------------------------------------------------------------------
% Now we need to load in the quantitative data:
%-------------------------------------------------------------------------------
% In each case, we want GOID, p-value

% --- Forest2017-TableS8-PathwayEnrichment_ReducedModel.xlsx
% (GOID (string) + corrected p-value)
resultsTables.ForestReduced = ImportForestReduced();

% --- Forest2017_TableS3-PathwayEnrichment_FullModel.xlsx
resultsTables.ForestFull = ImportForestFullModel();

% ---French2011 (nothing much here -- different enrichment for NE and OL)
% -(neurons vs oligodendrocytes)-
[resultsTables.French2011NE,resultsTables.French2011OE] = ImportFrench2011();

% ---French2015-ConsistentGOGroups.csv
% (GO-IDs as numeric)
resultsTables.French2015Consistent = ImportFrench2015Consistent();

% ---French2015-InconsistentGOGroups.csv
% (GO-IDs as numeric)
resultsTables.French2015Inconsistent = ImportFrench2015Inconsistent();

% ---Parkes et al.:
PCs = [1,2,5,9];
numPCs = length(PCs);
GOtoNumber = @(x)str2num(x(4:end));
for i = 1:numPCs
    fileName = sprintf('Parkes2017_PC%u.txt',PCs(i));
    ResultsTable = ReadInErmineJ(fileName);
    GOID = cellfun(GOtoNumber,ResultsTable.GOID);
    pValCorr = ResultsTable.pValCorr;
    resultsTables.(sprintf('ParkesPC%u',PCs(i))) = table(GOID,pValCorr);
end

% ---Tan2013-table-s6-david-200pos-transport.csv
% (results from DAVID)
resultsTables.Tan2013 = ImportTan2013();

% ---Vertes-rstb20150362supp1.xlsx
% (note that Vertes actually excluded many categories using exclude column,
% but we keep all) [excluded are very general terms, >1000 gene annotations,
% or those considered 'redundant'?]
resultsTables.Vertes2015_PLS1pos = ImportVertes2015('PLS1 pos');
resultsTables.Vertes2015_PLS2pos = ImportVertes2015('PLS2 pos');
resultsTables.Vertes2015_PLS3pos = ImportVertes2015('PLS3 pos');
resultsTables.Vertes2015_PLS1neg = ImportVertes2015('PLS1 neg');
resultsTables.Vertes2015_PLS2neg = ImportVertes2015('PLS2 neg');
resultsTables.Vertes2015_PLS3neg = ImportVertes2015('PLS3 neg');

% ---Whitaker:
resultsTables.WhitakerCompletePLS2pos = ImportWhitaker('Complete_PLS2pos');
resultsTables.WhitakerCompletePLS2neg = ImportWhitaker('Complete_PLS2neg');
resultsTables.WhitakerDiscoveryPLS2pos = ImportWhitaker('Discovery_PLS2pos');
resultsTables.WhitakerDiscoveryPLS2neg = ImportWhitaker('Discovery_PLS2neg');
resultsTables.WhitakerValidationPLS2pos = ImportWhitaker('Validation_PLS2pos');
resultsTables.WhitakerValidationPLS2neg = ImportWhitaker('Validation_PLS2neg');

% ---Fulcher:
% (connected/unconnected comparison) [done at a pairwise level: 'either']:
resultsTables.Fulcher2016conn = readtable('Fulcher2016_connectedUnconnected_BP_TableS1.csv');
% (rich+feed versus peripheral comparison) [done at a pairwise level]:
resultsTables.Fulcher2016rich = readtable('Fulcher2016_richFeederPeripheral_BP_TableS5.csv');

% ---Meijer:
resultsTables.Meijer2019stress = ImportMeijer();

%-------------------------------------------------------------------------------
% Ok, so now we have resultsTables from all data combined :-D
% Now we need to separate into human and mouse studies
% Read in annotations of studies:
fid = fopen('AnnotateStudies.csv');
S = textscan(fid,'%s%s');
fclose(fid);
% Match to names:
allTableNames = fieldnames(resultsTables);
mouseOrHuman = cell(length(allTableNames),1);
for i = 1:length(allTableNames)
    isHere = strcmp(S{1},allTableNames{i});
    if ~any(isHere)
        warning('No annotation found for %s',allTableNames{i})
        keyboard
        mouseOrHuman(i) = NaN;
    else
        mouseOrHuman(i) = S{2}(isHere);
    end
end
mouseOrHuman = categorical(mouseOrHuman);

%-------------------------------------------------------------------------------
% Save:
fileNameSave = fullfile('DataOutputs','LiteratureEnrichmentLoaded.mat');
if doSave
    save(fileNameSave,'resultsTables','mouseOrHuman');
end
fprintf(1,'Saved to %s\n',fileNameSave);


end

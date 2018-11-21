function [resultsTables,mouseOrHuman] = ImportLiteratureEnrichment(filterOnOurGenes,doSave)
% Loads in a table of enrichment results for all prior studies using GO
% enrichment
%-------------------------------------------------------------------------------

% Whether to only look at categories with annotations for genes in our set:
if nargin < 1
    filterOnOurGenes = true;
end
if nargin < 2
    doSave = true;
end

%-------------------------------------------------------------------------------
% Load in the file in which I've manually annotated GO terms:
manualDataFile = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/TableGOBPs.csv';
TableGOBPs = ImportGOBPs(manualDataFile);
emptyGO = cellfun(@isempty,TableGOBPs.GOCategory);
TableGOBPs = TableGOBPs(~emptyGO,:);
numManual = height(TableGOBPs);
fprintf(1,'Loaded manual annotations for %u GO Terms from %s\n',numManual,manualDataFile);

%===============================================================================
% MATCH CATEGORIES TO OUR BP INFORMATION
%===============================================================================
% Load in general gene info data (we just need it for the entrez_ids)
params = GiveMeDefaultParams();
params.e.sizeFilter = [1,1e5];
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

% ---Load in annoated GO Table---
% (Only look at categories with annotations for genes in our set)
if filterOnOurGenes
    restrictEntrez = geneInfo.entrez_id;
else
    restrictEntrez = [];
end
GOTerms = GetFilteredGOData(params.e.dataSource,params.e.processFilter,...
                                params.e.sizeFilter,restrictEntrez);
numGOCategories = height(GOTerms);

%-------------------------------------------------------------------------------
% Try to match GO categories by name:
fprintf(1,'Matching %u manually-compiled GO categories to full list of GO BPs by name\n',numManual);
matchMe = nan(numManual,1);
GOIDsMatches = nan(numManual,1);
for i = 1:numManual
    itsHere = find(strcmp(TableGOBPs.GOCategory{i},GOTerms.GOName));
    if isempty(itsHere)
        warning('No match found for ''%s''',TableGOBPs.GOCategory{i})
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
% Ignore some results that are elsewhere:
ignoreThese = {'ParkesPC1','ParkesPC2','ParkesPC5','ParkesPC9'};
theManualResultNames = setxor(theManualResultNames,ignoreThese);
numManualResults = length(theManualResultNames);
fprintf(1,'%u manual result tables\n',numManualResults);

for i = 1:numManualResults
    theData = TableGOBPs.(theManualResultNames{i});
    hasHits = ~isnan(theData);
    theGOIDs = GOIDsMatches(hasHits);
    theData = theData(hasHits);

    isValid = ~isnan(theGOIDs);
    GOID = theGOIDs(isValid);
    pValCorr = theData(isValid);

    hasPVals = any(theData < 1);
    if hasPVals
        resultsTables.(theManualResultNames{i}) = table(GOID,pValCorr);
    else
        pValCorr = zeros(size(pValCorr));
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
    fileName = sprintf('/Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Parkes2017/PC%u.txt',...
                        PCs(i));
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
% (rich+feed versus peripheral comparison) [also done at a pairwise level]:
[resultsTables.FulcherConnected,~,resultsTables.FulcherRichFeed] = ImportFulcher2016();

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
        mouseOrHuman(i) = NaN;
    else
        mouseOrHuman(i) = S{2}(isHere);
    end
end
mouseOrHuman = categorical(mouseOrHuman);

%-------------------------------------------------------------------------------
% Save?
fileNameSave = fullfile('DataOutputs','LiteratureEnrichmentLoaded.mat');
if doSave
    save(fileNameSave,'resultsTables','mouseOrHuman');
end
fprintf(1,'Saved to %s\n',fileNameSave);

%-------------------------------------------------------------------------------
function TableGOBP = ImportGOBPs(filename, startRow, endRow)
%IMPORTFILE Import numeric data from a text file as a matrix.
%   TABLEGOBPS1 = IMPORTFILE(FILENAME) Reads data from text file FILENAME
%   for the default selection.
%
%   TABLEGOBPS1 = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data from
%   rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   TableGOBP = importfile('TableGOBPs.csv', 2, 181);
%
%    See also TEXTSCAN.

% Auto-generated by MATLAB on 2017/10/27 16:06:22

%% Initialize variables.
delimiter = ',';
if nargin<=2
    startRow = 2;
    endRow = inf;
end

%% Read columns of data as text:
% For more information, see the TEXTSCAN documentation.
formatSpec = '%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%q%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

%% Close the text file.
fclose(fileID);

%% Convert the contents of columns containing numeric text to numbers.
% Replace non-numeric text with NaN.
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));

for col=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]
    % Converts text in the input cell array to numbers. Replaced non-numeric
    % text with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1)
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData(row), regexstr, 'names');
            numbers = result.numbers;

            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if numbers.contains(',')
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'))
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric text to numbers.
            if ~invalidThousandsSeparator
                numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch
            raw{row, col} = rawData{row};
        end
    end
end


%% Split data into numeric and string columns.
rawNumericColumns = raw(:, [2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34]);
rawStringColumns = string(raw(:, 1));


%% Replace non-numeric cells with NaN
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),rawNumericColumns); % Find non-numeric cells
rawNumericColumns(R) = {NaN}; % Replace non-numeric cells

%% Create output variable
TableGOBP = table;
TableGOBP.GOCategory = rawStringColumns(:, 1);
TableGOBP.ID = cell2mat(rawNumericColumns(:, 1));
TableGOBP.Forest2017blue = cell2mat(rawNumericColumns(:, 2));
TableGOBP.Forest2017green = cell2mat(rawNumericColumns(:, 3));
TableGOBP.Forest2017pink = cell2mat(rawNumericColumns(:, 4));
TableGOBP.Liu2017bioRxivChronicSZ = cell2mat(rawNumericColumns(:, 5));
TableGOBP.Liu2017bioRxivAutism = cell2mat(rawNumericColumns(:, 6));
TableGOBP.Liu2017bioRxivMiddleInsularAreaHCP = cell2mat(rawNumericColumns(:, 7));
TableGOBP.Liu2017bioRxivLeftHippocampusHCP = cell2mat(rawNumericColumns(:, 8));
TableGOBP.KunchevaCluster1REVIGOTop = cell2mat(rawNumericColumns(:, 9));
TableGOBP.KunchevaCluster2REVIGOTop = cell2mat(rawNumericColumns(:, 10));
TableGOBP.KunchevaCluster3REVIGOTop1InTopList2 = cell2mat(rawNumericColumns(:, 11));
TableGOBP.Rommetop100ConsensusPathDB = cell2mat(rawNumericColumns(:, 12));
TableGOBP.Rommetop540BFPANTHERrawPisSig = cell2mat(rawNumericColumns(:, 13));
TableGOBP.Rommetop200PANTHERrawp = cell2mat(rawNumericColumns(:, 14));
TableGOBP.Rommetop100PANTHERrawp = cell2mat(rawNumericColumns(:, 15));
TableGOBP.ParkesPC1 = cell2mat(rawNumericColumns(:, 16));
TableGOBP.ParkesPC2 = cell2mat(rawNumericColumns(:, 17));
TableGOBP.ParkesPC5 = cell2mat(rawNumericColumns(:, 18));
TableGOBP.ParkesPC9 = cell2mat(rawNumericColumns(:, 19));
TableGOBP.Mills2017bioRxivBPNoDistCorr = cell2mat(rawNumericColumns(:, 20));
TableGOBP.Mills2017bioRxivDistCorr = cell2mat(rawNumericColumns(:, 21));
TableGOBP.Ji2015primaryInjectionDrivenPID = cell2mat(rawNumericColumns(:, 22));
TableGOBP.Ji2015allInjectionDrivenAID = cell2mat(rawNumericColumns(:, 23));
TableGOBP.Fahkry2015 = cell2mat(rawNumericColumns(:, 24));
TableGOBP.Rubinov2015S3topcategories = cell2mat(rawNumericColumns(:, 25));
TableGOBP.Rubinov2015S4bottom25 = cell2mat(rawNumericColumns(:, 26));
TableGOBP.RichiardiPanther = cell2mat(rawNumericColumns(:, 27));
TableGOBP.TanPositive100isSig = cell2mat(rawNumericColumns(:, 28));
TableGOBP.TanNegative100isSig = cell2mat(rawNumericColumns(:, 29));
TableGOBP.French2011outgoing = cell2mat(rawNumericColumns(:, 30));
TableGOBP.French2011incoming = cell2mat(rawNumericColumns(:, 31));
% TableGOBP.French2015inconsistent = cell2mat(rawNumericColumns(:, 32));
% TableGOBP.French2011FrontN = cell2mat(rawNumericColumns(:, 33));

end

end

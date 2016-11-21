function writeErmineJFile(whatData,geneMeasures,theGeneEntrez,columnName)
% Writes a tab-delimited input file for ermineJ
%
% Ben Fulcher, 2014-10-14
% ------------------------------------------------------------------------------

if nargin < 4
    columnName = 'geneMeanP';
end

%
% whatData = 'machineLearn';
%
% % Open the file
% switch whatData
% case 'machineLearn'
%     fileName = 'ermineJInputFile_massU.txt';
%     % geneMeasures = pvalues(:,6); %meanP;
%     % theGeneEntrez = {theGeneStruct.gene_acronym};
% case 'meanEnrich'
%     fileName = 'ermineJInputFile_Enrich.txt';
%     % geneMeasures = meanP;
%     % theGeneEntrez = {theGeneStruct.filtered.gene_acronym};
% end

fileName = sprintf('ermineJInputFile_%s.txt',whatData);

fid = fopen(fileName,'w');

numGenes = length(theGeneEntrez);

% Case we have meanP and theGeneStruct
% ------------------------------------------------------------------------------

% 1. The header
fprintf(fid,'%s\t%s\n','Gene',columnName);

% 2. The gene list with p-values.
for i = 1:numGenes
    fprintf(fid,'%u\t%.8g\n',theGeneEntrez(i),geneMeasures(i));
    % fprintf(fid,'%s\t%f\n',theGeneStruct(i).gene_acronym,meanP(i));
end

fclose(fid);

% ------------------------------------------------------------------------------
% Display result:

fprintf(1,'\nWrote %s to file for ermineJ.\n\n',fileName);

end

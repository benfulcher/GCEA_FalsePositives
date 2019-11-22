whatSpecies = 'mouse';
structFilter = 'none';
numMaps = 10000;
whatSurrogate = 'spatialLag';
customSurrogate = '';

%-------------------------------------------------------------------------------
% Load in surrogate enrichment results:
fileNameIn = sprintf('SurrogateGOTables_%u_%s_%s_%s.mat',numMaps,whatSpecies,whatSurrogate,customSurrogate);
load(fileNameIn,'GOTableGeneric','surrogatePVals');
% surrogatePVals is the GO category x random map matrix

%-------------------------------------------------------------------------------
% Load in the random maps used:
% (from ComputeAllCategoryNulls; should be same as used for SurrogateEnrichment,
% loading from LoadMeG)
switch whatSpecies
case 'mouse'
    if strcmp(structFilter,'cortex')
        fprintf(1,'Spatial maps for mouse cortex\n');
        dataFileSurrogate = 'mouseCortexSurrogate_N20000_rho8_d040.csv';
    else
        fprintf(1,'Spatial maps for mouse whole brain\n');
        dataFileSurrogate = 'mouseSurrogate_N20000_rho8_d040.csv';
    end
case 'human'
    fprintf(1,'Spatial maps for human cortex\n');
    dataFileSurrogate = 'humanSurrogate_N20000_rho8_d02000.csv';
end
nullMaps = dlmread(dataFileSurrogate,',',1,1);

%-------------------------------------------------------------------------------
% Get real gene-expression data (to compute correlation scores with)
params = GiveMeDefaultParams(whatSpecies);
[geneDataReal,geneInfoReal,structInfoReal] = LoadMeG(params.g);
numGenesReal = height(geneInfoReal);

%-------------------------------------------------------------------------------
% Compute corrected p-vals:
pValCorr = zeros(size(surrogatePVals));
for j = 1:size(surrogatePVals,2)
    pValCorr(:,j) = mafdr(surrogatePVals(:,j),'BHFDR','true');
end

%-------------------------------------------------------------------------------
% Now we can count number of significant categories for each random map:
numSigFDR = sum(pValCorr < 0.05);
[~,ix] = sort(numSigFDR,'descend');

% Output some examples of the top-significance ones:
numToOutput = 10;
baseName = 'mouse_spatialLag_map_';
for i = 1:numToOutput
    index = ix(i);
    theMap = nullMaps(:,index);

    % Compute gene scores (from SurrogateEnrichment)
    geneScores = zeros(numGenesReal,1);
    for j = 1:numGenesReal
        geneScores(j) = corr(theMap,geneDataReal(:,j),'type','Spearman','rows','pairwise');
    end

    % Save the gene scores:
    fid = fopen(sprintf('%s%u_geneScores.csv',baseName,index),'w','n');
    for j = 1:numGenesReal
        fprintf(fid,'%u, %g\n',geneInfoReal.entrez_id(j),geneScores(j));
    end
    fclose(fid);

    % Save the GO significance summary (FDR = 5%):
    fid_FDR = fopen(sprintf('%s%u_GOsig.csv',baseName,index),'w','n');
    % Save the GO significance summary (raw p-value < 0.05):
    fid_RAW = fopen(sprintf('%s%u_GOsigRAW.csv',baseName,index),'w','n');

    isSigFDR = find(pValCorr(:,index) < 0.05);
    isSigRAW = find(surrogatePVals(:,index) < 0.05);

    sigCategoriesFDR = length(isSigFDR);
    sigCategoriesRAW = length(isSigRAW);
    fprintf(1,'%u FDR-significant, %u RAW-significant\n',sigCategoriesFDR,sigCategoriesRAW);

    % FDR:
    sigPValFDR = pValCorr(isSigFDR,index);
    [~,ip] = sort(sigPValFDR,'ascend');
    isSigFDR = isSigFDR(ip);
    for j = 1:sigCategoriesFDR
        fprintf(fid_FDR,'%s | %g\n',GOTableGeneric.GOName{isSigFDR(j)},pValCorr(isSigFDR(j),index));
    end
    fclose(fid_FDR);

    % RAW:
    sigPVal = surrogatePVals(isSigRAW,index);
    [~,ip] = sort(sigPVal,'ascend');
    isSigRAW = isSigRAW(ip);
    for j = 1:sigCategoriesRAW
        fprintf(fid_RAW,'%s | %g\n',GOTableGeneric.GOName{isSigRAW(j)},surrogatePVals(isSigRAW(j),index));
    end
    fclose(fid_RAW);
end

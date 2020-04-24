function SpatialScoringCategories(whatSpecies,whatStructFilt)
% See if we can get more nuanced spatial scoring going
%-------------------------------------------------------------------------------
if nargin < 1
    whatSpecies = 'human';
end
if nargin < 2
    whatStructFilt = 'cortex';
end
params = GiveMeDefaultParams(whatSpecies,whatStructFilt);
%-------------------------------------------------------------------------------

% Load gene-expression data:
params.g.normalizationGene = 'zscore';
params.g.normalizationRegion = 'zscore';
[geneData,geneInfo,structInfo] = LoadMeG(params.g);
numGenes = height(geneInfo);

% Get GO categories:
GOTable = GiveMeGOData(params,geneInfo.entrez_id);
numGenes = height(geneInfo);
numGOCategories = height(GOTable);

% Get pairwise distances:
distMat = GiveMeDistanceMatrix(params);
getUpperDiag = @(x) x(triu(true(size(x)),+1));
distUpper = getUpperDiag(distMat);

%-------------------------------------------------------------------------------
% Look at spatial scale and strength of each GO category:
d0 = nan(numGOCategories,1);
A = nan(numGOCategories,1);
B = nan(numGOCategories,1);
R = nan(numGOCategories,1);
negRho = nan(numGOCategories,1);

parfor j = 1:numGOCategories
    hereIam = ismember(geneInfo.entrez_id,GOTable.annotations{j});
    geneData_j = geneData(:,hereIam);
    G = corr(geneData_j','rows','pairwise');
    gData = getUpperDiag(G);

    negRho(j) = -corr(distUpper,gData,...
                            'type','Spearman',...
                            'rows','pairwise');

    % Fit exponential to CGE:
    if strcmp(whatSpecies,'mouse')
        s = fitoptions('Method','NonlinearLeastSquares',...
                    'StartPoint',[1,0,0.7],'Lower',[0,-0.5,0],'Upper',[2,1,10]);
    else
        s = fitoptions('Method','NonlinearLeastSquares',...
                    'StartPoint',[1,0,0.01],'Lower',[0,-0.5,0],'Upper',[2,1,10]);
    end
    f = fittype('A*exp(-x*n) + B','options',s);
    try
        % [f_handle,Stats,c] = GiveMeFit(distUpper,getUpperDiag(G),'exp')
        [c,stats] = fit(distUpper,gData,f);
        % f_handle = @(x) c.A.*exp(-x*c.n) + c.B;
        % f = figure('color','w');
        % hold('on')
        % plot(distUpper,getUpperDiag(G),'.k')
        % dRange = linspace(0,max(distUpper),100);
        % plot(dRange,f_handle(dRange),'r','LineWidth',2)
        d0(j) = 1/c.n; % spatial scale of transcriptional autocorrelation in this category
        A(j) = c.A;
        B(j) = c.B;
        R(j) = stats.rsquare;
        % fprintf(1,'Category %u/%u succeeded\n',j,numGOCategories);
    catch
        % fprintf(1,'Category %u/%u failed\n',j,numGOCategories);
    end
end

%-------------------------------------------------------------------------------
% Annotate and save:
GOTable.A_fitted = A;
GOTable.B_fitted = B;
GOTable.d0_fitted = d0;
GOTable.R2fit = R;
GOTable.negRho = negRho
fileName = sprintf('CategorySpatialScoring_%s.mat',whatSpecies);
save(fullfile('DataOutputs',fileName),'GOTable');
fprintf(1,'Saved results to %s\n',fileName);

%-------------------------------------------------------------------------------
% Sort by R2 of the exponential fit:
[~,ix] = sort(GOTable.R2fit,'descend');
GOTable = GOTable(ix,:);

%-------------------------------------------------------------------------------
% Save it to csv file:
%-------------------------------------------------------------------------------
newTable = table(GOTable.GOName,GOTable.GOIDlabel,GOTable.GOID,...
                    GOTable.A_fitted,GOTable.B_fitted,GOTable.d0_fitted,...
                    GOTable.R2fit,GOTable.negRho);
fileOut = fullfile('SupplementaryTables',...
        sprintf('CategorySpatialScoring_%s-%s.csv',whatSpecies,whatStructFilt));
writetable(newTable,fileOut,'Delimiter',',','QuoteStrings',true);
fprintf(1,'Saved category spatial embedding scores to %s\n',fileOut);


% (now can load in other analysis functions to see whether these parameters
% are informative of FPSR characteristics)...
end

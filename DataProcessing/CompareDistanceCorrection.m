% Quick script to investigate

%-------------------------------------------------------------------------------
[gScore_1,geneEntrezIDs_1] = DoMe(false);
[gScore_2,geneEntrezIDs_2] = DoMe(true);


f = figure('color','w');
plot(gScore_1,gScore_2,'.k');
title(whatEdgeProperty)
xlabel('no distance correction')
ylabel('with distance correction')

%-------------------------------------------------------------------------------
function [gScore,geneEntrezIDs] = DoMe(correctDistance)
    % Just the casual 50 input arguments:
    whatEdgeProperty = 'distance';
    connectomeType = 'Oh-brain';
    absType = 'pos'; % 'pos','neg','abs' -> e.g., pos -> coexpression contribution increases with the statistic
    corrType = 'Spearman'; % {'Spearman','Pearson'};
    normalizationGene = 'none';
    normalizationRegion = 'none'; % {'none','zscore'}
    pThreshold = 0.05;
    energyOrDensity = 'energy';
    pValOrStat = 'stat';
    thresholdGoodGene = 0.5;


    % Compute edge-level statistics:
    [edgeData,regionStruct] = GiveMeEdgeStat(connectomeType,pThreshold,...
                            whatEdgeProperty,false);

    % Load in our gene data, properly processed:
    [geneData,geneInfo,structInfo] = LoadMeG({normalizationGene,...
                        normalizationRegion},energyOrDensity);

    % Check gene data matches connectome data
    if ~all([regionStruct.id]'==structInfo.id)
        % Take subset
        [~,ia,ib] = intersect([regionStruct.id]',structInfo.id,'stable');
        geneData = geneData(ib,:);
        structInfo = structInfo(ib,:);
        fprintf(1,'Gene data matched to subset of %u Allen regions\n',length(ib));
    end

    % Correct distance?
    if correctDistance
        C = load('Mouse_Connectivity_Data.mat','Dist_Matrix');
        distanceRegressor = C.Dist_Matrix{1,1};
        fprintf(1,'Regressing ipsilateral distances\n');
    else
        distanceRegressor = []; % just compute normal correlations
    end

    % Compute gene scores:
    % (sometimes entrez IDs change -- e.g., when matching to distance results)
    [gScore,geneEntrezIDs] = GiveMeGCC(edgeData,geneData,geneInfo.entrez_id,...
                                    corrType,distanceRegressor,...
                                    absType,thresholdGoodGene,pValOrStat);
end

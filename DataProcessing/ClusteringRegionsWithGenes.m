function ClusteringRegionsWithGenes(geneData,geneInfo,structInfo)

% ------------------------------------------------------------------------------
% Perform the clustering of regions using gene expression data:
% ------------------------------------------------------------------------------
% switch ClusterHow
% case 'linkage'
%     numGenesHere = size(dataMatrix,2);
%     R = BF_pdist(dataMatrix,distanceMetric);
%     links = linkage(R,linkageMethod);
%
%     % Compute cluster assignments at a given number of clusters:
%     Ci = cluster(links,'maxclust',numClusters);
%
% case 'kmeans'
%     numReplicates = 20;
%     % Now cluster at the given number of clusters:
%     numGenesHere = size(dataMatrix_noNaN,2);
%     fprintf(1,'k-means clustering to k = %u (%u replicates) using %u genes\n', ...
%                                     numClusters,numReplicates,numGenesHere);
%     [Ci,Centroids,sumd] = kmeans(dataMatrix_noNaN,numClusters,...
%                             'distance','sqEuclidean',...
%                             'replicates',numReplicates, ...
%                             'start','sample');
%
% otherwise
%     error('Unknown clustering method ''%s''',ClusterHow);
% end


% ------------------------------------------------------------------------------
% Plot:
% ------------------------------------------------------------------------------
figure('color','w'); box('on');
myColorMap = [flipud(BF_getcmap('blues',9,0)); BF_getcmap('reds',9,0)];
colormap(myColorMap)

clusterSizes = arrayfun(@(x)sum(Ci==x),1:numClusters);

% Loop over clusters
for i = 1:numClusters
    fCii = find(Ci==i);
    numMembers = length(fCii);

    % Strip of expression in this region:
    subplot(1,numClusters,i); box('on')
    % subplot(2,numClusters,numClusters+i); box('on')

    GenExp = dataMatrix(Ci==i,:);

    % Within each cluster, re-order regions using hierarchical clustering:
    ord_Ci = Cluster_Reorder(GenExp,distanceMetric,linkageMethod);

    % Reorder using clustered permutation of gene expression:
    ColPerm = G.GeneExpData_norm_cl.([energyOrDensity '_iy']);

    % Re-order the indices and the gene expression matrix:
    fCii = fCii(ord_Ci);
    GenExp = GenExp(ord_Ci,ColPerm(keepCol));

    pcolor([GenExp, zeros(size(GenExp,1),1); zeros(1,size(GenExp,2)+1)]);
    shading flat

    % Add regions:
    for j = 1:numMembers
        colorHere = rgbconv(G.RegionStruct(fCii(j)).color_hex_triplet);
        rectangle('Position',[-numGenesHere/8,j,numGenesHere/8+1,1], ...
                    'FaceColor',colorHere,'EdgeColor',colorHere)
    end
    % Add text:
    for j = 1:numMembers
        text(-numGenesHere/9,j+0.5,G.RegionStruct(fCii(j)).acronym,'FontSize',9)
    end

    % Rescale to have ~same row height for all regions:
    gcaPosition = get(gca,'Position');
    gcaPosition(4) = gcaPosition(4)*numMembers/max(clusterSizes);
    set(gca,'Position',gcaPosition);

    xlim([-numGenesHere/8,numGenesHere+1])
    ylabel('Regions')
    xlabel('Genes')
    set(gca,'YTick',[])

    title(sprintf('Module %u: Expression %s',i,energyOrDensity))

    % ------------------------------------------------------------------------------
    % Text output
    % ------------------------------------------------------------------------------
    fprintf(1,'\n------------------------------------------------------------------------------\n');
    fprintf(1,'Group %u\n',i);
    fprintf(1,'------------------------------------------------------------------------------\n');
    fCii = find(Ci==i);
    for j = 1:length(fCii)
        fprintf(1,'%s (%s)\n',G.RegionStruct(fCii(j)).acronym,...
                                G.RegionStruct(fCii(j)).MajorRegionName);
    end

end

end

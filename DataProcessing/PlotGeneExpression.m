function PlotGeneExpression(geneData,geneInfo,structInfo,doNormalize,structScores,geneScores)

if nargin < 4
    doNormalize = false;
end
if nargin < 5
    structScores = [];
end
if nargin < 6
    geneScores = [];
end

numStructs = height(structInfo);
numGenes = size(geneData,2);

if doNormalize
    geneData = BF_NormalizeMatrix(geneData,'scaledSigmoid');
end

if ~isempty(structScores)
    structScores = BF_NormalizeMatrix(structScores,'maxmin');
else
    structScores = ones(numStructs,1);
end

% ------------------------------------------------------------------------------
% Plot:
% ------------------------------------------------------------------------------
figure('color','w'); box('on');
myColorMap = [flipud(BF_getcmap('blues',9,0));1,1,1; BF_getcmap('reds',9,0)];
colormap(myColorMap)

if ~isempty(geneScores)
    subplot(5,1,1)
    plot(geneScores,'k')
    subplot(5,1,2:5)
end

% Re-order the indices and the gene expression matrix:
% pcolor([geneData, zeros(size(geneData,1),1); zeros(1,size(geneData,2)+1)]);
% shading flat
imagesc(geneData)

ratNum = 2.2;

% Add regions:
for j = 1:numStructs
    if ismember(structInfo.Properties.VariableNames,'color_hex_triplet');
        colorHere = rgbconv(structInfo.color_hex_triplet{j});
    else
        colorHere = 'k';
    end
    rectangle('Position',[-(numGenes/ratNum+1)*structScores(j),-0.5+j,(numGenes/ratNum+1)*structScores(j),1], ...
                'FaceColor',colorHere,'EdgeColor',colorHere)
end
% Add text:
for j = 1:numStructs
    text(-numGenes/(ratNum*1.2),j,structInfo.acronym(j),'FontSize',9)
end

xlim([-numGenes/ratNum,numGenes+1])
ylabel('Structures')
xlabel('Genes')
set(gca,'YTick',[])

end

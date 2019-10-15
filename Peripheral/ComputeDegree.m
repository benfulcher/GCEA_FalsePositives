function [k,structInfo] = ComputeDegree(params,doBinarize)
% Compute the degree across areas of a given brain parcellation
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams(params);
end
if nargin < 2
    doBinarize = true;
end

%-------------------------------------------------------------------------------
% Compute the binary degree for each brain region:
[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(params.c.connectomeSource,...
                                            params.c.pThreshold,...
                                            doBinarize,...
                                            params.c.whatWeightMeasure,...
                                            params.c.whatHemispheres,...
                                            params.c.structFilter);
k = sum(A_bin,1)' + sum(A_bin,2);

%-------------------------------------------------------------------------------
% Match to gene expression data:
[geneData,geneInfo,structInfo] = LoadMeG(params.g);

% Match structures:
switch params.humanOrMouse
case 'human'
    [~,ia,ib] = intersect(regionAcronyms,structInfo.ROI_ID,'stable');
case 'mouse'
    [~,ia,ib] = intersect(regionAcronyms,structInfo.acronym,'stable');
end
k = k(ia);
structInfo = structInfo(ib,:);

% [A_bin,geneData,structInfo,keepStruct] = filterStructures(params.c.structFilter,structInfo,A_bin,geneData);

end

function [k,structInfo] = ComputeDegree(params,doBinarize)

if nargin < 1
    params = GiveMeDefaultParams(params);
end
if nargin < 2
    doBinarize = true;
end
%-------------------------------------------------------------------------------

[A_bin,regionAcronyms,adjPVals] = GiveMeAdj(params.c.connectomeSource,...
                                            params.c.pThreshold,...
                                            doBinarize,...
                                            params.c.whatWeightMeasure,...
                                            params.c.whatHemispheres,...
                                            params.c.structFilter);

[geneData,geneInfo,structInfo] = LoadMeG(params.g);
[A_bin,geneData,structInfo,keepStruct] = filterStructures(params.c.structFilter,structInfo,A_bin,geneData);
k = sum(A_bin,1)' + sum(A_bin,2);

end

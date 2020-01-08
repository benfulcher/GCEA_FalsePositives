function GOTable = GiveMeGOData(params,entrezIDs)
% Wraps params vector into requirements for Matlab enrichment package,
% including the function `GetFilteredGOData'
%-------------------------------------------------------------------------------
if nargin < 1
    params = GiveMeDefaultParams();
end
if nargin < 2
    entrezIDs = [];
end

%-------------------------------------------------------------------------------
% Get GO information:
GOTable = GetFilteredGOData(params.e.dataSource,params.e.processFilter,...
                                params.e.sizeFilter,entrezIDs);

%-------------------------------------------------------------------------------
% Filter to a maximum number of annotations for category
% (for investigating null scaling for a fixed category size)
if isfield(params.e,'sizeFix') && ~isempty(params.e.sizeFix)
    warning('Taking a maximum of %u annotations per GO category :-O',params.e.sizeFix)
    GOTable = FixedSizeGO(GOTable,params.e.sizeFix);
end


end

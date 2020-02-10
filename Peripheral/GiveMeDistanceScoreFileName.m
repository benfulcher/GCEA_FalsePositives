function fileName = GiveMeDistanceScoreFileName(params)
% Get filename for saving/loading distance score data
%-------------------------------------------------------------------------------

fileName = sprintf('spatialACscores_%s-%s_%s.mat',params.humanOrMouse,...
            params.structFilter,params.gcc.pValOrStat);

end

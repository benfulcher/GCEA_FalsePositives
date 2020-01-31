function fileName = GiveMeDistanceScoreFileName(params)
% Get filename for saving/loading distance score data
%-------------------------------------------------------------------------------

fileName = sprintf('dScores_%s-%s_%s_%s_abs-%s_conn%u.mat',params.humanOrMouse,...
                        params.structFilter,params.gcc.whatCorr,params.gcc.pValOrStat,...
                        params.gcc.absType,params.gcc.onlyConnections);

end

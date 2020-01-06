function dataFileSurrogate = FindNullFile(whatSpecies,structFilter,numNullSamples)
% Setting for the ensemble enrichment parameters
%-------------------------------------------------------------------------------

if nargin < 1
    whatSpecies = 'mouse';
end
if nargin < 2
    structFilter = 'all';
end
if nargin < 3
    numNullSamples = 20000;
end

%-------------------------------------------------------------------------------
% In the case of some custom spatially constrained ensemble, can specify the file
% containing the info:
switch whatSpecies
case 'mouse'
    if strcmp(structFilter,'cortex')
        fprintf(1,'Spatial maps for mouse cortex\n');
        switch numNullSamples
        case 20000
            dataFileSurrogate = 'mouseCortexSurrogate_N20000_rho8_d040.csv';
        case 40000
            dataFileSurrogate = 'mouseCortexSurrogate_N40000_rho8_d0270.csv';
        end
    else
        fprintf(1,'Spatial maps for mouse whole brain\n');
        switch numNullSamples
        case 20000
            dataFileSurrogate = 'mouseSurrogate_N20000_rho8_d040.csv';
        case 40000
            dataFileSurrogate = 'mouseSurrogate_N40000_rho8_d078.csv';
        end
    end
case 'human'
    fprintf(1,'Spatial maps for human cortex\n');
    switch numNullSamples
    case 20000
        dataFileSurrogate = 'humanSurrogate_N20000_rho8_d02000.csv';
    case 40000
        dataFileSurrogate = 'humanSurrogate_N40000_rho8_d03500.csv';
    end
end

end

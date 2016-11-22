function RunErmineJ(inputFileName)
% Runs ermineJ on an input file:

% cf. http://erminej.chibi.ubc.ca/help/tutorials/erminej-cli/
%
% -a: sets the annotation file
% -b: sets bigger is better for gene scores
% -c: sets the gene set (class file): i.e., GO XML
% -d sets the data directory
% -e sets the column for the scores in your gene score file
% -g sets how to deal with scores for replicate genes
% -i number of iterations
% -j include gene symbols for all gene sets
% -M multiple test correction method
% -m method for computing raw class statistics for GSR
% -n method for computing gene set significance
% -o output filename
% -s gene score file
% -x maximum class size (100)
% -y minimum class size


% Set the ermineJ home directory:
% homeDir = '/home/benfulcher/Monash076/Ben/MouseConnectome/ermineJ/ermineJ-3.0.2';
% [status,cmdOut] = system('ERMINEJ_HOME=/Users/benfulcher/Downloads/ermineJ-3.0.2');

% Get settings for writing out to file:
shellScriptPath = '/Users/benfulcher/Downloads/ermineJ-3.0.2/bin/ermineJ.sh';
propertiesFilePath = which('ermineJBP.properties')
inputFilePath = which(inputFileName);
outputFile = fullfile(pwd,'tmp_ermineJ.out')

% Set the JAVA_HOME variable:
setenv('JAVA_HOME','/Library/Java/JavaVirtualMachines/jdk1.8.0_73.jdk/Contents/Home');

% Construct the command:
command = sprintf('%s -C %s -b -j -s %s -o %s',shellScriptPath,propertiesFilePath,inputFilePath,outputFile)

% Execute the command:
[status,cmdOut] = system(command);

% Read in the data:
[GOName,GOID,pval,corr_pval,numGenes,geneMembers] = ReadInErmineJ(outputFile);

% Write out to screen:
isSig = (corr_pval < 0.05);
fisSig = find(isSig);
numSig = length(fisSig);
for i = 1:numSig
    fprintf(1,'%s (%s) [p=%.2f]\n',GOName(fisSig(i)),GOID(fisSig(i)),corr_pval(fisSig(i)));
end

% $ERMINEJ_HOME/bin/ermineJ.sh -C $ERMINEJ_HOME/ermineJBP.properties -b -j -s
% $ERMINEJ_HOME/ermineJ.data/ermineJInputFile_GGblock_17349genes_k44_norm_energy_reciprocalConnected_distance_alldiv-expFitAll_tstat.txt
% -o ermineJ_reciprocalConnectedBP.out

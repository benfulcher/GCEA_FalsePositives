function [GSR10MconnectedEitherBP5,resultsTableRichFeed,resultsTableRichFeedEither] = ImportFulcher2016();

    %% Import data from text file.
    % Script for importing data from the following text file:
    %
    %    /Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Fulcher/GSR_10M_connectedEither_BP_5.csv
    %
    % To extend the code to different selected data or a different text file,
    % generate a function instead of a script.

    % Auto-generated by MATLAB on 2017/11/10 16:54:26

    %% Initialize variables.
    filename = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Fulcher/GSR_10M_connectedEither_BP_5.csv';
    delimiter = '\t';
    startRow = 2;

    %% Read columns of data as text:
    % For more information, see the TEXTSCAN documentation.
    formatSpec = '%s%s%[^\n\r]';

    %% Open the text file.
    fileID = fopen(filename,'r');

    %% Read columns of data according to the format.
    % This call is based on the structure of the file used to generate this
    % code. If an error occurs for a different file, try regenerating the code
    % from the Import Tool.
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

    %% Close the text file.
    fclose(fileID);

    %% Convert the contents of columns containing numeric text to numbers.
    % Replace non-numeric text with NaN.
    raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
    for col=1:length(dataArray)-1
        raw(1:length(dataArray{col}),col) = mat2cell(dataArray{col}, ones(length(dataArray{col}), 1));
    end
    numericData = NaN(size(dataArray{1},1),size(dataArray,2));

    for col=[1,2]
        % Converts text in the input cell array to numbers. Replaced non-numeric
        % text with NaN.
        rawData = dataArray{col};
        for row=1:size(rawData, 1)
            % Create a regular expression to detect and remove non-numeric prefixes and
            % suffixes.
            regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
            try
                result = regexp(rawData(row), regexstr, 'names');
                numbers = result.numbers;

                % Detected commas in non-thousand locations.
                invalidThousandsSeparator = false;
                if numbers.contains(',')
                    thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                    if isempty(regexp(numbers, thousandsRegExp, 'once'))
                        numbers = NaN;
                        invalidThousandsSeparator = true;
                    end
                end
                % Convert numeric text to numbers.
                if ~invalidThousandsSeparator
                    numbers = textscan(char(strrep(numbers, ',', '')), '%f');
                    numericData(row, col) = numbers{1};
                    raw{row, col} = numbers{1};
                end
            catch
                raw{row, col} = rawData{row};
            end
        end
    end


    %% Create output variable
    GSR10MconnectedEitherBP5 = table;
    GSR10MconnectedEitherBP5.GOID = cell2mat(raw(:, 1));
    GSR10MconnectedEitherBP5.pValCorr = cell2mat(raw(:, 2));

%===============================================================================
%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Fulcher/GSR_10M_richfeed_BP_10.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/11/10 16:58:21

%% Initialize variables.
filename = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Fulcher/GSR_10M_richfeed_BP_10.csv';
delimiter = '\t';
startRow = 2;

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
GSR10MrichfeedBP10 = table(dataArray{1:end-1}, 'VariableNames', {'GOCategory','Description','Numberofgenes','ermineJpvalue'});


%% Create output variable
resultsTableRichFeed = table();
GOtoNumber = @(x)str2num(x(4:end));
resultsTableRichFeed.GOID = cellfun(GOtoNumber,GSR10MrichfeedBP10.GOCategory);
resultsTableRichFeed.pValCorr = GSR10MrichfeedBP10.ermineJpvalue;


%===============================================================================
%===============================================================================
%===============================================================================
%% Import data from text file.
% Script for importing data from the following text file:
%
%    /Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Fulcher/GSR_10M_richfeedEither_BP_10.csv
%
% To extend the code to different selected data or a different text file,
% generate a function instead of a script.

% Auto-generated by MATLAB on 2017/11/10 17:04:56

%% Initialize variables.
filename = '/Users/benfulcher/GoogleDrive/Work/CurrentProjects/GeneExpressionEnrichment/DataSets/Fulcher/GSR_10M_richfeedEither_BP_10.csv';
delimiter = '\t';
startRow = 2;

%% Format for each line of text:
%   column1: text (%s)
%	column2: text (%s)
%   column3: double (%f)
%	column4: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%s%f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to the format.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'TextType', 'string', 'HeaderLines' ,startRow-1, 'ReturnOnError', false, 'EndOfLine', '\r\n');

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Create output variable
GSR10MrichfeedEitherBP10 = table(dataArray{1:end-1}, 'VariableNames', {'GOCategory','Description','Numberofgenes','ermineJpvalue'});

%% Create output variable
resultsTableRichFeedEither = table();
GOtoNumber = @(x)str2num(x(4:end));
resultsTableRichFeedEither.GOID = cellfun(GOtoNumber,GSR10MrichfeedEitherBP10.GOCategory);
resultsTableRichFeedEither.pValCorr = GSR10MrichfeedEitherBP10.ermineJpvalue;

end
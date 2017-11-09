function structInfo = GiveMeHCPNames()
% Loads in Simon's data on HCP names, and outputs a table
%
% % (see Table 1 of the following supplementary info for more information on parcels:
% https://images-nature-com.ezproxy.lib.monash.edu.au/full/nature-assets/nature/journal/v536/n7615/extref/nature18933-s3.pdf)
%-------------------------------------------------------------------------------

fid = fopen('ROInames_HCP.txt');
S = textscan(fid,'%u%s');
fclose(fid);
ID = S{1};
acronym = S{2};

numROIs = length(ID);

% (assuming all left and we know all are cortical):
isLeft = true(numROIs,1);
isCortex = true(numROIs,1);

% Generate a structure information table:
structInfo = table(ID,acronym,isLeft,isCortex);

end

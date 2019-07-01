function GOTable = ImportMeijer()

meijerResults = readtable('Meijer_etal_bioRxiv_2019.xlsx');

GOtoNumber = @(x)str2num(x(4:end));
GOID = cellfun(GOtoNumber,meijerResults.GOTerm);

pValCorr = meijerResults.FDRQ_value;

GOTable = table(GOID,pValCorr);

end

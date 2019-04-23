# MouseEdge
MouseEdge is a repository for obtaining enrichment signatures from spatial maps in human and mouse.

## Data processing

### Literature enrichment signatures
The script for importing enrichment data reported in existing literature is `ImportLiteratureEnrichment`.
It uses data in `/LiteratureEnrichmentData`.

First type of annotations are manually-curated, from studies that noted enrichment results in-text with no supplementary files for full results: `TableGOBPs.csv`.

The second type are using scripts (in `/DataProcessing/IndividualEnrichmentImportScripts/`) to directly process data provided as supplementary material from the following studies:
* `WhitakerReformatted.xlsx`: Whitaker, K. J. et al. Adolescence is associated with genomically patterned consolidation of the hubs of the human brain connectome. Proc. Natl. Acad. Sci. USA 113, 201601745–9110 (2016).
* `Vertes-rstb20150362supp1.xlsx`: Vértes, P. E. et al. Gene transcription profiles associated with inter-modular hubs and connection distance in human functional magnetic resonance imaging networks. Phil. Trans. Roy. Soc. B 371, 20150362 (2016).
* `Fulcher2016_connectedUnconnected_BP_TableS1.csv`, `Fulcher2016_richFeederPeripheral_BP_TableS5.csv`: Fulcher, B. D. & Fornito, A. A transcriptional signature of hub connectivity in the mouse connectome. Proc. Natl. Acad. Sci. USA 113, 1435–1440 (2016).
* `Tan2013-table-s6-david-200pos-transport.csv`: Tan, P. P. C., French, L. & Pavlidis, P. Neuron-Enriched Gene Expression Patterns are Regionally Anti-Correlated with Oligodendrocyte-Enriched Patterns in the Adult Mouse and Human Brain. Front. Psychiat. 7, (2013).
* `Parkes2017_PC1.txt`, `Parkes2017_PC2.txt`, `Parkes2017_PC5.txt`, `Parkes2017_PC9.txt`: 1.	Parkes, L., Fulcher, B. D., Yücel, M. & Fornito, A. Transcriptional signatures of connectomic subregions of the human striatum. Genes, Brain and Behavior 25, 1176–663 (2017).


## Analysis

### Enrichment signatures of spatially-correlated null maps

#### Generating surrogate maps
First generate pairwise distance matrices for the regions in human cortex and mouse brain:
```matlab
SaveOutDistanceMatrices
```
This generates `mouseDistMat.csv` and `humanDistMat.csv`.
These files can be used as input to python code by [Josh Burt et al.](https://github.com/benfulcher/surrogateMaps).
For example, here:
```
python3 GenerateMapsFixed.py
```

The outputs, `mouseSurrogate_N10000_rho8_d040.csv` and `mouseSurrogate_rho10.csv`, are spatially correlated null maps that can be visualized in Matlab. For example:
```matlab
VisualizeSurrogateMaps('mouse');
```

#### Conventional enrichment analysis on each null map

The following computes for 1000 of the null maps generated above:
```matlab
SurrogateEnrichment('mouse',1000,'spatialLag');
```
This generates results stored in the file `SurrogateGOTables_1000_mouse_spatialLag.mat` in `DataOutputs`.

#### Analyzing enrichment signatures of spatially correlated null maps

The results computed above can be read in and processed as a GO Table:
```matlab
GOTableNull = SurrogateEnrichmentProcess(whatSpecies,numMaps,whatSurrogate);
```

### Enrichment signatures of random-gene null maps

### Visualizing transcriptional data

Clustered row x gene expression matrices can be plotted for mouse:
```matlab
PlotExpressionMatrix('mouse')
```

And human:
```matlab
PlotExpressionMatrix('human')
```


### Non-specific spatial effects
These analyses look at quantifying nonspecific spatial patterning of gene-expression maps.

```matlab
DistanceConfoundResults()
```

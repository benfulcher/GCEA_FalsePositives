# MouseEdge
MouseEdge is a repository for obtaining enrichment signatures from spatial maps in human and mouse.

## Raw data import

The directory `/HumanData` requires three data files:
* `100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortex_ROI_NOdistCorrSurfaceANDEuclidean.mat`
* `100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortexSEPARATE_ROI_NOdistCorrSurfaceANDEuclidean.mat`
* `100DS360scaledRobustSigmoidNSGDSQC1Lcortex_ROI_NOdistCorrSurface.mat`

The directory `/MouseData` requires:
* `AllenGeneDataset_19419.mat`

## Enrichment data

This relies on a toolbox for matlab-based GO enrichment, which can be installed by cloning:
```
git clone git@github.com:benfulcher/GeneEnrichment.git
```
It is expected to be in the home directory, or otherwise accessible in the path.
This can be modified using `GiveMeFile('EnrichmentToolbox');`

It should be run to get the following data files:
* `GOTerms_BP.mat`
* `GOAnnotationDirect-mouse-biological_process-Prop.mat`
* `GOAnnotationDirect-human-biological_process-Prop.mat`

## Data processing

### Literature enrichment signatures

Information about enrichment results reported in published studies can be imported and processed by running
```matlab
ImportLiteratureEnrichment;
```
Results are saved as `LiteratureEnrichmentLoaded.mat`.
All data is read in from the `/LiteratureEnrichmentData`.

First type of annotations are manually-curated, from studies that noted enrichment results in-text with no supplementary files for full results: `TableGOBPs.csv`:

The second type are using scripts (in `/DataProcessing/IndividualEnrichmentImportScripts/`) to directly process data provided as supplementary material from the following studies:
* `WhitakerReformatted.xlsx`: Whitaker, K. J. et al. Adolescence is associated with genomically patterned consolidation of the hubs of the human brain connectome. Proc. Natl. Acad. Sci. USA 113, 201601745–9110 (2016).
* `Vertes-rstb20150362supp1.xlsx`: Vértes, P. E. et al. Gene transcription profiles associated with inter-modular hubs and connection distance in human functional magnetic resonance imaging networks. Phil. Trans. Roy. Soc. B 371, 20150362 (2016).
* `Fulcher2016_connectedUnconnected_BP_TableS1.csv`, `Fulcher2016_richFeederPeripheral_BP_TableS5.csv`: Fulcher, B. D. & Fornito, A. A transcriptional signature of hub connectivity in the mouse connectome. Proc. Natl. Acad. Sci. USA 113, 1435–1440 (2016).
* `Tan2013-table-s6-david-200pos-transport.csv`: Tan, P. P. C., French, L. & Pavlidis, P. Neuron-Enriched Gene Expression Patterns are Regionally Anti-Correlated with Oligodendrocyte-Enriched Patterns in the Adult Mouse and Human Brain. Front. Psychiat. 7, (2013).
* `Parkes2017_PC1.txt`, `Parkes2017_PC2.txt`, `Parkes2017_PC5.txt`, `Parkes2017_PC9.txt`: 1.	Parkes, L., Fulcher, B. D., Yücel, M. & Fornito, A. Transcriptional signatures of connectomic subregions of the human striatum. Genes, Brain and Behavior 25, 1176–663 (2017).


## Analysis

### Enrichment signatures of random spatial maps

This code computes the enrichment across 1000 actual independent random number samples:
```matlab
SurrogateEnrichment('mouse',1000,'randomUniform');
```

#### Analyzing enrichment signatures of null maps

The results computed above can be read in and processed as a GO Table for:

* Spatial lag model:
```matlab
GOTableNull = SurrogateEnrichmentProcess('mouse',1000,'spatialLag');
```

* Spatially random model:
```matlab
SurrogateEnrichment('mouse',5000,'randomUniform','');
```

* Spatially random model (plus coordinated shuffling of genes through space) [should be equivalent to previous]:
```matlab
SurrogateEnrichment(‘mouse’,10000,’randomUniform’,’coordinatedSpatialShuffle’);
```

* Spatially random model (plus independent shuffling of space, separately per gene) [should be no signal---a real null of correlated noise with noise]:

SurrogateEnrichment(‘mouse’,10000,’randomUniform’,’independentSpatialShuffle’);




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

This code computes the enrichment across 1000 null maps generated above:
```matlab
SurrogateEnrichment('mouse',1000,'spatialLag');
```
This generates results stored in the file `SurrogateGOTables_1000_mouse_spatialLag.mat` in `DataOutputs`.

#### Analyzing enrichment signatures of spatially correlated null maps

The results computed above can be read in and processed as a GO Table:
```matlab
GOTableNull = SurrogateEnrichmentProcess('mouse',1000,'spatialLag');
```




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

# EnrichmentNulls
EnrichmentNulls is a repository for obtaining enrichment signatures from spatial maps in human and mouse.

## Raw data import

The directory `/HumanData` requires three data files:
* `100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortex_ROI_NOdistCorrSurfaceANDEuclidean.mat`
* `100DS220scaledRobustSigmoidNSGDSQC1LcortexSubcortexSEPARATE_ROI_NOdistCorrSurfaceANDEuclidean.mat`
* `100DS360scaledRobustSigmoidNSGDSQC1Lcortex_ROI_NOdistCorrSurface.mat`

The directory `/MouseData` requires:
* `AllenGeneDataset_19419.mat`

## Enrichment data

This relies on a toolbox for Matlab-based GO enrichment, which can be installed by cloning:
```bash
git clone git@github.com:benfulcher/GeneEnrichment.git
```
It is expected to be in the home directory, or otherwise accessible in the path.
This can be modified using `GiveMeFile('EnrichmentToolbox');`

Please follow the instructions for that repository to compute (or download) the following required data files:
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
All data is read in from the `LiteratureEnrichmentData` directory.

First type of annotations are manually-curated, from studies that noted enrichment results in-text with no supplementary files for full results: `TableGOBPs.csv`:

The second type are using scripts (in `/DataProcessing/IndividualEnrichmentImportScripts/`) to directly process data provided as supplementary material from the following studies:
* `WhitakerReformatted.xlsx`: Whitaker, K. J. et al. Adolescence is associated with genomically patterned consolidation of the hubs of the human brain connectome. _Proc. Natl. Acad. Sci. USA+ **113**, 201601745–9110 (2016).
* `Vertes-rstb20150362supp1.xlsx`: Vértes, P. E. et al. Gene transcription profiles associated with inter-modular hubs and connection distance in human functional magnetic resonance imaging networks. _Phil. Trans. Roy. Soc. B_ **371**, 20150362 (2016).
* `Fulcher2016_connectedUnconnected_BP_TableS1.csv`, `Fulcher2016_richFeederPeripheral_BP_TableS5.csv`: Fulcher, B. D. & Fornito, A. A transcriptional signature of hub connectivity in the mouse connectome. _Proc. Natl. Acad. Sci. USA_ **113**, 1435–1440 (2016).
* `Tan2013-table-s6-david-200pos-transport.csv`: Tan, P. P. C., French, L. & Pavlidis, P. Neuron-Enriched Gene Expression Patterns are Regionally Anti-Correlated with Oligodendrocyte-Enriched Patterns in the Adult Mouse and Human Brain. Front. Psychiat. 7, (2013).
* `Parkes2017_PC1.txt`, `Parkes2017_PC2.txt`, `Parkes2017_PC5.txt`, `Parkes2017_PC9.txt`: 1.	Parkes, L., Fulcher, B. D., Yücel, M. & Fornito, A. Transcriptional signatures of connectomic subregions of the human striatum. _Genes, Brain and Behavior_ **25**, 1176–663 (2017).

## Precomputing

Note that batch job scripts for all bulk computations are in the `BatchComputing` directory.
There is a file for all mouse-related analyses: `batchAllMouseAnalyses.sh`, and for all human analyses: `batchAllHumanAnalyses.sh`.

### Intra-category coexpression

The within-category coexpression metric computation is also precomputed:

```matlab
IntraCorrelationByCategory('mouse','geneShuffle',20000,'VE1',true)
```

This saves the results as `Intra_mouse_geneShuffle_VE1_20000.mat`.

### Ensemble-based nulls

Null distributions for all GO categories is done using the companion Matlab package for both conventional gene-score enrichment and ensemble-based enrichment.

Here is an example, computing a null distribution for each GO category in mouse and human according to both `'randomMap'` and `'spatialLag'` null models (specified using the `'customEnsemble'` setting):

```matlab
species = {'mouse','human'};
numSpecies = length(species);
nullTypes = {'randomMap','customEnsemble'};
numNullTypes = length(nullTypes);
for i = 1:numSpecies
    % Load in all of the default parameters:
    params = GiveMeDefaultParams(species{i});
    for j = 1:numNullTypes
        params.e.whatEnsemble = nullTypes{j};

        % Wrapper for running ensemble-based nulls using appropriate gene-expression data:
        NullComputation(params);
    end
end
```

Null information is saved as a cell in a new column in the GO category table.
Each null distribution is stored in a table that is saved to a `.mat` file.
The name of this file is set in `params.e.fileNameOut`.

For example, this one is 40000 null samples for mouse data using the `'randomMap'` null ensemble:
* `PhenotypeNulls_40000_mouse_randomMap_Spearman_mean.mat`


## Analysis

### Enrichment signatures of random spatial maps

#### Precomputing
This code computes the enrichment across 1000 actual independent random number samples:
```matlab
SurrogateEnrichment('mouse',1000,'randomUniform');
```

This pre-computation is required for the analyses presented here (cf. in `batchAllHumanAnalyses` and `batchAllMouseAnalyses`).
The empty second input uses the default number of maps.

```matlab
% Spatially random model (plus independent shuffling of space, separately per gene) [should be no signal---a real null of correlated noise with noise]:
SurrogateEnrichment('mouse',[],'randomUniform','independentSpatialShuffle');
% Spatially random model:
SurrogateEnrichment('mouse',[],'randomUniform','');
% Spatial lag model:
SurrogateEnrichment('mouse',[],'spatialLag','');
% (and similarly for 'human')
```

There is also the (irrelevant) spatially random model (plus coordinated shuffling of genes through space) [equivalent to not doing a coordinatedSpatialShuffle]:
```matlab
SurrogateEnrichment('mouse',[],'randomUniform','coordinatedSpatialShuffle');
```

The computed results are saved a `.mat` files and can be read in and processed as a GO Table using `SurrogateEnrichmentProcess`.

#### Analysis
Plotting distributions of FPSE across GO categories for the three null cases:
```matlab
NullEnrichmentTogether('mouse',[],true)
NullEnrichmentTogether('human',[],true)
```

Outputting a table, and understanding some statistics of FPSE in mouse and human:
```matlab
FPSRTable();
```

Investigating the overlap between literature annotations and FPSE as histograms:
```matlab
OverlapLitFPSR('mouse')
OverlapLitFPSR('human')
```

### Specific GO categories

You can zoom into specific GO categories using:

```
PlotCategoryNullCompare
```

#### Generating surrogate maps

##### Matlab

```matlab
numMaps = 40000;
plotSummary = true;
GenerateSpatialEnsemble('mouse','all',plotSummary,numMaps)
GenerateSpatialEnsemble('mouse','cortex',plotSummary,numMaps)
GenerateSpatialEnsemble('human','cortex',plotSummary,numMaps)
```

##### Python (old)
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

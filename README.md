# MouseEdge
MouseEdge is a repository for obtaining enrichment signatures from spatial maps in human and mouse.

## Data processing
### Literature enrichment signatures
Data is in `/LiteratureEnrichmentData`.

First type of annotations are manually-curated, from studies that noted enrichment results in-text with no supplementary files for full results: `TableGOBPs.csv`.

The second type are using scripts (in `/DataProcessing/IndividualEnrichmentImportScripts/`) to directly process data provided as supplementary material from the following studies:
* `WhitakerReformatted.xlsx`: Whitaker, K. J. et al. Adolescence is associated with genomically patterned consolidation of the hubs of the human brain connectome. Proc. Natl. Acad. Sci. USA 113, 201601745–9110 (2016).
* `Vertes-rstb20150362supp1.xlsx`: Vértes, P. E. et al. Gene transcription profiles associated with inter-modular hubs and connection distance in human functional magnetic resonance imaging networks. Phil. Trans. Roy. Soc. B 371, 20150362 (2016).
* `Fulcher`: Fulcher, B. D. & Fornito, A. A transcriptional signature of hub connectivity in the mouse connectome. Proc. Natl. Acad. Sci. USA 113, 1435–1440 (2016).
* 


## Analysis
### Non-specific spatial effects
These analyses look at quantifying nonspecific spatial patterning of gene expression maps.

```matlab
DistanceConfoundResults()
```

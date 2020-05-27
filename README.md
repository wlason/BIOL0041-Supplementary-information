# BIOL0041-Supplementary-information

This GitHub repository provides supplementary information for BIOL0041 project **APOBEC-linked dormancy and exhaustion in cancer**
## Main directory
* _PLOS-Biology-template.pdf_ - example manuscript showing editing required for PLOS Biology
## Research paper
### Contains supplementary figures and files for the research paper
#### R Markdown scripts

* _PHATE and Fisher analysis.Rmd_ - script used to produce correlation plots in Figure 1 and conduct all gene enrichment analyses
* _PHATE cancer and control.Rmd_ - script used to produce clustering heatmaps in Figure 2
* _Hypoxia analysis.Rmd_ - script used to produce panels in Figure 4
* _Immune infiltration analysis.Rmd_ - script used to produce panels in Figure 4
* _sigs_Myriad_mutect2.R_ - script used to calculate mutational signatures for Figure 1 and Figure 4 using high-performance computing

#### Supplementary files

* __S01__ - table summarising evidence for dormancy and exhaustion gene signatures used in the analysis
* __S02-S13__ - Fisher's exact test results for COSMIC mutations in k-means clusters created using various genesets and following analysis of pathway annotation and gene ontology
* __S14__ - Fisher's exact test results for hypoxia and all programmes clustering factoring in the predicted effect of mutation (only deleterious mutations were used)
* __S15__ - Bindea immune infiltration estimates as PHATE maps and corresponding quantification boxplots

## Laboratory report
### Contains supplementary figures and files for the laboratory report
#### R Markdown scripts

* **_BIOL0041-R-script.Rmd_ - laboratory notebook with all analyses performed over the course of the project**
* _MAF_Myriad.R_ - scrip used to download single nucleotide variation data using high-performance computing
* _sigs_Myriad.R_ - script used to calculate mutational signatures using high-performance computing
* _ConsensusClusterPlus.R_ - script used to calculate optimal cluster assignments using high-performance computing
* _FisherTest.R_ and _FisherTestAdvanced.R_ - function to perform enrichment analysis on clusters identified by ConsensusClusterPlus using high-performance computing

#### Supplementary files

* __S01-S18__ - Fisher's exact test results for COSMIC mutations in clusters created identified using hierarchical clustering, consensus hierarchical clustering, and PHATE together with the following analysis of pathway annotation and gene ontology
* __S19__ - Fisher's exact test results for hypoxia and all programmes clustering factoring in the predicted effect of mutation (only deleterious mutations were used)
* __S20__ - raw mutational signature contributions calculated using various pipelines and COSMIC v3 mutational signature definitions
* __S21__ - benchmarking plots for ConsensusClusterPlus run on k in 1-10
* __S22__ - PHATE maps for every programme in (APOBEC, dormancy, hypoxia) clustered on dormancy, exhaustion, and all programmes
* __S23__ - by-tissue correlation plots for APOBEC, dormancy and exhaustion programmes and ConsensusTME cell scores
* __S24__ - Bindea immune infiltration estimates as PHATE maps and corresponding quantification boxplots

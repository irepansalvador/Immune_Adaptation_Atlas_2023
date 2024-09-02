# Immune Adaptation Atlas 2023
This repository contains the notebooks related to the manuscript "Quantifying adaptive evolution of the human immune cell landscape" by Salvador-Martínez, Murga-Moreno, et al.

Below, we specify the contents of each subfoder.

## Developmental

We provide Jupyter notebooks to obtain, for each developmental compartment, the list of Differentially expressed genes for each cell type.
The annotated data for each compartment was downloaded in the AnnData format from https://developmental.cellatlas.io/fetal-immune.
We used the Python scanpy package to read and process the h5ad file.

## Adult

We provide Jupyter notebooks to obtain, for each adult compartment, the list of Differentially expressed genes for each cell type.
The annotated data for each compartment was downloaded in the AnnData format from https://www.tissueimmunecellatlas.org/.
We used the Python scanpy package to read and process the h5ad file.

## Lung 

We provide Jupyter notebooks to obtain, for each adult compartment, the list of Differentially expressed genes for each cell type.
The annotated data for each compartment was downloaded in the AnnData format from https://5locationslung.cellgeni.sanger.ac.uk/.
We used the Python scanpy package to read and process the h5ad file.
We also provide the script to produce the heatmap with the spatial loadings based on the deconvolution of the Spatial transcriptomics data.

## Macrophages

Macrophage activation DE genes were obtained from the preprint by Panousis et al (2023), available in https://doi.org/10.1101/2023.05.29.542425 as Supplementary material.
We provide a R markdown file where we further process these DE genes. Specifically, from these DE gene lists we only considered genes with a log2fc > 1 and with a p-value < 0.05.
For each stimulus, we further categorised the genes as “early” (found exclusively in DE gene list after 6h), “late” (found exclusively in DE gene list after 24h) or “sustained” (found at both time points).

## Phylostratigraphy

We obtained the Phylostrata (PS) values for human genes from the study of Litman and Stein (2019; https://doi.org/10.1053/j.seminoncol.2018.11.002), where the authors determine consensus ages for human protein-coding and noncoding genes, using publicly available databases. Based on the data provided by Litman and Stein, we generated the file "Human_Gene_Ages.csv" that contains PS values for 27,974 human genes. We provide Rmd files to run the PS enrichment test using a Hypergeometric test. The test was performed using the `phyper()` function of the R stats package.

## GO_terms_2024

We provide the scripts to perform GO term analyses using the Biological Process (BP) ontology for each set of DE genes used for the ABC-MK test, using the R packages cluster-Profiler (version 4.2.2), enrichplot (version 1.14.2), GO.db and org.Hs.eg.db (BioConductor version 3.14).

## abcmk

We provide the scripts to run progressive bootstrapping, ABC-MK inferences and FDR.

## annotations

Folder containing polymorphic data from 1000GP and human fixations from Human-Chimp-Gorilla-Orangutan orthologs alignemnts, Hg38 gene coordinates, mammals orthologs list and confounding factors to run progressive resampling.

# pTcells: Cell-specific regulatory feature exctraction from public scRNA-seq data

### Table of contents:
* Biological framework
* Pipeline overview
* Installation and execution
* Tutorial
* Bibliography

### Biological framework: 
pTcells uses Differential Expression Analysis (DEA) tools and the Ensembl database to exctract regulatory elements from those genes that are overexpressed in T cells against . By contrasting T cells to another cell type of interest 

ContexTF uses the manually curated transcription factor database from model organisms to perform a homology-based search for putative transcription factors in unannotated prokaryotic genomes. It combines the results of three different homology searches: sequence, structure and genomic context. It is widely assumed that proteins with high sequence and, especially, structural similarity have high probability of performing the same function (Krissinel, 2007). For this reason, most bioinformatic tools take advantage of these characteristic to annotate homologous proteins throughout species. In prokaryotic genomes, where most TFs are found in context with the proteins that they regulate, the conservation of the genes near the regulatory protein could also be used to annotate proteins in prokaryotic organisms. ContexTF integrates a classical search of homologous candidates and provides a snapshot of the genomic contexts of the involved proteins to evaluate the conservation and similarity of the flanking genes of the candidates.

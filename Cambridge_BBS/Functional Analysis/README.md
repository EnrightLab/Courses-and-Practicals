# Functional Analysis of Differential Expression Data

Differential Gene expression analysis usually yields lists of genes from most up-regulated to most down-regulated.

In this practical we will briefly explore using such tools to search for biological commonalities in our datasets.

We will work with the genelist obtained from the mRNA Seq data in the last practical.

You should download the genelist file from here to your Desktop:
mir_210_protein_coding_list.txt


## Gene Ontology Analysis

The mission of the GO Consortium is to develop a comprehensive, computational model of biological systems, ranging from the molecular to the organism level, across the multiplicity of species in the tree of life. The Gene Ontology (GO) knowledgebase is the world’s largest source of information on the functions of genes. This knowledge is both human-readable and machine-readable, and is a foundation for computational analysis of large-scale molecular biology and genetics experiments in biomedical research.

The ontology is split into three categories for describing genes

* Cellular Compopnent
* Biological Process
* Moleclar Function

Let's load our genelist into Cytoscape and explore Gene Ontology Enrichments using the BinGO Package.

File -> Import -> Network -> File -> miR_210_protein_coding_list.txt
![import](import.png,"Logo Title Text 1")

Make sure you select the log FC column as node data and the first column as the source Node



## GeneMania


# Functional Analysis of Differential Expression Data

Differential Gene expression analysis usually yields lists of genes from most up-regulated to most down-regulated.

In this practical we will briefly explore using such tools to search for biological commonalities in our datasets.

We will work with the genelist obtained from the mRNA Seq data in the last practical.

You should download the genelist file from here to your Desktop:
mir_210_protein_coding_list.txt

## Sylamer Analysis of Genelists for microRNA signatures

The purpose of our miR-210 overexpression experiment was to see if an over-expressed microRNA would influence gene-expression in MCF7 cells. To try and understand the physiological role of miR-210 regulation in this system.

We will use [Sylamer](http://wwwdev.ebi.ac.uk/enright-dev/sylarray2/) to take the gene Names from the genelists we produced. These are sorted according to log2 Fold change. This algorithm will search for 6nt motifs that correlate with expression changes. We are hoping to see a motif for the seed match to miR-210.

![sylamer](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional%20Analysis/sylamer.png)

* Choose upload gene list
** Then pick mir_210_short_list.txt
* Then select *homo sapiens* as the *species*
* Then choose *external gene name* as the *Input Gene Ids* 
* Then click on *Launch*

## Gene Ontology Analysis

The mission of the [GO Consortium](http://geneontology.org) is to develop a comprehensive, computational model of biological systems, ranging from the molecular to the organism level, across the multiplicity of species in the tree of life. The Gene Ontology (GO) knowledgebase is the world’s largest source of information on the functions of genes. This knowledge is both human-readable and machine-readable, and is a foundation for computational analysis of large-scale molecular biology and genetics experiments in biomedical research.

The ontology is split into three categories for describing genes

* Cellular Compopnent
* Biological Process
* Moleclar Function

Let's load our genelist into Cytoscape and explore Gene Ontology Enrichments using the BinGO Package.

File -> Import -> Network -> File -> miR_210_protein_coding_list.txt
![import](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional%20Analysis/import.png)

Make sure you select the log FC column as node data and the first column as the source Node



## GeneMania


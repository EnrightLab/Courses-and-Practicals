# Functional Analysis of Differential Expression Data

**Anton Enright (aje39@cam.ac.uk)  Stephanie Wenlock (scw65@cam.ac.uk) **

## Overview 

Differential Gene expression analysis usually yields lists of genes from most up-regulated to most down-regulated.

In this practical we will briefly explore using such tools to search for biological commonalities in our datasets.

We will work with the genelist obtained from the mRNA Seq data in the last practical.

You should download the genelist file from here to your Desktop:
[mir_210_short_list.txt](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/mir_210_short_list.txt)

## Sylamer Analysis of Genelists for microRNA signatures

The purpose of our miR-210 overexpression experiment from last week was to see if an over-expressed microRNA would influence gene-expression in MCF7 cells. To try and understand the physiological role of miR-210 regulation in this system.

We will use [Sylamer](http://wwwdev.ebi.ac.uk/enright-dev/sylarray2/) to take the gene Names from the genelists we produced. These are sorted according to log2 Fold change. This algorithm will search for 6nt motifs that correlate with expression changes. We are hoping to see a motif for the seed match to miR-210.

Information about Sylamer is in this short presentation:
[PDF](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/bbs_short_talk.pdf)


* Launch [Sylamer](http://wwwdev.ebi.ac.uk/enright-dev/sylarray2/)
* Choose *Start Analysis Now*
* Choose upload gene list
  * Then pick mir_210_short_list.txt
* Then select *homo sapiens (hsa)* as the *species*
* Then choose *HGNC_Symbol* as the *Input Gene Ids* 
* Then click on *Launch*

This will produce a Sylamer analysis of 6nt motifs from matched 3' UTR sequences from our sorted genelist from the miR-210 experiment.

Reference: **Stijn van Dongen, Cei Abreu-Goodger & Anton J. Enright; Nature Methods(2008)**

## MicroRNA Seed level analysis

The miR-210 microRNA has the following sequence:
```
AGCCCCAGCCCACCGCACACUG
```

its reverse complementary seed match would be expected to be: ```CTGGGGC```

A binding site would hence look something like:
```
AGCCCCAGCCCACCGCACACUG
 |||||||
NCGGGGTCNNNNNNNNNNNNNN
```

We should see enrichments of ```CTGGGGC``` motifs if the microRNA has caused a significant effect on mRNA expression.

## Gene Ontology Analysis

We can also explore the functional annotations of the same genelist. We can arbitrarily select up-regulated and down-regulated genes to test for functional commonalities using the Gene ontology (GO).

The mission of the [GO Consortium](http://geneontology.org) is to develop a comprehensive, computational model of biological systems, ranging from the molecular to the organism level, across the multiplicity of species in the tree of life. The Gene Ontology (GO) knowledgebase is the world’s largest source of information on the functions of genes. This knowledge is both human-readable and machine-readable, and is a foundation for computational analysis of large-scale molecular biology and genetics experiments in biomedical research.

The ontology is split into three categories for describing genes

* Cellular Compopnent
* Biological Process
* Moleclar Function

Let's load our genelist into Cytoscape and explore Gene Ontology Enrichments using the BinGO Package.

Download the full file here:
[mir_210_protein_coding_list.txt](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/mir_210_protein_coding_list.txt)

*File -> Import -> Network -> File -> miR_210_protein_coding_list.txt*

![import](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/import.png)

Make sure you select the log FC column as node data and the first column as the source Node, i.e. **Green Document Icon**.


The genes will load as an unconnected network. Make sure the gene list is sorted on log2Fold change.

We will now select genes from the list. Highlight the top 50 or so up or down-regulated genes.
1. Right click then choose *Select nodes from selected rows*
2. Then choose *Apps* -> *BinGO*
![cytoscape](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/genelist.png)
3. Enter *upregulated* or *downregulated* as the **Cluster Name**
4. Choose *Homo Sapiens* as the species and *GO Full* as the ontology file.
5. Click *Start BinGO*
6. You can make a hiearchical layout of the results by choosing:
  *Layout* -> *Hierarchical Layout*

![bingo1](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/bingo.png)

![bingo2](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/bingo_result.png)

## GeneMania

Gene Mania is a probabalistic network of gene functions according to large-scale genomics datasets. It is an interesting alternative to human curated datasets such as *GO*. 

1. Go to the [GeneMania](https://genemania.org) website. 
2. Paste in gene names from either up-regulated or down-regulated genes.
3. Download and open in a spreadsheet the [mir_210_protein_coding_list.txt](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/mir_210_protein_coding_list.txt).
4. Select about 100 either up or down-regulated genes from either end of the list and paste them into the search box.
5. It will automaticaly parse and process the genes and build a functional network.

![genemania](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/GeneMania.png)

## Reactome

[Reactome](http://www.reactome.org) is a slightly different view of function that is based on human curated biological pathways.

1. On Reactome choose *Analyze Data*
2. Click on *Choose File* and the hit *Continue*
3. Project to human is fine, clicl *Analyze!*

![reactome](https://github.com/EnrightLab/Courses-and-Practicals/blob/master/Cambridge_BBS/Functional_Analysis/reactome.png)

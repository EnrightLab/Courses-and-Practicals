<<<<<<< HEAD

# WTAC Advanced Course - Small RNA Analysis

## A conserved RNA regulates miRNA turnover and animal behavior through a near-perfect miRNA site

Angelo Bitetti<sup>1,2,3,\*</sup>, Allison C Mallory<sup>1,2,3,\*</sup>, Claudia Carrieri<sup>4</sup>, 
Elisabetta Golini<sup>5</sup>, Hector Carreño Gutierrez<sup>6</sup>, Emerald Perlas<sup>7</sup>, Yuvia A. Pérez-Rico<sup>1,2,3,8</sup>, 
Glauco P. Tocchini-Valentini<sup>5</sup>, Anton J. Enright<sup>9</sup>, William H. J. Norton<sup>6</sup>, 
Silvia Mandillo<sup>5</sup>, Dónal O'Carroll<sup>4</sup> and Alena Shkumatava<sup>1,2,3</sup>	

*These authors contributed equally to this work.

<sub>1. Institut Curie, 26 Rue d’Ulm, 75248 Paris Cedex 05, France.</sub><br>
<sub>2. CNRS UMR3215, 75248 Paris Cedex 05, France.</sub><br>
<sub>3. INSERM U934, 75248 Paris Cedex 05, France..</sub><br>
<sub>4. MRC Centre for Regenerative Medicine, Institute for Stem Cell Research, School of Biological Sciences, University of Edinburgh, 5 Little France Drive, Edinburgh, UK.</sub><br>
<sub>5. Consiglio Nazionale delle Ricerche, Istituto di Biologia Cellulare e Neurobiologia, European Mouse Mutant Archive-Infrafrontier-International Mouse Phenotyping Consortium, I-00015 Monterotondo Scalo (Rome), Italy.</sub><br>
<sub>6. Department of Neuroscience, Psychology and Behaviour, University of Leicester, Leicester, UK.</sub><br>
<sub>7. European Molecular Biology Laboratory (EMBL), Mouse Biology Unit, Monterotondo Scalo, Italy.</sub><br>
<sub>8. INSERM, U900, 75248 Paris Cedex 05, France.</sub><br>
<sub>9. European Molecular Biology Laboratory (EMBL), European Bioinformatics Institute, Wellcome Trust Genome Campus, Hinxton, Cambridge, CB10 1SD, United Kingdom.</sub><br>

## Code used for small RNA analyses
The codebase describes the small RNA analysis and was authored by Anton Enright (aje@ebi.ac.uk) March-May 2017.

* [The R/BioConductor code and results (markdown) for small RNA analysis](alena_new_data_mar_2017.md)
* The [Rmarkdown](alena_new_data_mar_2017.rmd) for small RNA analysis
* The [final table](mouse_results.txt) (tab-delimited) of normalised miRNA counts and statistical results for each miR

### Initial Data
* [The raw count-table derived from the chimiRa mapping tool](mouse_counts_mar_2017.txt)
* [The experimental design (pdata) file](pdata_mar_2017.txt)
* [An accessory table containing small RNA lengths extracted from the chimiRa mapping tool](length_tables_mouse_mar_2017.txt)

## Version Information

[ChimiRa](http://www.ebi.ac.uk/research/enright/software/chimira) was used to trim adapters and map reads to microRNA precursors<br>
R version used was: _R 3.3.1_ (2016-06-21) -- "Bug in Your Hair"<br>
BioConductor version used was: _Bioconductor version 3.4_<br><br>

Dependencies: *RColorBrewer, gplots, DESeq2, reshape2, ggplot2*<br>
=======
Part II - BBS - Bioinformatics 
===============================
![Cambridge](/images/cambridge.jpg)

2019 Part II - BBS - Bioinformatics - Differential Gene Expression
-------------------------------------------------------------------

**Instructors: Anton Enright, Jack Monahan

|Anton Enright (aje39@cam.ac.uk) | Jack Monahan (monahanj@ebi.ac.uk) | |
|---------------------------|------------------------------------|------------------------------------|

***

<img src="../images/cambridge.jpg" align="right" width="150">

_Department of Pathology,  
Tennis Court Road,  
Cambridge, UK._  

***

Presentation:
------------
1. [Anton's Combined Presentations from the course](Anton_Presentations_WTAC_2018.pdf)

List of Practicals:
------------------

1. [Introduction to R/BioConductor](Intro_R/Intro_R_Practical.md)

2. Basics of SmallRNA Sequencing Data Analysis
   * [SmallRNA Sequencing - Basic Analysis and QC](small_RNA_seq/Practical_1/Practical_1.md)
   * [SmallRNA Sequencing - Mapping and Statistical analysis of count data](small_RNA_seq/Practical_2/Practical_2.md)

3. Worked Example for microRNA Target Analysis - MiR-210 in Breast Cancer
   * [SmallRNA Sequencing - Breast Cancer miR210 Samples - Participant Driven](miR_210_Experiment/small_RNASeq/small_RNASeq.md)
   * [Mini mapping practical using HISAT2](miR_210_Experiment/mini_mapping/)
   * [mRNA Sequencing - Breast Cancer miR210 Samples](miR_210_Experiment/mRNA_Seq/mRNA_Seq.md)
   
4. Downstream Analysis of microRNA data
   * Target prediction from Differential Gene Expression Data
   * **De novo** prediction of microRNAs from Small RNA Seq Data
   
5. Downstream Analysis of Differential Gene Expression Data
  * Ontology Analysis
  * Pathway Analysis
  * Functional Enrichment
  * Gene Set Enrichment Analysis
***

Useful resources:

* [Enright Lab](http://www.ebi.ac.uk/research/enright)
* [R Cran](https://cran.r-project.org/)
* [R Studio](http://www.rstudio.com/)
* [BioConductor](http://www.bioconductor.org)
* [DESeq2 Manual](http://bioconductor.org/packages/release/bioc/html/DESeq2.html)
* [ChimiRa Webserver](http://wwwdev.ebi.ac.uk/enright-dev/chimira/)
* [miRnovo Webserver](http://wwwdev.ebi.ac.uk/enright-dev/mirnovo/)
* [miRBase](http://www.mirbase.org)
* [Legacy Course Material, pre 2017](http://wwwdev.ebi.ac.uk/enright-srv/courses)
>>>>>>> 06843df9683f23ca871f9a97a4182743443b7ab0


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
* [The raw count-table derived from a full bowtie mappign](mouse_counts_mar_2017.txt)
* [The experimental design (pdata) file](pdata_mar_2017.txt)
* [An accessory table containing small RNA lengths extracted from the chimiRa mapping tool](length_tables_mouse_mar_2017.txt)

## Version Information

Bowtie2 was used (version 2.5.1) 
miRBase (release 22.1)
R version used was: _R 4.1.1 (2021-08-10) (2016-06-21)_  -- "Kick Things"<br>
BioConductor version used was: _Bioconductor version 3.4_<br><br>

Dependencies: *RColorBrewer, gplots, DESeq2, reshape2, ggplot2, rtsne*<br>

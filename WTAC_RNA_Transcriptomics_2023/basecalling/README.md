Wellcome Trust Advanced Courses - RNA Transcriptomics
===============================
![WTAC](../../images/acsc_logo.png)

![Wellcome](../../images/wellcome_logo.png)
![Cambridge](../../images/cambridge.jpg)

Basecalling Nanopore Data with Guppy
-------------------------------------------------------------------

**Instructors: Anton Enright, Stephanie Wenlock**

|Anton Enright (aje39@cam.ac.uk) | Stephanie Wenlock (scw65@cam.ac.uk)|
|---------------------------|------------------------------------|
|<img src="../../images/anton.jpg" height="150">|<img src="../../images/steph.png" height="150">|

***

<img src="../../images/cambridge.jpg" align="right" width="150">

_Department of Pathology,  
Tennis Court Road,  
Cambridge, UK._  

***

Guppy is a basecaller that turns raw signal data from a MinION or PromethION or similar nanopore device into nucleotide sequences. Guppy is a complex deep-learning AI software that tries to evaluate the most likely *k-mer* of 5nt that generated a particular signal. Because of this it is extremely slow.

# GPU Accelleration

Because running guppy on a CPU is extremely low, many people find it easier to use the very fast instruction set and processing power (matrix compute) of a high-end graphics card. To do this you need a computer with a proper installation of a fast card with a CUDA instruction set version >6.1.

Compatible GPUs include:

* GeForce RTX 4090	8.9
* GeForce RTX 4080	8.9
* GeForce RTX 4070 Ti	8.9
* GeForce RTX 3090 Ti	8.6
* GeForce RTX 3090	8.6
* GeForce RTX 3080 Ti	8.6
* GeForce RTX 3080	8.6
* GeForce RTX 3070 Ti	8.6
* GeForce RTX 3070	8.6
* Geforce RTX 3060 Ti	8.6
* Geforce RTX 3060	8.6
* GeForce GTX 1650 Ti	7.5
* TITAN RTX	7.5
* Geforce RTX 2080 Ti	7.5
* Geforce RTX 2080	7.5
* Geforce RTX 2070	7.5
* Geforce RTX 2060	7.5
* TITAN V	7.0
* TITAN Xp	6.1
* TITAN X	6.1
* GeForce GTX 1080 Ti	6.1
* GeForce GTX 1080	6.1
* GeForce GTX 1070 Ti	6.1
* GeForce GTX 1070	6.1
* GeForce GTX 1060	6.1
* GeForce GTX 1050	6.1

Running guppy on a compatible gpu will make it 1000s of times faster.

# Obtaining Guppy 

Guppy is downloaded from the Oxford Nanopore Community Site (https://community.nanoporetech.com)
We will be using the CPU version of guppy today.

# Running Guppy on Direct RNA Sequencing Data.

First we want to open up a terminal on our computer and navigate to where our FAST5 output folder is stored. We will be working with Group1s data in this section.

``` r
cd /home/manager/Desktop/course_data/basecalling/fast5/group1
ls  
```

This should now display the folder where the minION results are stored:
```
20230625_1833_MN38030_FAX28374_68a4d9b2
```

# Folder contents
lets ask to list all files in that folder and its subfolders:
``` r
ls -Ral 20230625_1833_MN38030_FAX28374_68a4d9b2

```


```
total 17120
drwxrwxrwx  13 aje  staff      416 27 Jun 09:02 .
drwxrwxrwx   4 aje  staff      128 28 Jun 07:45 ..
-rwxrwxrwx   1 aje  staff      233 26 Jun 12:53 barcode_alignment_FAX28374_68a4d9b2_0c72d8dd.tsv
drwxrwxrwx  10 aje  staff      320 27 Jun 09:02 fast5
-rwxrwxrwx   1 aje  staff      618 26 Jun 12:53 final_summary_FAX28374_68a4d9b2_0c72d8dd.txt
drwxrwxrwx   3 aje  staff       96 27 Jun 09:02 other_reports
-rwxrwxrwx   1 aje  staff   412975 26 Jun 12:53 pore_activity_FAX28374_68a4d9b2_0c72d8dd.csv
-rwxrwxrwx   1 aje  staff  1131588 26 Jun 12:53 report_FAX28374_20230625_1836_68a4d9b2.html
-rwxrwxrwx   1 aje  staff   326937 26 Jun 12:53 report_FAX28374_20230625_1836_68a4d9b2.json
-rwxrwxrwx   1 aje  staff   467401 26 Jun 12:53 report_FAX28374_20230625_1836_68a4d9b2.md
-rwxrwxrwx   1 aje  staff      169 26 Jun 12:53 sample_sheet_FAX28374_20230625_1836_68a4d9b2.csv
-rwxrwxrwx   1 aje  staff  6352484 26 Jun 12:53 sequencing_summary_FAX28374_68a4d9b2_0c72d8dd.txt
-rwxrwxrwx   1 aje  staff    53004 26 Jun 12:53 throughput_FAX28374_68a4d9b2_0c72d8dd.csv

20230625_1833_MN38030_FAX28374_68a4d9b2/fast5:
total 3831984
drwxrwxrwx  10 aje  staff        320 27 Jun 09:02 .
drwxrwxrwx  13 aje  staff        416 27 Jun 09:02 ..
-rwxrwxrwx   1 aje  staff  180067860 25 Jun 19:55 FAX28374_68a4d9b2_0c72d8dd_0.fast5
-rwxrwxrwx   1 aje  staff  205274607 25 Jun 21:05 FAX28374_68a4d9b2_0c72d8dd_1.fast5
-rwxrwxrwx   1 aje  staff  283476842 25 Jun 22:20 FAX28374_68a4d9b2_0c72d8dd_2.fast5
-rwxrwxrwx   1 aje  staff  213899906 25 Jun 23:47 FAX28374_68a4d9b2_0c72d8dd_3.fast5
-rwxrwxrwx   1 aje  staff  311989482 26 Jun 01:41 FAX28374_68a4d9b2_0c72d8dd_4.fast5
-rwxrwxrwx   1 aje  staff  238818340 26 Jun 04:59 FAX28374_68a4d9b2_0c72d8dd_5.fast5
-rwxrwxrwx   1 aje  staff  435969939 26 Jun 12:06 FAX28374_68a4d9b2_0c72d8dd_6.fast5
-rwxrwxrwx   1 aje  staff   92466960 26 Jun 12:53 FAX28374_68a4d9b2_0c72d8dd_7.fast5

20230625_1833_MN38030_FAX28374_68a4d9b2/other_reports:
total 9504
drwxrwxrwx   3 aje  staff       96 27 Jun 09:02 .
drwxrwxrwx  13 aje  staff      416 27 Jun 09:02 ..
-rwxrwxrwx   1 aje  staff  4864891 26 Jun 11:49 pore_scan_data_FAX28374_68a4d9b2_0c72d8dd.csv
```

You should see 8 FAST5 files that contain the compressed raw signal data. Each of these files contains 4000 reads. A big nanopore experiment may produce many thousands of these files. This will be the input for guppy.

# Running Guppy

``` r
../../ont-guppy-cpu/bin/guppy_basecaller -i 20230625_1833_MN38030_FAX28374_68a4d9b2/fast5 -s basecalled --kit SQK-RNA002 --flowcell FLO-MIN106

```

This will detect the FAST5 files and apply basecalling according to our kit **SQK-RNA002** and the flowcell type **FLO-MIN106**.

Finished sequences will be placed in the *basecalled* folder that we specified.


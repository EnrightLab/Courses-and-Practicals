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
|<img src="../../images/anton.jpg" height="150">|<img src="../images/steph.png" height="150">|

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



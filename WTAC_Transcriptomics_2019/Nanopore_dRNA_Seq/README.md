# Nanopore Direct RNA sequencing analysis

Jack Monahan & Anton Enright, EMBL-EBI

27th June 2019



## Why and where to use dRNA-Seq

#### Cons

- Lower yield than cDNA seq => 500Mb to 1.5Gb per flowcell 
- Slighly higher error ~ 11%
- Requires a massive quantity of polyA+ RNA (or of a specific target)



#### Pros

- Library preparation much simpler, fewer steps => less biaised

- Longer reads than cDNA

  ![](pictures/Slide17.png) 

- Better exons connectivity

![](pictures/exon_align.png)

* PolyAs can be detected and measured (hopefully)

  ![](pictures/Slide09.png)

* (Some) RNA modifications can modify the signal

  ![](pictures/Slide34.png)

   



## ONT multiFast5/Fast5 file format

MINKnow generates files containing the raw intensity signal in [HDF5 format](https://support.hdfgroup.org/HDF5/). The latest Fast5 format contains multiple reads per file.

![](pictures/HDF5.jpeg)



Files can be explored using **HDFview**

**Multifast5 containing raw data only**

![](pictures/fast5_pre.png)




**Example of raw signal**

![](pictures/fast5_raw.png)

![](pictures/Raw1.png)

![](pictures/Raw2.png)

![](pictures/Raw3.png)



## Useful tools for dRNA-Seq ONT analysis

### Basecalling

* ~~Albacore~~ (available to ONT community, development ceased)
* Guppy (available to ONT community) 
* [Scrappie](https://github.com/nanoporetech/scrappie)
* [Flappie](https://github.com/nanoporetech/flappie)  = Flip-flop basecaller, for DNA and cDNA
* [Chiron](https://github.com/haotianteng/chiron) = Community alternative

See basecaller comparison => https://github.com/rrwick/Basecalling-comparison

### Quality control

* [pycoQC](https://github.com/a-slide/pycoQC) = ONT data QC from sequencing summary file
* [NanoPack](https://github.com/wdecoster/nanopack) = Suite of tools to QC and process raw ONT data

### RNA alignment

* [**Minimap2** ](https://github.com/lh3/minimap2)
* [LAST](http://last.cbrc.jp)
* [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
* [STAR](https://github.com/alexdobin/STAR)

### Read Polishing

* [Nanopolish event align](https://nanopolish.readthedocs.io/en/latest/)
* [Tombo resquiggle](https://nanoporetech.github.io/tombo/)

### DNA/RNA modification detection
* Nanopolish
* [Tombo detect_modifications](https://nanoporetech.github.io/tombo/)
* [Nanopolish call-methylation](https://nanopolish.readthedocs.io/en/latest/quickstart_call_methylation.html)


# Mini-Practical

1. Get you data !

   ```bash
   cd ~/Desktop/course_data/nanopore_dRNA_Seq/datasets/
   ```

   ```
   tar xvf ${Sample}.tar.gz
   ```

2. Inspect reads with the HDFView GUI

   ```
   cd ~/Desktop/course_data/nanopore_dRNA_Seq/datasets/
   ```

   ```bash
   hdfview
   ```

     

3. Basecall your data with Albacore

   ```bash
   cd ~/Desktop/course_data/nanopore_dRNA_Seq/
  
   ```

   ```bash
   ./ont-guppy-cpu/bin/guppy_basecaller --help
   
   ./ont-guppy-cpu/bin/guppy_basecaller --print_workflows
   
   ./ont-guppy-cpu/bin/guppy_basecaller -i wt1 -s wt1_basecalls --flowcell FLO-MIN106 --kit SQK-RNA002 -q 0 --enable_trimming true --trim_strategy rna --reverse_sequence true --pt_scaling --qscore_filtering 0
   
   ./ont-guppy-cpu/bin/guppy_basecaller -i wt1 -s scr1_basecalls --flowcell FLO-MIN106 --kit SQK-RNA002 -q 0 --enable_trimming true --trim_strategy rna --reverse_sequence true --pt_scaling --qscore_filtering 0
   
   ```

   *Flowcell and Kit information can be found in the fast5 files*

   *With your sample data it should take around 10 mins*
   


4. QC the basecalled files with pycoQC

   https://github.com/a-slide/pycoQC

   ```bash
   pycoQC -f guppy/${Sample}/sequencing_summary.txt -o ${Sample}.pycoQC.html
   ```
   

5. Align reads against the transcriptome or the genome with Minimap2

   https://github.com/lh3/minimap2

   Merge reads
   ```bash
   cat guppy/${Sample}/pass/*.fastq > ${Sample}.fastq
   ```

   *Spliced alignment against genome*
   ```bash
   minimap2 -ax splice -uf -k 14 -L -t 8 ../references/Mus_musculus_genome.fa.gz ${Sample}.fastq | samtools view -bh -F 2308 | samtools sort -o reads.bam
   ```

    *Unspliced alignment against transcriptome*

   ```bash
   minimap2 -ax map-ont -L -t 8 ../references/Mus_musculus_transcriptome.fa.gz ${Sample}.fastq | samtools view -bh -F 2308 | samtools sort -o transcriptome.bam
   ```

   

6. Visualise genome-aligned reads with IGV

   https://software.broadinstitute.org/software/igv/download

   *Index reads first for visualization*

   ```bash
   samtools index reads.bam
   ```

7. Generate transcript counts with Salmon

   ```bash

   salmon quant --noErrorModel -p 4 -t ../references/Mus_musculus_transcriptome.fa.gz -l U -a transcriptome.bam  -o salmon/$Sample
   ```

8. Inspect transcript counts

   ```bash

   less salmon/$Sample/quant.sf
   ```



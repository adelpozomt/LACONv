# LACONv



LACONv has been developed to characterize copy number variations (i.e. deletion and duplication) in samples sequenced with short-read targeted kits.

CNVs are identified on exons within the region of interest using their normalized read depth. LACONv does not require a reference sample set, as each sample plays the role of a control for the remaining set.

For each exon it is calculated a dosage value and its corresponding z-value with respect to the dosage distribution of the 'reference' set. Cutoffs are determined under the assumption that the samples carry very rare CNVs and the whole sample set is not homogeneous, although it is possible to include relatives within the set of samples to be tested. However, the inclusion of samples with suspected recurrent CNVs is not recommended.

The chromosomal sex of each sample must be determined to divide the samples into two groups in order to characterize CNVs on X and Y chromosomes. If there is no possibility of perform sexing, CNVs on chrX and chrY are not calculated.

Also, if the size of sample set is too scarce or sample read depth is very poor, CNVs can not be characterized. Samples with low read depth are not considered as well as exons whose depth along the sample set is under some threshold.

The name of the samples and exons that are filtered out are registered in a log file in order to get track of the process.



## Input data



LACONv needs a list of indexed BAM/CRAM files that can be provided in a text file (one path per line) or by command line. 

Input configuration file must be included. Check the format in [this](https://github.com/adelpozomt/LACONv/blob/17cf3414767c3b04a9b9966737088a8a096a4fef/cfg/example.cfg) example file. 

Also, it must be provided a bed file with the genomic intervals that want to be tested. BED file must be formatted before running the LACONv script in this way:



```python
python preprocess_analysis_bed.py --b <input bed file> --gtf <input annotated gft> --o <output path>
```



## Output files



LACONv returns for each sample the following files:



* *.all_intervals.bed, that includes all the exons or intervals of the input bed file with their dosage and z-scores values
* *.cnv.bed, that corresponds to the intervals that have been characterized as DEL or DUP
* *_CNV.xlsx that it is an excel table with many sheets that explains in a complete way the CNV analysis results.



Additionally,  the tool returns a bed file with the intervals that have not been considered for the analysis and a .cfg file with a summary of the parameters configured for the CNV detection.



## Running example



```python
python CNV_analysis.py --f <bam1 file> --f <bam2 file> --o <output path> --cfg <config.cfg> --b <analysis bed> --clean
```



## Dependencies



python 2.7. Modules:

- numpy
- scipy
- pybedtools
- pandas



External tools:

- mosdepth >= 0.2.6
- samtools >= 1.3.1
- GenomeAnalysisTK v3.3-0



## 

# Illumina read trimming and QC

## Requirements

### Data needed
The commands below assume that paired-end Illumina reads are in separate *.fastq files that have been named using the convention "SampleName_R1" and "SampleName_R2".

### Software needed (versions indicated are those used in the 2022 paper)
* `trimmomatic` [v0.39] (http://www.usadellab.org/cms/index.php?page=trimmomatic)
* `fastQC` [v0.11.7] (https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)

_Dependening on how these programs are installed, it may be necessary to install other binaries or programs as well. Please reference the indicated documentation for details._  
</br>

## Procedure
*The commands below assume that the necessary software has been invoked at each step; the commands used to invoke each software program will depend on the mode of installation and individual system setup and so are not included here.* 

### Step 1: Use Trimmomatic to remove Illumina adapters, trim the 5' and 3' ends using a quality cutoff of 3, further trim each read using a sliding window of 4 bp and a quality cutoff of 25, and then filter out any reads less than 40 bp in length
*The filtered Illumina reads are used for assembly polishing and are thus saved in a new sample-specific working subdirectory named "/Polishing/Illumina_trim/". This directory must be created prior to running the commands below.*

*Also, the appropriate Illumina Adapter File needs to be saved in an accessible directory and the path name modified accordingly.*
```
cd /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim

trimmomatic PE -phred33 -summary StrainName_trimmomatic_summary.txt /Path/RawIlluminaReads/SampleName_R1.fastq /Path/RawIlluminaReads/SampleName_R2.fastq -baseout StrainName_trimmed.fastq ILLUMINACLIP:/Path/IlluminaAdatperFile/NexteraPE-PE.fa:2:30:10:6:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:25 MINLEN:40
```
The Illumina Nextera Adapter file used in the 2021 paper can be found [here](NexteraPE-PE.fa).

### Step 2: Use FastQC to assess quality of trimmed Illumina reads
```
cd /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/
mkdir FastQC

fastqc -o fastqc/ ./StrainName_trimmed_1P.fastq
fastqc -o fastqc/ ./StrainName_trimmed_2P.fastq
```

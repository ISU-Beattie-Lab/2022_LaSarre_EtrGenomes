# PacBio read trimming and QC

## Requirements

### Data needed
The commands below assume that the PacBio reads are in legacy *.bax.h5 format and that reads for a single sample (e.g., strain) have been demultiplexed into a sample-specific directory.

### Software needed (versions indicated are those used in the 2022 paper)
* smrtlink [v7.0.1] (https://www.pacb.com/wp-content/uploads/SMRT_Link_Installation_v701.pdf)
* Nanoplot [v1.38.1] (https://github.com/wdecoster/NanoPlot)

_Depending on how these programs are installed, it may be necessary to install other binaries or programs as well. Please reference the indicated documentation for details._  
</br>

## Procedure
*The commands below assume that the necessary software has been invoked at each step; the commands used to invoke each software program will depend on the mode of installation and individual system setup and so are not included here.* 

### Step 1: Use smrtlink to filter PacBio reads for quality and length
*The filtered PacBio reads are used for assembly polishing and are thus saved in a new sample-specific working subdirectory named "/Polishing/PBreads/". This directory must be created prior to running the commands below.*
```
cd /Path/WorkingDirectory/StrainName/Polishing/PBreads

for j in p0.1 p0.2 p0.3; do bax2bam -o SampleName_PacBio_${j} /Path/SampleReadDirectory/*_${j}.bax.h5 --subread  --pulsefeatures=DeletionQV,DeletionTag,InsertionQV,IPD,MergeQV,SubstitutionQV,PulseWidth,SubstitutionTag; done

for j in p0.1 p0.2 p0.3; do dataset create --type SubreadSet --name SampleName_${j} SampleName_${j}_subreadset.xml *_${j}.subreads.bam; done

for j in p0.1 p0.2 p0.3; do dataset filter SampleName_${j}_subreadset.xml SampleName_${j}_subreadset_filtered.xml 'rq>80' 'length>1000'; done

dataset merge SampleName_all_subreadset_filtered.xml SampleName_p0.1_subreadset_filtered.xml SampleName_p0.2_subreadset_filtered.xml SampleName_p0.3_subreadset_filtered.xml

dataset consolidate SampleName_all_subreadset_filtered.xml SampleName_subreads_filtered.bam SampleName_subreads_filtered.xml

bam2fastq -o SampleName_subreads_filtered SampleName_subreads_filtered.bam
```

### Step 2: Use NanoPlot to assess quality of filtered PacBio reads
```
cd /Path/WorkingDirectory/StrainName/Polishing/PBreads/
mkdir ReadQC

NanoPlot --fastq ./SampleName_subreads_filtered.fastq.gz -o ./ReadQC/SampleName_subreads_filtered_nanoplot --title SampleName_subreads_filtered --loglength
```

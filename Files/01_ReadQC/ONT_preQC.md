# Nanopore read trimming and QC

## Requirements

### Data needed
* Demultiplexed raw ONT reads in *.fastq.gz format

### Software needed (versions indicated are those used in the 2022 paper)
* Porechop [v0.2.4] (https://github.com/rrwick/Porechop)
* Nanofilt [v2.8.0] (https://github.com/wdecoster/nanofilt)
* Nanoplot [v1.38.1] (https://github.com/wdecoster/NanoPlot)


_Depending on how these programs are installed, it may be necessary to install other binaries or programs as well. Please reference the indicated documentation for details._  
</br>

## Procedure
*The commands below assume that the necessary software has been invoked at each step; the commands used to invoke each software program will depend on the mode of installation and individual system setup and so are not included here.* 

### Step 1: Concatenate all ONT reads for a given strain into a single file to facilitate processing and assembly and save in strain-specific working directory

```
zcat /Path/SampleReadDirectory/FAQ*.fastq.gz | gzip -c > /Path/WorkingDirectory/StrainName/StrainName.fastq.gz
```

### Step 2: Working within the strain-specific working directory, use Porechop to remove ONT adapters
```
cd /Path/WorkingDirectory/StrainName/

porechop -i StrainName.fastq.gz > StrainName_ONTtrim_pass_all.fastq
gzip StrainName_ONTtrim_pass_all.fastq
```

### Step 3: Use NanoFilt to remove reads shorter than 1000 bp
```
gunzip -c StrainName_ONTtrim_pass_all.fastq | NanoFilt -l 1000 | gzip > ./StrainName_ONTtrimfilt.fastq.gz
```

### Step 4: Assess the resulting read quality using NanoPlot
```
NanoPlot --fastq StrainName_ONTtrimfilt.fastq.gz -o /Path/WorkingDirectory/StrainName/read_QC/StrainName_ONTtrimfilt_nanoplot --title StrainName_ONTtrimfilt --loglength
```
</br>

# Polishing of consensus assembly using Nanopore and PacBio reads

## Requirements

Before running this step, you should have generated a consensus assembly of your genome using Trycycler, as described [here](../02_Assembly/Trycycler_reconciliation.md).

### Data needed
* Trimmed and filtered ONT reads in *.fastq.gz format
* Filtered PacBio reads in *.fastq.gz format

*If you have Illumina reads rather than PacBio reads, please use the polishing pipeline described [here](Polishing_ONTplusIllumina.md).*

### Software needed (indicated versions are those used for the 2022 paper)
* `Medaka` [v1.4.1] (https://github.com/nanoporetech/medaka)
* `Racon` [v1.6.1] (https://github.com/isovic/racon)
* `smrtlink` [v7.0.1] (https://www.pacb.com/wp-content/uploads/SMRT_Link_Installation_v701.pdf)
* `Circlator` [v1.5.5] (https://sanger-pathogens.github.io/circlator/)

_Dependening on how these programs are installed, it may be necessary to install other binaries or programs as well. Please reference the indicated documentation for details._ 
</br>

## Procedure
*The commands below assume that the necessary software has been invoked at each step; the commands used to invoke each software program will depend on the mode of installation and individual system setup and so are not included here.*  
</br>
### Step 1: Polish the consensus assembly with ONT reads using Racon
```
cd /Path/WorkingDirectory/StrainName/Polishing/

bwa index Path/WorkingDirectory/StrainName/trycycler/StrainName_ONTtrimfilt_5assemb_consensus.fasta 

minimap2 -a -x map-ont /Path/WorkingDirectory/StrainName/trycycler/StrainName_ONTtrimfilt_5assemb_consensus.fasta /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz | samtools sort -O SAM -o StrainName_ONTtrimfilt_5assemb_consensus.minialign.sorted.sam

racon -m 8 -x 6 -g -8 -w 500 -u /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz StrainName_ONTtrimfilt_5assemb_consensus.minialign.sorted.sam /Path/WorkingDirectory/StrainName/trycycler/StrainName_ONTtrimfilt_5assemb_consensus.fasta > StrainName_ONTtrimfilt_5assemb_racon1.fasta
```
</br>

### Step 2:  Repeat Racon polishing a second time
```
cd /Path/WorkingDirectory/StrainName/Polishing/

bwa index /Path/WorkingDirectory/StrainName/Polishing/StrainName_ONTtrimfilt_5assemb_racon1.fasta 

minimap2 -a -x map-ont /Path/WorkingDirectory/StrainName/Poloishing/StrainName_ONTtrimfilt_5assemb_racon1.fasta /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz | samtools sort -O SAM -o StrainName_ONTtrimfilt_5assemb_racon1.minialign.sorted.sam

racon -m 8 -x 6 -g -8 -w 500 -u /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz StrainName_ONTtrimfilt_5assemb_racon1.minialign.sorted.sam /Path/WorkingDirectory/StrainName/Polishing/StrainName_ONTtrimfilt_5assemb_racon1.fasta > StrainName_ONTtrimfilt_5assemb_racon2.fasta
```
</br>

### Step 3: Polish the Racon-polished assembly once with ONT reads using Medaka
*The specified basecaller model (indicated with "-m") matches the ONT reads used in the 2021 paper but may need to be modified for other ONT read sets.*
```
cd /Path/WorkingDirectory/StrainName/Polishing/

medaka_consensus -i /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz -d StrainName_ONTtrimfilt_5assemb_racon2.fasta -o medaka -m r941_min_hac_g507 -t 12
```
</br>

### Step 4: Use Arrow to polish the ONT-polished assembly using PacBio reads (following read mapping using pbmm2; both programs are part of SMRT Link)
```
cd /Path/WorkingDirectory/StrainName/Polishing/
mkdir Arrow && cd Arrow

pbmm2 index ../StrainName_ONTtrimfilt_5assemb_racon2medaka.fasta StrainName_ONTtrimfilt_5assemb_racon2medaka.mmi

pbmm2 align StrainName_ONTtrimfilt_5assemb_racon2medaka.mmi ../PBreads/StrainName_subreads_filtered.bam | samtools sort > StrainName_ONTtrimfilt_5assemb_racon2medaka_pbmm2align.bam

pbindex StrainName_ONTtrimfilt_5assemb_racon2medaka_pbmm2align.bam

variantCaller -j8 --algorithm=arrow --log-file ./StrainName_ONTtrimfilt_5assemb_arrow1_log.txt -r ../StrainName_ONTtrimfilt_5assemb_racon2medaka.fasta -o ./StrainName_ONTtrimfilt_5assemb_arrow1.gff -o ../StrainName_ONTtrimfilt_5assemb_arrow1.fasta StrainName_ONTtrimfilt_5assemb_racon2medaka_pbmm2align.bam
```  
</br>

### Step 5: Use Circulator to adjust contig start sites
```
cd /Path/WorkingDirectory/StrainName/Polishing/
mkdir circlator && cd circlator

circlator fixstart ../StrainName_ONTtrimfilt_5assemb_arrow1.fasta ./StrainName_ONTtrimfilt_5assemb_arrow1_fix
	
cp ./StrainName_ONTtrimfilt_5assemb_arrow1_fix.fasta ../StrainName_ONTtrimfilt_5assemb_arrow1_fix.fasta
```   
</br>

### Step 6: Use Arrow to do a final polishing of the rotated contigs using PacBio reads
```
cd /Path/WorkingDirectory/StrainName/Polishing/

samtools faidx StrainName_ONTtrimfilt_5assemb_arrow1_fix.fasta

cd Arrow

pbmm2 index ../StrainName_ONTtrimfilt_5assemb_arrow1_fix.fasta StrainName_ONTtrimfilt_5assemb_arrow1_fix.mmi

pbmm2 align StrainName_ONTtrimfilt_5assemb_arrow1_fix.mmi ../PBreads/StrainName_subreads_filtered.bam | samtools sort > StrainName_ONTtrimfilt_5assemb_arrow1_fix_pbmm2align.bam

pbindex StrainName_ONTtrimfilt_5assemb_arrow1_fix_pbmm2align.bam

variantCaller -j8 --algorithm=arrow --log-file ./StrainName_ONTtrimfilt_5assemb_arrow2_log.txt -r ../StrainName_ONTtrimfilt_5assemb_arrow1_fix.fasta -o ./StrainName_ONTtrimfilt_5assemb_Polished.gff -o ../StrainName_ONTtrimfilt_5assemb_Polished.fasta StrainName_ONTtrimfilt_5assemb_arrow1_fix_pbmm2align.bam
```  
</br>

## The FINAL polished assembly will be located at:
`/Path/WorkingDirectory/StrainName/Polishing/StrainName_ONTtrimfilt_5assemb_Polished.fasta`
# Polishing of consensus assembly using Nanopore and PacBio reads

## Requirements

Before running this step, you should have generated a consensus assembly of your genome using Trycycler, as described [here](../02_Assembly/Trycycler_reconciliation.md).

### Data needed
* Trimmed and filtered ONT reads in *.fastq.gz format
* Trimmed and filtered Illumina reads in *.fastq.gz format

*If you have PacBio reads rather than Illumina reads, please use the polishing pipeline described [here](Polishing_ONTplusPacBio.md).*

### Software needed (indicated versions are those used for the 2022 paper)
* `Medaka` [v1.4.1] (https://github.com/nanoporetech/medaka)
* `Racon` [v1.6.1] (https://github.com/isovic/racon)
* `MaSuRCA` [v4.0.5] (https://github.com/alekseyzimin/masurca)
* `NextPolish` [v1.4.0] (https://github.com/Nextomics/NextPolish)
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

### Step 4: Polish the ONT-polished assembly with Illumina reads once using POLCA (part of MaSuRCA)
```
cd /Path/WorkingDirectory/StrainName/Polishing/
mkdir POLCANextPolish && cd POLCANextPolish

polca.sh -a /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt_5assemb_racon2medaka.fasta -r '/Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_1P.fastq /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_2P.fastq /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_1U.fastq /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_2U.fastq'
```

</br>

### Step 5: Polish the POLCA-polished assembly with Illumina reads once using NextPolish
```
cd /Path/WorkingDirectory/StrainName/Polishing/POLCANextPolish

ls /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_1P.fastq /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_2P.fastq /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_1U.fastq /Path/WorkingDirectory/StrainName/Polishing/Illumina_trim/StrainName_trimmed_2U.fastq > sgs.fofn

nextPolish NextPolish1.run.conf

cp /Path/WorkingDirectory/StrainName/Polishing/POLCANextPolish/Round1/genome.nextpolish.fasta /Path/WorkingDirectory/StrainName/Polishing/StrainName_ONTtrimfilt_5assemb_POLCANext.fasta
```

The `NextPolish1.run.conf` configuration file can be found [here](NextPolish1.run.conf).  
</br>


### Step 6: Use Circulator to adjust contig start sites
```
cd /Path/WorkingDirectory/StrainName/Polishing/
mkdir Circlator && cd Circlator

circlator fixstart ../StrainName_ONTtrimfilt_5assemb_POLCANext.fasta ./StrainName_ONTtrimfilt_5assemb_POLCANext_fix

cp ./StrainName_ONTtrimfilt_5assemb_POLCANext_fix.fasta ../StrainName_ONTtrimfilt_5assemb_POLCANext_fix.fasta
``` 
</br>


### Step 7: Use NextPolish to do a final polishing of the rotated contigs using Illumina reads
```
cd /Path/WorkingDirectory/StrainName/Polishing/POLCANextPolish

nextPolish NextPolish2.run.conf

cp /Path/WorkingDirectory/StrainName/Polishing/POLCANextPolish/Round2/genome.nextpolish.fasta /Path/WorkingDirectory/StrainName/Polishing/StrainName_ONTtrimfilt_5assemb_Polished.fasta
```
The `NextPolish2.run.conf` configuration file can be found [here](NextPolish2.run.conf).


</br>

## The FINAL polished assembly will be located at:
`/Path/WorkingDirectory/StrainName/Polishing/StrainName_ONTtrimfilt_5assemb_Polished.fasta`
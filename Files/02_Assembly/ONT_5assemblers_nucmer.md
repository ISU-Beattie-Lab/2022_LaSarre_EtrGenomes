# ONT read assembly using five assemblers and identification of over-circularized contigs

## Requirements

### Data needed
* Trimmed and filtered ONT reads saved as a single .fastq.gz file in a strain-specific working directory

### Software needed (indicated versions are those used in the 2022 paper)
* NextDenovo [v2.5.0] (https://github.com/Nextomics/NextDenovo)
* NextPolish [v1.4.0] (https://github.com/Nextomics/NextPolish)
* Flye [v2.9] (https://github.com/fenderglass/Flye)
* Raven [v1.6.1] (https://github.com/lbcb-sci/raven)
* miniasm [v0.3] (https://github.com/lh3/miniasm)
* minipolish [v0.1.3] (https://github.com/rrwick/Minipolish)
* Unicycler (https://github.com/rrwick/Unicycler)
* MUMmer [v3.23] (https://mummer.sourceforge.net/)



_Dependening on how these programs are installed, it may be necessary to install other binaries or programs as well. Please reference the indicated documentation for details._  
</br>

## Procedure
*The commands below assume that the necessary software has been invoked at each step; the commands used to invoke each software program will depend on the mode of installation and individual system setup and so are not included here.*  
</br>

### Step 1: Assemble the trimmed and filtered ONT reads using five different assemblers
*Within each strain-specific directory, a new assembler-specific subdirectory is created to store the output of each assembler.*   
</br>

### Assembler 1: NextDenovo/NextPolish
*The two configuration files need to be saved within the NextDenovo subdirectory once created.*
```
cd /Path/WorkingDirectory/StrainName/
mkdir NextDenovo && cd NextDenovo

ls /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz > input.fofn
nextDenovo ND.run.cfg

mkdir NextPolish && cd NextPolish
ls /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz > lgs.fofn
nextPolish NP.run.cfg
```

*The ND.run.cfg configuration file for NextDenovo can be found [here](ND.run.cfg) (should be modified before use).*  
*The NP.run.cfg configuration file for NextPolish can be found [here](NP.run.cfg) (should be modified before use).*

</br>

### Assembler 2: Flye
```
cd /Path/WorkingDirectory/StrainName/
mkdir Flye && cd Flye

flye --nano-raw /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz --out-dir . -g 5m
```  
</br>

### Assembler 3: Raven
```
cd /Path/WorkingDirectory/StrainName/
mkdir Raven && cd Raven

raven /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz --graphical-fragment-assembly StrainName_ONTtrimfilt_raven.gfa

awk '/^S/{print ">"$2"\n"$3}' StrainName_ONTtrimfilt_raven.gfa > StrainName_ONTtrimfilt_raven_all.fasta
```  
</br>

### Assembler 4: miniasm/minipolish
```
cd /Path/WorkingDirectory/StrainName/
mkdir miniasm && cd miniasm

minimap2/minimap2 -x ava-ont -t8 /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz | gzip -1 > StrainName_ONTtrimfilt_minimap.paf.gz

miniasm/miniasm -f /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz StrainName_ONTtrimfilt_minimap.paf.gz > StrainName_ONTtrimfilt_miniasm.gfa

minipolish /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz StrainName_ONTtrimfilt_miniasm.gfa > StrainName_ONTtrimfilt_miniasmpolish.gfa

awk '/^S/{print ">"$2"\n"$3}' StrainName_ONTtrimfilt_miniasmpolish.gfa | fold > StrainName_ONTtrimfilt_miniasmpolish.fasta
```  
</br>

### Assembler 5: Unicycler
```
cd /Path/WorkingDirectory/StrainName/
mkdir Unicycler && cd Unicycler

unicycler-runner.py -l /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz -o .
```  
</br>

### Step 2: Assess contig over-circularization (i.e., terminal redundancy) by self-alignment with nucmer (part of MUMmer) ###
*The output from each assembler is a multi-fasta file, but each contig needs to be aligned to itself separately; thus, before running nucmer, each multi-fasta file is split so that each contig is a separate file.*
```
cd /Path/WorkingDirectory/StrainName/
mkdir nucmer && cd nucmer

mkdir Flye && cd Flye
awk '/^>/{s=++d".fasta"} {print > s}' /Path/WorkingDirectory/StrainName/Flye/assembly.fasta
for j in $(ls *.fasta); do nucmer --maxmatch --nosimplify --prefix=$j $j $j; done
for k in $(ls *.delta); do mummerplot --prefix=$k --png $k; show-coords -r $k > $k.coords; done

cd ..
mkdir miniasm && cd miniasm
awk '/^>/{s=++d".fasta"} {print > s}' /Path/WorkingDirectory/StrainName/miniasm/StrainName_ONTtrimfilt_miniasmpolish.fasta
for l in $(ls *.fasta); do nucmer --maxmatch --nosimplify --prefix=$l $l $l; done
for m in $(ls *.delta); do mummerplot --prefix=$m --png $m; show-coords -r $m > $m.coords; done

cd ..
mkdir Raven && cd Raven
awk '/^>/{s=++d".fasta"} {print > s}' /Path/WorkingDirectory/StrainName/Raven/StrainName_ONTtrimfilt_raven_all.fasta
for n in $(ls *.fasta); do nucmer --maxmatch --nosimplify --prefix=$n $n $n; done
for o in $(ls *.delta); do mummerplot --prefix=$o --png $o; show-coords -r $o > $o.coords; done

cd ..
mkdir Unicycler && cd Unicycler
awk '/^>/{s=++d".fasta"} {print > s}' /Path/WorkingDirectory/StrainName/Unicycler/assembly.fasta
for p in $(ls *.fasta); do nucmer --maxmatch --nosimplify --prefix=$p $p $p; done
for q in $(ls *.delta); do mummerplot --prefix=$q --png $q; show-coords -r $q > $q.coords; done

cd ..
mkdir NextDenovo && cd NextDenovo
awk '/^>/{s=++d".fasta"} {print > s}' /Path/WorkingDirectory/StrainName/NextDenovo/NextPolish/01_rundir/genome.nextpolish.fasta
for r in $(ls *.fasta); do nucmer --maxmatch --nosimplify --prefix=$r $r $r; done
for s in $(ls *.delta); do mummerplot --prefix=$s --png $s; show-coords -r $s > $s.coords; done
```  

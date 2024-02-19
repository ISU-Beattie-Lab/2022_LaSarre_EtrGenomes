# Trycycler reconciliation of trimmed assemblies to generate consensus assembly

## Requirements

Before performing this step you need to have trimmed each of your independent assemblies to remove any terminal redundancy from over-circularized contigs, as described [here](Remove_terminal_redundancy.md). Performing Trycycler reconciliation without trimming your contigs first is likely to cause problems during contig clustering.

### Software needed (indicated version is that used in the 2022 paper)
* `Trycycler` (v0.4.1) and all external tool dependencies (see https://github.com/rrwick/Trycycler/wiki/Software-requirements).  
</br>

## General notes

The reconciliation instructions here are based entirely on the Trycyler documentation available at https://github.com/rrwick/Trycycler/wiki. 

**_Please read the Trycyler documentation fully before performing the steps outlined below._**  
</br>

## Procedure

### Step 1:  Copy all trimmed assemblies (.fasta format) for a given genome into a single directory ```/Path/StrainName/assemblies/``` (in accordance with Trycycler recommendations).  
</br>

### Step 2: Run the following Trycycler steps for each genome:

#### 1. Cluster
```
trycycler cluster --assemblies /Path/WorkingDirectory/StrainName/assemblies/*.fasta --reads /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz --out_dir /Path/WorkingDirectory/StrainName/trycycler
```
*It may be necessary to discard select clusters at this stage; see Trycycler documentation for details.*  
</br>

#### 2. Reconcile
```
trycycler reconcile --max_indel_size 500 --reads /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz --cluster_dir /Path/WorkingDirectory/StrainName/trycycler/cluster_XXX
```
*Separately reconcile each good cluster by modifying the above command to replace "XXX" with the cluster ID.*  
</br>

#### 3. MSA
```
trycycler msa --cluster_dir trycycler/cluster_XXX
```
*Separately perform MSA on each reconciled cluster by modifying the above command to replace "XXX" with the cluster ID.*  
</br>

#### 4. Partition
```
trycycler partition --reads /Path/WorkingDirectory/StrainName/StrainName_ONTtrimfilt.fastq.gz --cluster_dirs /Path/WorkingDirectory/StrainName/trycycler/cluster_*
```  
</br>

#### 5. Consensus
```
trycycler consensus --cluster_dir /Path/WorkingDirectory/StrainName/trycycler/cluster_XXX
```
*Separately perform consensus on each partitioned cluster by modifying the above command to replace "XXX" with the cluster ID.*  
</br>


### Step 3: Concatenate the consensus sequence from each cluster into a single file

```
cat /Path/WorkingDirectory/StrainName/trycycler/cluster_*/7_final_consensus.fasta > /Path/WorkingDirectory/StrainName/trycycler/StrainName_ONTtrimfilt_5assemb_consensus.fasta
```
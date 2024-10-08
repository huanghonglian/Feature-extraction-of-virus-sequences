# Feature-extraction-of-virus-sequences
a universal algorithm for analyzing virus genome NGS data


### Data preparation

Before starting the project, we need to prepare the sequencing files in fastq format in the path *./data/fastqraw* . We also need the reference sequence file and the sample information file in *./public/* .  

Setting the file path:
```shell
project="/home/feature_extraction4seq/"#change to your path
fq=${project}/data/fastqraw/
public=${project}/public/
result=${project}/result/
mkdir ${project}/data/
mkdir $result
mkdir $fq

```


### Unzip sequence file

Place the FASTQ files in the *./data/fastqraw/* folder. Batch decompression of paired-end sequencing data files is performed as follows
```shell
cd ${project}/data/fastqraw/
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	gunzip ${project}/data/fastqraw/${id}_2.fastq.gz
	gunzip ${project}/data/fastqraw/${id}_1.fastq.gz
done
```

### Quality control-fastqc
The original  sequencing data was subjected to quality control using FastQC software, including assessing base quality, sequencing depth, GC content, and adapter content.

```shell
outdir=${project}//data/fastqc
mkdir $outdir
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	fastqc $fq/${id}_1.fastq $fq/${id}_2.fastq -o $outdir
done
```


### Cutadapt
If the results of quality control are unsatisfactory, further removal of sequence adapters is required. The trim_galore software was used to remove adapter sequences from both read1 and read2 of the paired-end sequencing data. The results are saved at *./data/clean/*. The quality control results after adapter trimming are saved in *./data/fastqc_clean*.

```shell
outdir_clean=${project}//data/fastqc_clean
mkdir $outdir_clean
```

Paired-End
```shell
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
  # QC
  #fastqc ${fq}/${id}_1.fastq  ${fq}/${id}_2.fastq -o 
  
  # cut adapters
   trim_galore --cores 2 -q 20 \
           --phred33 --stringency 3 --length 20 -e 0.1 \
           --paired ${fq}/${id}_1.fastq  ${fq}/${id}_2.fastq \
           --gzip -o ${project}/data/clean/
   
  #QC again
  fastqc  ${project}/data/clean/${id}_1_val_1.fq.gz  ${project}/data/clean/${id}_2_val_2.fq.gz -o $outdir_clean
done
```
Single-read
```shell
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
  # QC
  #fastqc ${fq}/${id}_1.fastq  ${fq}/${id}_2.fastq -o 
  
  # cut adapters
   trim_galore --cores 2 -q 20 \
           --phred33 --stringency 3 --length 20 -e 0.1 \
             ${fq}/${id}.fastq \ 
           --gzip -o ${project}/data/clean/
   
  #QC again
  fastqc  ${project}/data/clean/${id}_trimmed.fq.gz -o $outdir_clean
done
```
Unzip file
```shell
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	gunzip ${project}/data/clean/${id}_2_val_2.fq.gz
	gunzip ${project}/data/clean/${id}_1_val_1.fq.gz
done
```

### Reference sequence index construction
The Burrows-Wheeler Aligner (BWA) software was used to build an index file for the reference sequence.

```shell
cd ${public}/reference/
bwa index -a bwtsw ref_seq.fasta 
```


### Merge paired-end reads
The Fast Length Adjustment of Short Reads (FLASH) software was used to merge the quality-controlled short sequences from Read1 and Read2, generating longer fragment sequences. 

```shell
mkdir ${project}/data/fastqbind/
cd ${project}/data/fastqbind/
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	flash ${project}/data/clean/${id}_2_val_2.fq ${project}/data/clean/${id}_2_val_2.fq -p 33 -r 250 -f 500 -s 100 -o ${id}
done
cd ${project}
```




### Sequence alignment-bwa

For successfully merged reads, BWA was then used to align them with the reference sequence to extract the aligned sequences. For those unmerged reads from Read1 and Read2, BWA was used to align Read1 and Read2 separately with the reference sequence, and SAMtools software was subsequently conducted to extract the sequences that aligned at both ends. If there were duplicate regions in the sequences aligned at both ends, the duplicates were removed and the sequences were merged into one. On the contrary, Read1 and Read2 sequences were kept. 

BWA for successfully merged reads, and the results are saved in *./result/bwasam/*:
```shell
mkdir $result/bwasam/
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	bwa mem -M ${public}/reference/ref_seq.fasta ${project}/data/fastqbind/${id}.extendedFrags.fastq > $result/bwasam/${id}.sam
done
```

BWA for unmerged reads, and the results are saved in *./result/unbindsam/*.:
```shell
mkdir $result/unbindsam/
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	echo $id
	bwa mem -M ${public}/reference/ref_seq.fasta ${project}/data/fastqbind/${id}.notCombined_1.fastq ${project}/data/fastqbind/${id}.notCombined_2.fastq > $result/unbindsam/${id}_unbind.sam
done
```

Single-read
```shell
mkdir $result/bwasam/
tail -n +2 ${public}/sample/metadata.txt | cut -f 1 | while read id
do
	bwa mem -M ${public}/reference/ref_seq.fasta ${project}/data/clean/${id}_trimmed.fq > $result/bwasam/${id}.sam
done
```


### Integration of sequence matching results
This process will generate multiple intermediate files. The final results, which include the mapping information of reads to the reference sequence for each sample, will be saved in *./result/finaldata/seq/*.


```shell
# Extract the bwa results from *./result/bwasam/*
python3 deal.py
# Extract the bwa results from *./result/unbindsam/*. Skip if single-ended sequencing.
python3 dealunbind.py
python3 againunbind.py
#Integration of the extraction results
python3 dealall.py
```


### Sequence feature calculation
For sequence feature calculation, the nucleotide-level mutation and Shannon entropy algorithms were selected. The results will separately output feature matrices for sequencing depths of [10,30,50,75,100,150,200,250,300,350,400,450,500], saved in *./feature/data/sampledata*.

```shell
mkdir $project/feature/
#Calculate sequence features for each sample
python3 feature.py
#Merge sequence features from all samples.
python3 bind.py
```


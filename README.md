# Feature-extraction-of-virus-sequences
a universal algorithm for analyzing virus genome NGS data

### Data preparation

Before starting the project, we need to prepare the sequencing files in fastq format in the path *./data/fastqraw* . We also need the reference sequence file and the sample information file in *./public/* .  

Setting the file path:
```shell
project="/home/feature_extraction4seq/"
fq=${project}/data/
public=${project}/public/
result=${project}/result/
```


### Unzip sequence file

Batch decompression of paired-end sequencing data files is performed as follows
```shell
cd ${project}/data/fastqraw/
tail -n +2 ${project}/data/metadata.txt | cut -f 1 | while read id
do
	gunzip ${project}/data/fastqraw/${id}_2.fastq.gz
	gunzip ${project}/data/fastqraw/${id}_1.fastq.gz
done
```
or:
```shell
cd ${project}/data/fastqraw/
for i in *.fastq.gz;
do
	gunzip $i
done
```

### Quality control-fastqc
The original  sequencing data was subjected to quality control using FastQC software, including assessing base quality, sequencing depth, GC content, and adapter content.

```shell

outdir=${project}//data/fastqraw/fastq
log=${project}//data/fastqraw/log
for i in $fq/*.fastq;
do
  sample=${i%_clean*}
  sample=${sample##*/}
  fastqc -o $outdir $i > $log/$sample.log
done;
```


### Cutadapt
The Cutadapt software was used to remove adapter sequences from both read1 and read2 of the paired-end sequencing data. The results are saved at *./data/clean/*

```shell
mkdir ${fq}/clean/
tail -n +2 ${public}/sample/sampleinfo.txt | cut -f 1 | while read id
do
  # QC
  fastqc ${fq}/${id}_1.fq.gz  ${fq}/${id}_2.fq.gz
  
  # cut adapters
   trim_galore --cores 2 -q 20 \
           --phred33 --stringency 3 --length 20 -e 0.1 \
           --paired ${fq}/${id}.1.fq.gz  ${fq}/${id}.2.fq.gz \
           --gzip -o ${fq}/clean/
   
  #QC again
  fastqc  ${fq}/clean/${id}.*_1.fq.gz  ${fq}/clean/${id}.*_2.fq.gz
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
tail -n +2 ${public}/sample/sampleinfo.txt | cut -f 1 | while read id
do
	flash ${fq}${id}_1.fastq ${fq}${id}_2.fastq -p 33 -r 250 -f 500 -s 100 -o ${id}
done
```


### Sequence alignment-bwa
For successfully merged reads, BWA was then used to align them with the reference sequence to extract the aligned sequences. For those unmerged reads from Read1 and Read2, BWA was used to align Read1 and Read2 separately with the reference sequence, and SAMtools software was subsequently conducted to extract the sequences that aligned at both ends. If there were duplicate regions in the sequences aligned at both ends, the duplicates were removed and the sequences were merged into one. On the contrary, Read1 and Read2 sequences were kept. 

BWA for successfully merged reads, and the results are saved in *./result/bwasam/*:
```shell
mkdir $result/bwasam/
cd ${project}/data/fastqbind/
tail -n +2 ${public}/sample/sampleinfo.txt | cut -f 1 | while read id
do
	bwa mem -M ${public}/reference/sequence.fasta ${project}/data/fastqbind/${id}.extendedFrags.fastq > $result/bwasam/${id}.sam
done
```

BWA for unmerged reads, and the results are saved in *./result/unbindsam/*:
```shell
mkdir $result/unbindsam/
cd ${project}/data/fastqbind/
tail -n +2 ${public}/sample/sampleinfo.txt | cut -f 1 | while read id
do
	echo $id
	bwa mem -M ${public}/reference/sequence.fasta ${project}/data/fastqbind/${id}.notCombined_1.fastq ${project}/data/fastqbind/${id}.notCombined_2.fastq > $result/unbindsam/${id}_unbind.sam
done
```


### Integration of sequence matching results

```shell
# Extract the bwa results from *./result/bwasam/*
python3 deal.py
# Extract the bwa results from *./result/unbindsam/*
python3 dealunbind.py
python3 againunbind.py
#Integration of the extraction results
python3 dealall.py
```


### Sequence feature calculation

```shell
python3 myfeature.py
python3 bind.py
```

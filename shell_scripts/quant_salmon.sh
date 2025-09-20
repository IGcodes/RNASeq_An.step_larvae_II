#!/usr/bin/bash
#PBS -l nodes=1:ppn=10
#PBS -k doe

# Navigating to account on data node
cd /data/gunarathnai

# Set paths - Replace with the path to your Apptainer image
CONTAINER_PATH="/data/gunarathnai/Apptainer_containers/salmon_July2025.sif"

# Running for loop to iterate over the samples
for sID in `cat /data/gunarathnai/RNAseq_EthiopiaNUCI/Fastq_files/UCI_lab/UCIlab_samples.txt`

do

# Running indexing with Salmon
apptainer exec "$CONTAINER_PATH" salmon quant \
	-i /data/gunarathnai/RNAseq_EthiopiaNUCI/Ref_transcriptome/Anstep_UCI_V1.0_rna_index \
	-l A \
	-1 /data/gunarathnai/RNAseq_EthiopiaNUCI/Fastq_files/UCI_lab/${sID}_1.fq.gz \
	-2 /data/gunarathnai/RNAseq_EthiopiaNUCI/Fastq_files/UCI_lab/${sID}_2.fq.gz \
	-p 8 \
	--validateMappings \
	--gcBias \
	-o /data/gunarathnai/RNAseq_EthiopiaNUCI/quant_files/UCI_quants/${sID}_quant

done

## -i Transcriptome index
## -l tells salmon that it should automatically determine the library type of the sequencing reads (e.g. stranded vs. unstranded etc.)
## -1 and -2 specify paired read files
## -p number of threads
## --validateMappings Switch on selective alignment (default mode in current version)
## --gcBias To Salmon will enable it to learn and correct for fragment-level GC biases in the input data. Estimates a correction factor for systematic biases commonly present in RNA-seq data.
## -o defines output directory
## Ref1 - https://combine-lab.github.io/salmon/getting_started/
## Ref2 - https://salmon.readthedocs.io/en/latest/salmon.html
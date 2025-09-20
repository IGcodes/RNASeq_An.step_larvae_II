#!/usr/bin/bash
#PBS -l nodes=1:ppn=5
#PBS -k doe

# Navigating to the directory with Reference transcriptome
cd /data/gunarathnai/RNAseq_EthiopiaNUCI/Ref_transcriptome

# Running indexing with Salmon
apptainer exec /data/gunarathnai/Apptainer_containers/salmon_July2025.sif salmon index -t Anstep_UCI_V1.0_rna.fna.gz -i Anstep_UCI_V1.0_rna_index

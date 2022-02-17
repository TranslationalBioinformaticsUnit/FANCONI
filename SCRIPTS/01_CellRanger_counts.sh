#!/bin/bash

#SBATCH --time=24:00:00
#SBATCH -o /datos_2/FANCONI/SCRIPTS/jobs/%x_%j.out
#SBATCH -e /datos_2/FANCONI/SCRIPTS/jobs/%x_%j.err
#SBATCH --job-name=countfastq
#SBATCH -p allnodes
#SBATCH --nodes=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=2G
#SBATCH --nodelist=nodo04

# Perform quality control of sequencing reads

abort()
{
    echo >&2 '
***************
*** ABORTED ***
***************
'
    echo "An error occurred. Exiting..." >&2
    exit 1
}

set -e
trap 'abort' 0

# Get genome name

Work_dir=$1
FastqID=$2

hostname >&2
echo $0 >&2
printf START >&2; uptime >&2
date >&2


/opt_repository/cellranger-3.0.1/cellranger count --id=${FastqID}_job \
                                       --transcriptome=/datos_2/datos_temporal/datos/BMNiche/BMNiche_human/reference_genome_HUMAN/hs38/ \
                                       --fastqs=/datos_2/FANCONI/Public_data/${FastqID}/ \
                                       --sample=${FastqID}


printf END >&2; uptime >&2

trap : 0

echo >&2 '
************
*** DONE *** 
************
'
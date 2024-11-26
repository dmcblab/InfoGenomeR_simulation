#!/bin/bash

source $HOME\/miniforge3/etc/profile.d/conda.sh
conda activate tools
export PATH="/home/ylee4/miniforge3/envs/truvari/bin/:$PATH"
alias python='python3'


art_profiler_illumina PrecisionFDAv2 fastq fq.gz 40

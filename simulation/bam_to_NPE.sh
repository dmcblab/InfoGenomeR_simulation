#!/bin/bash

samtools view -h -F 2 -b simulated_sorted.bam >simulated_sorted_NPE.bam
samtools sort -n simulated_sorted_NPE.bam > simulated_sorted_NPE_sorted.bam
bamToFastq  -i simulated_sorted_NPE_sorted.bam  -fq NPE.fq1 -fq2 NPE.fq2

# Overview

- InfoGenomeR_simulation generates approximately 3000 germline SVs based on the phase 3 of the 1000 Genomes Project (by randomly selecting 1500 germline SVs from the reported SVs in phase 3 individuals
and randomly generating 1500 germline SVs) and approximately 200 somatic SVs based on 140 TCGA cancer genomes (by using 13 simple and complex operations).
- For each cancer genome, cancer purity is simulated at 60%, 75%, and 90%, and haplotype coverage is simulated at 3X, 5X, 10X, 15X, and 20X.
- Available for GRCh37 and GRCh38.

<p align="center">
    <img height="700" src="./doc/overview.png">
  </a>
</p>

# Requirements

- BWA-MEM
- ART
# Usage (generation)

- Selection of a 1000G individual and SV simulations.\
`./simulation/simulation_iteration_3_to_20_pre_simulation.sh`
- It generates 32 simulated genomes from simul1 to simul32. Select one of them, and move to the directory.\
`cd simul1`\
`./simulation_iteration_3_to_20_pre_simulation_selection.sh`
- In the folder, related files are generated as follows.
  - germline_SVs.results.*: a simulated germline genome (tab-seperated format; coordinate1, coordinate2, cumulative length, haplotype, and chromosome).
  - SVs.results.*: a simulated cancer genome (tab-seperated format; coordinate1, coordinate2, cumulative length, haplotype, and chromosome).
  - germ_simulated.fa: a simulated germline genome  in the fasta format, which is used for ART simulation.
  - simulated.fa: a simulated cancer genome in the fasta format, which is used for ART simulation.
  - true_SV_sets_dup_checked.refcoor.germline_initial_SVs_changed_size1000.int: a true germline SV set of DEL, DUP, and INV (>1kb and except ALU, LINE, and SVA insertion).
  - true_SV_sets_somatic_dup_checked.refcoor: a true somatic SV set.
  - germline_results.cnv.edited: a true germline CNV set.
  - somatic_results.cnv: a true somatic CNV set.
- Excecute the ART simulation (3X, 5X, 10X, 15X, and 20X haplotype coverage and 0.6, 0.75, and 0.9 cancer purity), and perform BWA-MEM mapping.\
`./simulation_iteration_art_3_to_20.sh`
- Releated folders and files are generated.
  - germ_simulated_sorted.bam: a simulated germline BAM (in the f3, f5, f10, f15, and f20 folders).
  - simulated_sorted.bam: a simulated somatic BAM (in p0.6, p0.75, and p0.9 in the f* folders).
# Simulated data
- NA12878-based simulated data for GRCh37 and GRCh38 are deposited in Zenodo.
  - copy_numbers.control: copy number segmentation profiles from the simulated germline genome by InfoGenomeR.
  - cp_norm: genome-binning read depths from the simulated cancer genome.
  - cp_norm_germ: genome-binning read depths from the simulated cancer genome.
  - coverage.txt: genome-binning read depths for JaBbA analysis.
  - delly.format: somatic SVs detected by DELLY2.
  - manta.format: somatic SVs detected by Manta.
  - novobreak.format: somatic SVs detected by novoBreak.
  - het_snps.format: heterozygous SNPs.
  - hom_snps.format: homozygous SNPs.
  - NPE.fq1 and NPE.fq2: non-propoerly paired reads.
  - simulated_genome: related files and folders generated by the simulation code. To generate the fasta files and BAM files (not required for InfoGenomeR inputs but optional), execute following scripts.\
  `perl results_to_fa.pl`\
  `perl germ_results_to_fa.pl`\
  `./simulation_iteration_art_3_to_20.sh`
   - InfoGenomeR_output: output files from InfoGenomeR.


# Snakemake install (beta)
```
curl -L https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -o Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda config --set channel_priority strict
conda activate snakemake
```
# InfoGenomeR simulation repository
```
git clone https://github.com/dmcblab/InfoGenomeR_simulation.git
InfoGenomeR_repo=${PWD}/InfoGenomeR_simulation
```
# Conda environment setting
```
snakemake --core all --use-conda InfoGenomeR_simulation_env
```
# Dataset download
```
snakemake --core all --use-conda InfoGenomeR_download
```
# InfoGenomeR simulation workflow (from an input SV list)
## Inputs
- SVs.txt
  - Large variants to be simulated.
  - tab separated: type, start, end, chromosome index.
  - should be sorted.
  - overlapping in a haplotype is not allowed.
    - supported types
      - del
      - del (arm_loss)
      - dup
      - aoh/loh
      - chr_loss
        - should occur lastly at the file, in the descending order of chromosome indices.
      - chr_gain
      - arm_gain
    - chromosome index is from 0 to 45
      - 0 to 22 for chr1-to-chrX hap1
      - 23 to 45 for chr1-to-chrX hap2
- Take an example in `examples/SVs.txt`
```
dup     33194399        76954315        0       0.2     0       dup
aoh     165490763       166930763       0       1       0       aoh
del     221706658       234994253       0       1       0       del
```
## Outputs
```
${workspace_dir}
|-- simulation_output
|   |-- germ_simulated.fa
|   |-- simulated.fa
|   `-- truth.vcf
|-- art_reads_output
|   |-- sample_simulated1.fq.gz
|   `-- sample_simulated2.fq.gz
```

## Workflow
### Make a workspace
```
# go to the InfoGenomeR repository.
cd ${InfoGenomeR_repo}

# make a workspace directory
workspace_dir=InfoGenomeR_workspace1
mkdir -p ${workspace_dir}

# link the reference in the workspace directory
ln -s ${PWD}/1000G_haplotypes ${workspace_dir}/1000G_haplotypes
```
### Simulation run
```
snakemake --core all --use-conda ${workspace_dir}/simulation_output --config subclonality=1.0 coverage=33 haplotype=NA12878 SV=SVs.txt
```

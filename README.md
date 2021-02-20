# Usage
- Selection of a 1000G individual and SV simulations.\
`./simulation/simulation_iteration_3_to_20_pre_simulation.sh`
- It generates 32 simulated genomes from simul1 to simul32. Select one of them, and move to the directory.\
`cd simul1`\
`./simulation_iteration_3_to_20_pre_simulation_selection.sh`
- In the folder, related files are generated as follows.
  - germline_SVs.results: a simulated germline genome (tab-seperated format; coordinate1, coordinate2, cumulative length, haplotype, and chromosome).
  - SVs.results: a simulated cancer genome (tab-seperated format; coordinate1, coordinate2, cumulative length, haplotype, and chromosome).
  - germ_simulated.fa: a simulated germline genome  in the fasta format, which is used for ART simulation.
  - simulated.fa: a simulated cancer genome in the fasta format, which is used for ART simulation.
  - true_SV_sets_dup_checked.refcoor.germline_initial_SVs_changed_size1000: a true germline SV set of DEL, DUP, and INV (>1kb and except ALU, LINE, and SVA insertion).
  - true_SV_sets_somatic_dup_checked.refcoor: a true somatic SV set.
  - germline_results.cnv.edited: a true germline CNV set.
  - somatic_results.cnv: a true somatic CNA set.
- Excecute the ART simulation (3X, 5X, 10X, 15X, and 20X haplotype coverage and 0.6, 0.75, and 0.9 cancer purity), and perform BWA-MEM mapping.\
`./simulation_iteration_art_3_to_20.sh`
- Releated folders and files are generated.
  - germ_simulated_sorted.bam: a simulated germline BAM (in the f3, f5, f10, f15, and f20 folders).
  - simulated_sorted.bam: a simulated somatic BAM (in p0.6, p0.75, and p0.9 in the f* folders).

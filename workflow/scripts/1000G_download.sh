#!/bin/bash
set -euo pipefail
DEST_DIR=$1
mkdir -p "$DEST_DIR"
cd $DEST_DIR

# Function to perform the download
do_work() {
  local id=$1
  echo "Starting download for chromosome $id"
  wget -q "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${id}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
  wget -q "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chr${id}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
  echo "Completed download for chromosome $id"
}

# Maximum number of parallel jobs (use all available processors)
MAX_JOBS=$(nproc --all)
SEMAPHORE="semaphore.tmp"

# Initialize the semaphore
mkfifo "$SEMAPHORE"
exec 3<>"$SEMAPHORE"  # Open file descriptor
for _ in $(seq "$MAX_JOBS"); do echo >&3; done

# Start multiple downloads with a limit on parallel jobs
for i in {1..22}; do
  # Wait for a slot to open
  read -u 3
  {
    do_work "$i"
    echo >&3  # Release the slot
  } &
done

# Handle chrX download separately
read -u 3
{
  echo "Starting download for chromosome X"
  wget -q "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
  wget -q "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi"
  echo "Completed download for chromosome X"
  echo >&3  # Release the slot
} &

# Wait for all background jobs to complete
wait

# Clean up
exec 3>&-
rm -f "$SEMAPHORE"

# rename the chrX vcf file for consistency with the others
ln -s 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz  1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
ln -s 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi  1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi

echo "All downloads completed successfully into $DEST_DIR."


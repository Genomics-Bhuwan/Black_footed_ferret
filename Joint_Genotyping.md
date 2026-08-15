#!/bin/bash
#SBATCH --job-name=BFF_JointGenotype
#SBATCH --output=logs/joint_genotype_%A.out
#SBATCH --error=logs/joint_genotype_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=125
#SBATCH --mem=105G
#SBATCH --time=68:00:00
#SBATCH --partition=shared
#SBATCH --account=bio260092
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL


module purge
module load gcc/11.2.0
module load samtools/1.12
module load gatk/4.1.8.1


REF="/anvil/scratch/x-bbist/Black_footed_ferret/Reference_Genome_assembly_BFF/8_Samples/GCF_022355385.1_MUSNIG.SB6536_genomic.fna"

WORKDIR="/anvil/scratch/x-bbist/Black_footed_ferret/Reference_Genome_assembly_BFF/8_Samples/BFF_Preprocessing/GATK/GVCFs"

OUTDIR="/anvil/scratch/x-bbist/Black_footed_ferret/Reference_Genome_assembly_BFF/8_Samples/BFF_Preprocessing/GATK/GVCFs/Joint_Genotyping_Hard_Filtering"

mkdir -p ${OUTDIR}
mkdir -p logs


COMBINED_GVCF="${OUTDIR}/all_samples_combined.g.vcf.gz"
RAW_VCF="${OUTDIR}/all_samples_joint_genotyped.vcf.gz"
FILTERED_VCF="${OUTDIR}/all_samples_joint_genotyped_filtered.vcf.gz"


echo "=============================================="
echo "Black-footed ferret GATK Joint Genotyping"
echo "Start: $(date)"
echo "=============================================="


echo ">>> Combining GVCFs..."

gatk --java-options "-Xmx100G" CombineGVCFs \
   -R ${REF} \
   --variant ${WORKDIR}/SRR35734985.g.vcf.gz \
   --variant ${WORKDIR}/SRR35734986.g.vcf.gz \
   --variant ${WORKDIR}/SRR35734987.g.vcf.gz \
   --variant ${WORKDIR}/SRR35734988.g.vcf.gz \
   --variant ${WORKDIR}/SRR35734989.g.vcf.gz \
   --variant ${WORKDIR}/SRR35734991.g.vcf.gz \
   --variant ${WORKDIR}/SRR35734992.g.vcf.gz \
   -O ${COMBINED_GVCF}


echo ">>> CombineGVCFs finished"


echo ">>> Running GenotypeGVCFs..."

gatk --java-options "-Xmx100G" GenotypeGVCFs \
   -R ${REF} \
   -V ${COMBINED_GVCF} \
   -O ${RAW_VCF}


echo ">>> Joint genotyping finished"


echo ">>> Running VariantFiltration..."

gatk --java-options "-Xmx100G" VariantFiltration \
   -R ${REF} \
   -V ${RAW_VCF} \
   --filter-expression "DP < 5" \
   --filter-name "DP5" \
   --filter-expression "FS > 60.0" \
   --filter-name "FS60" \
   --filter-expression "SOR > 3.0" \
   --filter-name "SOR3" \
   --filter-expression "MQ < 40.0" \
   --filter-name "MQ40" \
   --filter-expression "MQRankSum < -12.5" \
   --filter-name "MQRankSum-12.5" \
   --filter-expression "ReadPosRankSum < -8.0" \
   --filter-name "ReadPosRankSum-8" \
   -O ${FILTERED_VCF}


echo "=============================================="
echo "Finished successfully"
echo "Output:"
echo ${FILTERED_VCF}
echo "End: $(date)"
echo "=============================================="



#########################################################################################
##########################################################################################
#!/bin/bash
#SBATCH --job-name=BFF_GATK_HardFilter
#SBATCH --output=logs/gatk_hard_filter_%A.out
#SBATCH --error=logs/gatk_hard_filter_%A.err
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=100
#SBATCH --mem=80G
#SBATCH --time=50:00:00
#SBATCH --partition=shared
#SBATCH --account=bio260092
#SBATCH --mail-user=bistbs@miamioh.edu
#SBATCH --mail-type=BEGIN,END,FAIL

# ── 1. Clean Environment and Load Modules ─────────────────────────────────────
module purge
module load gcc/11.2.0
module load samtools/1.12
module load gatk/4.1.8.1

# ── 2. Set Up Input/Output Paths ──────────────────────────────────────────────
REF="/anvil/scratch/x-bbist/Black_footed_ferret/Reference_Genome_assembly_BFF/8_Samples/GCF_022355385.1_MUSNIG.SB6536_genomic.fna"

INPUT_VCF="/anvil/scratch/x-bbist/Black_footed_ferret/Reference_Genome_assembly_BFF/8_Samples/BFF_Preprocessing/GATK/GVCFs/Joint_Genotyping_Hard_Filtering/all_samples_joint_genotyped.vcf.gz"

OUTDIR="/anvil/scratch/x-bbist/Black_footed_ferret/Reference_Genome_assembly_BFF/8_Samples/BFF_Preprocessing/GATK/GVCFs/Joint_Genotyping_Hard_Filtering/Klaus_bcftools_variant_filtration"

mkdir -p "${OUTDIR}"
mkdir -p logs

# Intermediate temp files
TEMP_STEP1="${OUTDIR}/temp_step1_biallelic_autosomes.vcf.gz"
TEMP_STEP2="${OUTDIR}/temp_step2_masked_genotypes.vcf.gz"

# Final Output
OUTPUT_VCF="${OUTDIR}/all_samples_GATK_Klaus_filtered.vcf.gz"

echo "=========================================================="
echo "Starting Dedicated GATK-only Hard Filtering Pipeline"
echo "Input:  ${INPUT_VCF}"
echo "Output: ${OUTPUT_VCF}"
echo "Start:  $(date)"
echo "=========================================================="

# ── 3. Step 1: Subset to Biallelic Autosomal SNPs ─────────────────────────────
echo ">>> Step 1: Extracting biallelic autosomal SNPs..."

gatk --java-options "-Xmx35G" SelectVariants \
  -R "${REF}" \
  -V "${INPUT_VCF}" \
  -select-type SNP \
  --restrict-alleles-to BIALLELIC \
  -L NC_081557.1 -L NC_081558.1 -L NC_081559.1 -L NC_081560.1 -L NC_081561.1 -L NC_081562.1 \
  -L NC_081563.1 -L NC_081564.1 -L NC_081565.1 -L NC_081566.1 -L NC_081567.1 -L NC_081568.1 \
  -L NC_081569.1 -L NC_081570.1 -L NC_081571.1 -L NC_081572.1 -L NC_081573.1 -L NC_081574.1 \
  -O "${TEMP_STEP1}"


# ── 4. Step 2: Apply Site Filters & Convert Bad Genotypes to Missing ──────────
echo ">>> Step 2: Applying hard filters and masking low GQ/DP genotypes..."

gatk --java-options "-Xmx70G" VariantFiltration \
  -R "${REF}" \
  -V "${TEMP_STEP1}" \
  --filter-expression "QD < 2.0" --filter-name "QD2" \
  --filter-expression "FS > 60.0" --filter-name "FS60" \
  --filter-expression "SOR > 3.0" --filter-name "SOR3" \
  --filter-expression "MQ < 40.0" --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5" --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0" --filter-name "ReadPosRankSum-8" \
  --genotype-filter-expression "DP < 5" --genotype-filter-name "DP5" \
  --genotype-filter-expression "GQ < 20" --genotype-filter-name "GQ20" \
  --set-filtered-genotype-to-no-call \
  -O "${TEMP_STEP2}"


# ── 5. Step 3: Remove Filtered Sites & Enforce Max Missingness ────────────────
echo ">>> Step 3: Discarding flagged sites and removing high missingness..."

# --max-nocall-number 1 keeps only sites where at most 1 out of 7 samples is missing (./.)
gatk --java-options "-Xmx70G" SelectVariants \
  -R "${REF}" \
  -V "${TEMP_STEP2}" \
  --exclude-filtered \
  --max-nocall-number 1 \
  -O "${OUTPUT_VCF}"


# ── 6. Cleanup Temporary Files & Count Output ─────────────────────────────────
echo ">>> Cleaning up intermediate files..."
rm -f "${TEMP_STEP1}" "${TEMP_STEP1}.tbi"
rm -f "${TEMP_STEP2}" "${TEMP_STEP2}.tbi"

echo ">>> Calculating final SNP count..."
# Quick count that bypasses the VCF header lines
FINAL_COUNT=$(zgrep -v "^#" "${OUTPUT_VCF}" | wc -l)

echo "=========================================================="
echo "GATK Hard Filtering Completed Successfully!"
echo "Final Biallelic Autosomal SNPs remaining: ${FINAL_COUNT}"
echo "End: $(date)"
echo "=========================================================="

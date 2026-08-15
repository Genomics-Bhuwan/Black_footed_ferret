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

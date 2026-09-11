#!/bin/sh
#SBATCH --job-name=cowpea_varcall_4samps
#SBATCH --partition=batch          # ask your HPC docs for the right queue
#SBATCH --nodelist=compute05
#SBATCH --cpus-per-task=8
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=K.Mwangi@cgiar.org,ousodaniel@gmail.com

set -euo pipefail

# Setup: Proj Org; FS structure, Link resource directory (contain universal [shared] files: data, tools)
if [ ! -e ~/vacs-bioinfo ]; then
  ln -s /home/douso/vacs-bioinfo ~/vacs-bioinfo
  res_dir="$HOME/vacs-bioinfo/variant-calling"
  proj_dir="/var/scratch/global/$USER/projects/vacs-bioinfo/variant-calling"
  mkdir -p "${proj_dir}"/{alignments,annotation/snpeff_data,env,logs,qc/{fastq_raw,fastq_trim},raw_data/{fasta,fastq/trimmed,metadata,reference/{assembly,annotation}},scripts,variants}
else
  if [ -e ~/vacs-bioinfo ] && [ ! -L ~/vacs-bioinfo ]; then
    res_dir="$HOME/vacs-bioinfo/variant-calling"
    proj_dir="/var/scratch/global/$USER/vacs-bioinfo/variant-calling"
  fi
fi

# Setup: Load modules
module load miniforge3 2>/dev/null || echo "Your Sys Admin prefers a lean HPC, no system-wide conda/mamba support - figure it out" # check conda/mamba support
##mamba env create -f envs/requirements.yaml
##mamba activate varcall

# Ad-hoc setup for env: miniforge management
export MAMBA_ROOT_PREFIX="${HOME}"/vacs-bioinfo/.local/bin/miniforge3
source "${MAMBA_ROOT_PREFIX}"/etc/profile.d/mamba.sh # we will share a single environment - avoid redundancy; resource efficiency
mamba activate /var/scratch/global/douso/vacs/varcall/envs # env identified by path rather than name

# Ad-hoc setup for env: R
export PATH=~/vacs-bioinfo/.local/bin/R/bin:$PATH

# Setup: Threads
pct_processor_to_use=75 # the percentage of your local compute (PC) resources to commit for the analysis
threads=${SLURM_CPUS_PER_TASK:-$(( ($(nproc) * pct_processor_to_use + 50) / 100 ))}

## Download reference
#cd "${proj_dir}"/raw_data/reference
#wget -r -np -nH --cut-dirs=4 -A "*.fna.gz,*.gff3.gz" \
#  https://data.legumeinfo.org/Vigna/unguiculata/genomes/IT97K-499-35.gnm1.vQnBW/
#gunzip *.gz

## Get FASTQ from SRA
#cd "${proj_dir}"/raw_data/fastq
#while read -r srr; do
#  prefetch "srr"
#  vdb-validate "$srr" | tee a "${proj_dir}/logs/vdb-validate-on-srr-dump.txt" || continue
#  fasterq-dump --split-files --threads $threads "$srr"
#  gzip "${srr}"_1.fastq "${srr}"_2.fastq
#done < "${proj_dir}"/raw_data/metadata/SRR_Acc_List.txt

#  Sub-sample accessions
n_samples=4
n_reads=2000000

shuf -n "${n_samples}" "${res_dir}/raw_data/metadata/SRR_Acc_List.txt" |
while read -r samp_id; do
    echo "Sampling $samp_id"
    seqtk sample -s100 "${res_dir}"/raw_data/fastq/${samp_id}_1.fastq.gz $n_reads | gzip -c > "${proj_dir}"/raw_data/fastq/${samp_id}_sub_1.fastq.gz
    seqtk sample -s100 "${res_dir}"/raw_data/fastq/${samp_id}_2.fastq.gz $n_reads | gzip -c > "${proj_dir}"/raw_data/fastq/${samp_id}_sub_2.fastq.gz && \
    echo "${samp_id}" >> "${proj_dir}"/raw_data/metadata/SRR_Acc_List_sub.txt # a file to log our sampled accessions
done

# Quality assessment
mkdir -p "${proj_dir}"/qc/fastqc_raw
fastqc -t $threads -o "${proj_dir}"/qc/fastqc_raw "${res_dir}"/raw_data/fastq/*_sub*.fastq.gz
multiqc "${res_dir}"/qc/fastqc_raw -o "${proj_dir}"/qc/fastqc_raw

# Trimming
mapfile -t samples < "${proj_dir}"/raw_data/metadata/SRR_Acc_List_sub.txt

counter0=0
for sample in "${samples[@]}"; do
  ((counter0 += 1))
  echo "Currently trimming sample ${sample}... (${counter0}/${#samples[@]})"
  fastp \
    -i "${proj_dir}"/raw_data/fastq/${sample}_sub_1.fastq.gz -I "${proj_dir}"/raw_data/fastq/${sample}_sub_2.fastq.gz \
    -o "${proj_dir}"/raw_data/fastq/trimmed/${sample}_sub_1.trim.fastq.gz -O "${proj_dir}"/raw_data/fastq/trimmed/${sample}_sub_2.trim.fastq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --thread $threads \
    --json "${proj_dir}"/qc/trimmed/${sample}_sub.fastp.json \
    --html "${proj_dir}"/qc/trimmed/${sample}_sub.fastp.html
done

# Quality assessment after trimming
mkdir -p "${proj_dir}"/qc/fastqc_trim
fastqc -t $threads -o "${proj_dir}"/qc/fastqc_trim "${proj_dir}"/raw_data/fastq/trimmed/*.fastq.gz
multiqc "${proj_dir}"/qc/fastqc_trim -o "${proj_dir}"/qc/fastqc_trim

# Indexing
cd "${proj_dir}"/raw_data/reference
bwa index -p "${proj_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa \
"${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz
samtools faidx -o "${proj_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz.fai \
"${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz # requires a bgzip-compressed ref
gatk CreateSequenceDictionary -R "${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz \
-O "${proj_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.dict

# Align
counter1=0
for sample in "${samples[@]}"; do
  ((counter1 += 1))
  echo "Currently aligning sample ${sample}... (${counter1}/${#samples[@]})"
  bwa mem -t $threads \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:${sample}_lib1" \
    "${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa \
    "${proj_dir}"/raw_data/fastq/trimmed/${sample}_sub_1.trim.fastq.gz "${proj_dir}"/raw_data/fastq/trimmed/${sample}_sub_2.trim.fastq.gz \
    | samtools sort -@ $threads -o "${proj_dir}"/alignments/${sample}_sub.sorted.bam -
  samtools index "${proj_dir}"/alignments/${sample}_sub.sorted.bam
done

# Deduplication: By marking method
counter=0
for sample in "${samples[@]}"; do
  ((counter += 1))
  echo "Currently marking duplicates for sample ${sample}... (${counter}/${#samples[@]})"
  gatk MarkDuplicates \
    -I "${proj_dir}"/alignments/${sample}_sub.sorted.bam \
    -O "${proj_dir}"/alignments/${sample}_sub.dedup.bam \
    -M "${proj_dir}"/alignments/${sample}_sub.dup_metrics.txt
  samtools index "${proj_dir}"/alignments/${sample}_sub.dedup.bam
done

# Post-alignment QC
for sample in "${samples[@]}"; do
  samtools flagstat "${proj_dir}"/alignments/${sample}_sub.dedup.bam > "${proj_dir}"/qc/${sample}_sub.flagstat.txt
  mosdepth --by 10000 -t $threads "${proj_dir}"/qc/${sample}_sub "${proj_dir}"/alignments/${sample}_sub.dedup.bam
done

# Variant calling
counter2=0
for sample in "${samples[@]}"; do
  ((counter2 += 1))
  echo "Currently variant calling sample ${sample}... (${counter2}/${#samples[@]})"
  gatk HaplotypeCaller \
    -R "${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz \
    -I "${proj_dir}"/alignments/${sample}_sub.dedup.bam \
    -O "${proj_dir}"/variants/${sample}_sub.g.vcf.gz \
    -ERC GVCF \
    -L Vu03    # Focal analysis: restrict to chromosome of interest; drop -L for the full genome
done

# Joint-genotyping
gatk GenomicsDBImport \
  $(for s in $(< "${proj_dir}"/raw_data/metadata/SRR_Acc_List_sub.txt); do echo -V "${proj_dir}"/variants/${s}_sub.g.vcf.gz; done) \
  --genomicsdb-workspace-path "${proj_dir}"/variants/genomicsdb_Vu03 \
  -L Vu03

gatk GenotypeGVCFs \
  -R "${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz \
  -V gendb://"${proj_dir}"/variants/genomicsdb_Vu03 \
  -O "${proj_dir}"/variants/cowpea_panel_sub.Vu03.vcf.gz

# Hard-filtering
gatk SelectVariants -V "${proj_dir}"/variants/cowpea_panel_sub.Vu03.vcf.gz --select-type-to-include SNP \
  -O "${proj_dir}"/variants/cowpea_panel_sub.Vu03.snps.vcf.gz
gatk SelectVariants -V "${proj_dir}"/variants/cowpea_panel_sub.Vu03.vcf.gz --select-type-to-include INDEL \
  -O "${proj_dir}"/variants/cowpea_panel_sub.Vu03.indels.vcf.gz

gatk VariantFiltration -V "${proj_dir}"/variants/cowpea_panel_sub.Vu03.snps.vcf.gz \
  --filter-expression "QD < 2.0"                 --filter-name "QD2" \
  --filter-expression "FS > 60.0"                 --filter-name "FS60" \
  --filter-expression "MQ < 40.0"                 --filter-name "MQ40" \
  --filter-expression "MQRankSum < -12.5"         --filter-name "MQRankSum-12.5" \
  --filter-expression "ReadPosRankSum < -8.0"     --filter-name "ReadPosRankSum-8" \
  --filter-expression "SOR > 3.0"                 --filter-name "SOR3" \
  -O "${proj_dir}"/variants/cowpea_panel_sub.Vu03.snps.filtered.vcf.gz

gatk VariantFiltration -V "${proj_dir}"/variants/cowpea_panel_sub.Vu03.indels.vcf.gz \
  --filter-expression "QD < 2.0"                 --filter-name "QD2" \
  --filter-expression "FS > 200.0"                 --filter-name "FS200" \
  --filter-expression "MQ < 40.0"                 --filter-name "MQ40" \
  --filter-expression "MQRankSum < -20.0"         --filter-name "MQRankSum-20" \
  --filter-expression "ReadPosRankSum < -8.0"     --filter-name "ReadPosRankSum-8" \
  --filter-expression "SOR > 3.0"                 --filter-name "SOR3" \
  -O "${proj_dir}"/variants/cowpea_panel_sub.Vu03.indels.filtered.vcf.gz

# VCF merging
gatk MergeVcfs -I "${proj_dir}"/variants/cowpea_panel_sub.Vu03.snps.filtered.vcf.gz \
               -I "${proj_dir}"/variants/cowpea_panel_sub.Vu03.indels.filtered.vcf.gz \
               -O "${proj_dir}"/variants/cowpea_panel_sub.Vu03.filtered.vcf.gz

# Normalise VCF
bcftools norm -f "${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz \
  -m -any "${proj_dir}"/variants/cowpea_panel_sub.Vu03.filtered.vcf.gz -Oz \
  -o "${proj_dir}"/variants/cowpea_panel_sub.Vu03.norm.vcf.gz

# Population-level filtering: filter
vcftools --gzvcf "${proj_dir}"/variants/cowpea_panel_sub.Vu03.norm.vcf.gz \
  --max-missing 0.8 --minDP 5 --maf 0.05 --min-alleles 2 --max-alleles 2 \
  --recode --recode-INFO-all --out "${proj_dir}"/variants/cowpea_panel_sub.Vu03.final

# Sanity QC: Transition/Transversion ratios
bcftools stats "${proj_dir}"/variants/cowpea_panel_sub.Vu03.final.recode.vcf | grep "ts/tv" > "${proj_dir}"/qc/cowpea_panel_sub.Vu03.final.recode.ts-tv-ratios.txt

# SanityQC: Relatedness / Redundancy
plink --vcf "${proj_dir}"/variants/cowpea_panel_sub.Vu03.final.recode.vcf --make-bed --allow-extra-chr --out "${proj_dir}"/variants/cowpea_plink
plink --bfile "${proj_dir}"/variants/cowpea_plink --allow-extra-chr --pca 10 --out "${proj_dir}"/variants/cowpea_pca

cat > "${proj_dir}"/scripts/vcf-kinship-pca-plot.R <<EOF
library(ggplot2)
args <- commandArgs(trailingOnly = TRUE)
proj_dir <- args[1]

pca <- read.table(
    file.path(proj_dir, "variants", "cowpea_pca.eigenvec")
)

pca_plot <- ggplot(pca, aes(V3, V4)) +
geom_point() +
labs(x = "PC1", y = "PC2", title = "Cowpea panel structure (chr Vu03 SNPs; Sub-sampled reads)")

ggsave(
    file.path(proj_dir, "variants", "cowpea_pca.png"),
    plot = pca_plot,
    width = 7,
    height = 5,
    dpi = 600
)
EOF

Rscript "${proj_dir}"/scripts/vcf-kinship-pca-plot.R "${proj_dir}"

# Variants annotation: Build custom SnpEff DB
mkdir -p "${proj_dir}"/annotation/snpeff_data/Vunguiculata_540_v1.2
cp "${res_dir}"/raw_data/reference/assembly/Vunguiculata_540_v1.2.fa.gz "${proj_dir}"/annotation/snpeff_data/Vunguiculata_540_v1.2/sequences.fa.gz
cp "${res_dir}"/raw_data/reference/annotation/Vunguiculata_540_v1.2.gene.gff3.gz "${proj_dir}"/annotation/snpeff_data/Vunguiculata_540_v1.2/genes.gff.gz
cp "${res_dir}"/raw_data/reference/annotation/Vunguiculata_540_v1.2.protein.fa "${proj_dir}"/annotation/snpeff_data/Vunguiculata_540_v1.2/protein.fa
cp "${res_dir}"/raw_data/reference/annotation/Vunguiculata_540_v1.2.cds.fa "${proj_dir}"/annotation/snpeff_data/Vunguiculata_540_v1.2/cds.fa

cat >> "${proj_dir}"/annotation/snpEff.config <<EOF
Vunguiculata_540_v1.2.genome : Vunguiculata_540_v1.2
EOF

snpEff build -gff3 -v Vunguiculata_540_v1.2 -c "${proj_dir}"/annotation/snpEff.config -dataDir "${proj_dir}"/annotation/snpeff_data

# Variants annotation: Annotate
snpEff -v Vunguiculata_540_v1.2 -c "${proj_dir}"/annotation/snpEff.config \
  "${proj_dir}"/variants/cowpea_panel_sub.Vu03.final.recode.vcf \
  > "${proj_dir}"/annotation/cowpea_panel_sub.Vu03.annotated.vcf

# Interpreting SnpEff annotation: Subset impactful annotations
bcftools view -i 'INFO/ANN ~ "HIGH" || INFO/ANN ~ "MODERATE"' \
  "${proj_dir}"/annotation/cowpea_panel_sub.Vu03.annotated.vcf > "${proj_dir}"/annotation/cowpea_panel_sub.Vu03.highmod.vcf

# Interpreting SnpEff annotation: Intersecting with candidate region
zgrep "Vigun03g220400" "${res_dir}"/raw_data/reference/annotation/Vunguiculata_540_v1.2.gene.gff3.gz | awk '$3=="gene"' \
  | awk 'BEGIN{OFS="\t"}{print $1,$4-1,$5,"Vigun03g220400"}' > "${proj_dir}"/annotation/candidate_gene.bed

bedtools intersect -a "${proj_dir}"/annotation/cowpea_panel_sub.Vu03.highmod.vcf -b "${proj_dir}"/annotation/candidate_gene.bed -header \
  > "${proj_dir}"/annotation/candidate_gene.variants.vcf

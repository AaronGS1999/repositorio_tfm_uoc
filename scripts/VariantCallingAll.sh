#!/bin/bash

# alineamientos (Picard style y .bam)
BAM_FILES=("")

#configuración
REFERENCE= #Genoma de referencia en fasta
TMP_DIR= # Directorio para almacenar archivos temporales
THREADS= # Nº de hilos de CPU
MEMORY= # Cantidad de memoria RAM

samtools faidx $REFERENCE

for FILE in "${BAM_FILES[@]}"
do
  # Marcar duplicados y ordenar
  gatk MarkDuplicatesSpark --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -I $FILE \
      -O ${FILE%.bam}_sorted_dedup.bam \
      -M ${FILE%.bam}_dedup_metrics.txt
  
  # Métricas de alineamiento
  gatk CollectAlignmentSummaryMetrics --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE \
      -I ${FILE%.bam}_sorted_dedup.bam \
      -O ${FILE%.bam}_alignment_metrics.txt

  # Métricas del tamaño de inserción
  gatk CollectInsertSizeMetrics --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -I ${FILE%.bam}_sorted_dedup.bam \
      -O ${FILE%.bam}_insert_metrics.txt \
      -H ${FILE%.bam}_insert_size_histogram.pdf
  
  # Variant calling con HaplotypeCaller
  gatk HaplotypeCaller --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE \
      -I ${FILE%.bam}_sorted_dedup.bam \
      -O ${FILE%.bam}_raw_variants.vcf

  # Extracción de SNPs e INDELs
  gatk SelectVariants -R $REFERENCE -V ${FILE%.bam}_raw_variants.vcf \
      -select-type SNP -O ${FILE%.bam}_raw_snps.vcf
  gatk SelectVariants -R $REFERENCE -V ${FILE%.bam}_raw_variants.vcf \
      -select-type INDEL -O ${FILE%.bam}_raw_indels.vcf

  # Filtrado de SNPs
  gatk VariantFiltration --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE -V ${FILE%.bam}_raw_snps.vcf \
      -O ${FILE%.bam}_filtered_snps.vcf \
      -filter-name "Depth" -filter "DP < 10" \
      -filter-name "Freq" -filter "AF < 0.35" \
      -filter-name "QD_filter" -filter "QD < 2.0" \
      -filter-name "FS_filter" -filter "FS > 60.0" \
      -filter-name "MQ_filter" -filter "MQ < 40.0" \
      -filter-name "SOR_filter" -filter "SOR > 4.0" \
      -filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
      -filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0"

  # Filtrado de INDELs
  gatk VariantFiltration --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE -V ${FILE%.bam}_raw_indels.vcf \
      -O ${FILE%.bam}_filtered_indels.vcf \
      -filter-name "Depth" -filter "DP < 10" \
      -filter-name "Freq" -filter "AF < 0.35" \
      -filter-name "QD_filter" -filter "QD < 2.0" \
      -filter-name "FS_filter" -filter "FS > 200.0" \
      -filter-name "SOR_filter" -filter "SOR > 10.0"

  # Selección de SNPs e INDELs filtrados para BQSR
  gatk SelectVariants --exclude-filtered \
      -V ${FILE%.bam}_filtered_snps.vcf \
      -O ${FILE%.bam}_bqsr_snps.vcf
  gatk SelectVariants --exclude-filtered \
      -V ${FILE%.bam}_filtered_indels.vcf \
      -O ${FILE%.bam}_bqsr_indels.vcf

  # Recalibración de calidad de bases
  gatk BaseRecalibrator -R $REFERENCE -I ${FILE%.bam}_sorted_dedup.bam \
      --known-sites ${FILE%.bam}_bqsr_snps.vcf \
      --known-sites ${FILE%.bam}_bqsr_indels.vcf \
      -O ${FILE%.bam}_recal_data.table
  gatk ApplyBQSR -R $REFERENCE -I ${FILE%.bam}_sorted_dedup.bam \
      -bqsr ${FILE%.bam}_recal_data.table \
      -O ${FILE%.bam}_recal_reads.bam

  # Variant calling post-recalibración
  gatk HaplotypeCaller --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE \
      -I ${FILE%.bam}_recal_reads.bam \
      -O ${FILE%.bam}_raw_variants_recal.vcf

  # Extracción de SNPs e INDELs recalibrados
  gatk SelectVariants -R $REFERENCE -V ${FILE%.bam}_raw_variants_recal.vcf \
      -select-type SNP -O ${FILE%.bam}_raw_snps_recal.vcf
  gatk SelectVariants -R $REFERENCE -V ${FILE%.bam}_raw_variants_recal.vcf \
      -select-type INDEL -O ${FILE%.bam}_raw_indels_recal.vcf

  # Filtrado de SNPs e INDELs recalibrados
  gatk VariantFiltration --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE -V ${FILE%.bam}_raw_snps_recal.vcf \
      -O ${FILE%.bam}_filtered_snps_final.vcf \
      -filter-name "Depth" -filter "DP < 10" \
      -filter-name "Freq" -filter "AF < 0.01" \

  gatk VariantFiltration --java-options "-Xmx$MEMORY -Djava.io.tmpdir=$TMP_DIR -XX:ParallelGCThreads=$THREADS" \
      -R $REFERENCE -V ${FILE%.bam}_raw_indels_recal.vcf \
      -O ${FILE%.bam}_filtered_indels_final.vcf \
      -filter-name "QD_filter" -filter "QD < 2.0" \
      -filter-name "FS_filter" -filter "FS > 200.0" \
      -filter-name "SOR_filter" -filter "SOR > 10.0"
done
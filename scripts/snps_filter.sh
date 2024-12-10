#!/bin/bash

# Alineamientos (Picard style y .bam)
# Solo necesrio para los nombres de los archivos
BAM_FILES=("")

# Procesamiento de cada archivo
for FILE in "${BAM_FILES[@]}"
do
  # Filtrar variantes con bcftools (solo aquellas con FILTER = PASS)
  bcftools view -f PASS ${FILE%.bam}_filtered_snps_final.vcf -o ${FILE%.bam}_filtered_snps_clean.vcf
done
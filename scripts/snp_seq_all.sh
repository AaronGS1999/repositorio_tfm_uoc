#!/bin/bash

# Directorios
ALIGNMENTS_DIR= # Ruta a los alineamientos
SNPS_DIR= # Ruta a los archivos .vcf a filtrar
OUTPUT_DIR_BASE= # Ruta donde se almacenarán los resultados

# Archivos
ALIGNMENTS=(${ALIGNMENTS_DIR}/*.bam)  # Obtiene todos los archivos .bam
SNPS_FILES=(${SNPS_DIR}/*_filtered_snps_clean_filtered_snps_clean.vcf.gz)

# Profundidades a analizar
DEPTHS=(10 50 100)

# Iterar sobre las profundidades
for DEPTH in "${DEPTHS[@]}"; do

  # Directorio de salida para la profundidad actual
  OUTPUT_DIR="${OUTPUT_DIR_BASE}/DP${DEPTH}"

  # Crear directorio de salida si no existe
  mkdir -p $OUTPUT_DIR

  # Paso 0: Crear índices para archivos VCF.gz si no existen
  for vcf in "${SNPS_FILES[@]}"; do
    if [ ! -f "${vcf}.csi" ]; then
      echo "Generando índice para $vcf..."
      bcftools index -c $vcf
    fi
  done

  # Paso 1: Combinar SNPs de todas las variedades en un archivo BED
  ALL_SNPS_BED=${OUTPUT_DIR}/all_snps_unique.bed
  > $ALL_SNPS_BED

  for vcf in "${SNPS_FILES[@]}"; do
    bcftools query -f '%CHROM\t%POS0\t%POS\n' $vcf >> $ALL_SNPS_BED
  done

  # Eliminar duplicados en el archivo BED
  sort -u $ALL_SNPS_BED > $ALL_SNPS_BED.tmp && mv $ALL_SNPS_BED.tmp $ALL_SNPS_BED

  # Paso 2: Calcular la profundidad combinada
  COMBINED_DEPTH=${OUTPUT_DIR}/combined_depth.txt
  samtools depth -a -b $ALL_SNPS_BED "${ALIGNMENTS[@]}" > $COMBINED_DEPTH

  # Paso 3: Filtrar posiciones con profundidad >= DEPTH en todas las variedades
  VALID_POSITIONS_BED=${OUTPUT_DIR}/valid_positions.bed
  awk -v depth=$DEPTH '{
    valid=1; 
    for (i=3; i<=NF; i++) {
      if ($i < depth) {
        valid=0; 
        break;
      }
    }
    if (valid) print $1 "\t" $2 "\t" $2+1
  }' $COMBINED_DEPTH > $VALID_POSITIONS_BED

  # Paso 4: Filtrar SNPs en cada archivo VCF
  for vcf in "${SNPS_FILES[@]}"; do
    BASENAME=$(basename $vcf .vcf.gz)
    OUTPUT_VCF=${OUTPUT_DIR}/${BASENAME}_filtered.vcf.gz
    bcftools view -R $VALID_POSITIONS_BED $vcf -o $OUTPUT_VCF -O z
  done

  echo "Filtrado completado para profundidad $DEPTH. Archivos generados en $OUTPUT_DIR"

done
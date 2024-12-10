#!/bin/bash

# Crear la carpeta 'unicos' si no existe
mkdir -p unicos

# Comprimir todos los archivos .vcf a .vcf.gz si no están comprimidos
for file in *_filtered_snps_clean.vcf; do
    if [[ -f "$file" ]]; then
        bgzip -c "$file" > "$file.gz"
    fi
done

# Crear índices .tbi para todos los archivos .vcf.gz
for file in *.vcf.gz; do
    if [[ -f "$file" ]]; then
        tabix -p vcf "$file"
    fi
done

# Obtener una lista de archivos comprimidos y el total de archivos
compressed_files=(*_filtered_snps_clean.vcf.gz)
num_files=${#compressed_files[@]}

# Ejecutar bcftools isec con -w desde 1 hasta num_files y guardar los SNPs únicos en 'unicos'
for ((i=1; i<=num_files; i++)); do
    output_name="unicos/$(basename "${compressed_files[$i-1]}" .vcf.gz)_filtered_snps_clean.vcf.gz"
    bcftools isec -n=1 -w"$i" "${compressed_files[@]}" | bgzip -c > "$output_name"
done

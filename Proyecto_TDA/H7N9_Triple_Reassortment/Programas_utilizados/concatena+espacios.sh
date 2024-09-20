#!/bin/bash

# Directorio donde se encuentran los archivos .fasta
directorio="./carpeta"

# Archivo de salida
archivo_salida="final.txt"

# Eliminar el archivo de salida si ya existe
if [ -f "$archivo_salida" ]; then
  rm "$archivo_salida"
fi

# Recorrer cada archivo .fasta en el directorio
for archivo in "$directorio"/*.fasta; do
  # Verificar si el archivo existe
  if [ -e "$archivo" ]; then
    # Extraer el texto del archivo y agregarlo al archivo de salida
    cat "$archivo" >> "$archivo_salida"
    
    # Agregar un salto de línea al final del texto concatenado
    echo "" >> "$archivo_salida"
  fi
done

echo "Concatenación completada. El resultado se encuentra en $archivo_salida."
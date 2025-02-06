import os

# Lista de archivos a procesar
archivos = [
    "segment1_aligned.txt",
    "segment2_aligned.txt",
    "segment3_aligned.txt",
    "segment4_aligned.txt",
    "segment5_aligned_mark.txt",
    "segment6_aligned.txt",
    "segment7_aligned_mark.txt",
    "segment8_aligned_mark.txt"
]

# Diccionario para almacenar las secuencias
secuencias = {}

# Procesar cada archivo
for archivo in archivos:
    with open(archivo, "r") as f:
        contenido = f.read()
        # Dividir el contenido por el delimitador ">"
        partes = contenido.split(">")[1:]
        for parte in partes:
            # Obtener el nombre de la secuencia y la secuencia en sí
            nombre, secuencia = parte.strip().split("\n", 1)
            # Agregar la secuencia al diccionario
            if nombre not in secuencias:
                secuencias[nombre] = [""] * len(archivos)
            secuencias[nombre][archivos.index(archivo)] = secuencia.replace("\n", "")

# Escribir las secuencias en el archivo final.txt
with open("final.txt", "w") as f:
    for nombre, secuencia_list in secuencias.items():
        f.write(f">{nombre}\n")
        f.write("".join(secuencia_list) + "\n")
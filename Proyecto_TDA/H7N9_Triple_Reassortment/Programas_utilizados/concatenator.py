import re

def extract_name(header):
    match = re.search(r'\((.+?)\)', header)
    if match:
        return match.group(1)
    return None

# Leer las secuencias de los archivos
ha_sequences = {}
try:
    with open('HA_alligned.txt', 'r') as file:
        header = None
        seq = ''
        for line in file:
            if line.startswith('>'):
                if header:
                    ha_sequences[extract_name(header)] = seq
                header = line.strip()
                seq = ''
            else:
                seq += line.strip()
        if header:
            ha_sequences[extract_name(header)] = seq
except FileNotFoundError:
    print("Error: archivo HA_alligned.txt no encontrado")

na_sequences = {}
try:
    with open('NA_alligned.txt', 'r') as file:
        header = None
        seq = ''
        for line in file:
            if line.startswith('>'):
                if header:
                    na_sequences[extract_name(header)] = seq
                header = line.strip()
                seq = ''
            else:
                seq += line.strip()
        if header:
            na_sequences[extract_name(header)] = seq
except FileNotFoundError:
    print("Error: archivo NA_alligned.txt no encontrado")

pb2_sequences = {}
try:
    with open('PB2_alligned.txt', 'r') as file:
        header = None
        seq = ''
        for line in file:
            if line.startswith('>'):
                if header:
                    pb2_sequences[extract_name(header)] = seq
                header = line.strip()
                seq = ''
            else:
                seq += line.strip()
        if header:
            pb2_sequences[extract_name(header)] = seq
except FileNotFoundError:
    print("Error: archivo PB2_alligned.txt no encontrado")

# Encontrar las secuencias que tienen el mismo nombre
common_names = set(ha_sequences.keys()) & set(na_sequences.keys()) & set(pb2_sequences.keys())

# Escribir las secuencias concatenadas en el archivo final.txt
try:
    with open('final.txt', 'w') as file:
        for name in common_names:
            file.write(name + '\n')
            file.write(ha_sequences[name] + na_sequences[name] + pb2_sequences[name] + '\n\n')
except IOError:
    print("Error: no se pudo escribir en el archivo final.txt")

for name in ha_sequences.keys():
    print(name)

for name in na_sequences.keys():
    print(name)

for name in pb2_sequences.keys():
    print(name)

import numpy as np

temp_name=""

#Lee el archivo FASTA y manda los nombres y secuencias a un diccionario
with open('archivo.fasta') as f:
  sequences = {}
  name = None
  for line in f:
    line = line.strip()
    if line.startswith('>'):  
      name = line[1:]
      if (temp_name == ""):
        temp_name=name       
      sequences[name] = ''
    else:
      sequences[name] += line
 
add_spaces_in=0
add_spaces_fin=0
        #Encuentra la cantidad de espacios iniciales "-" para saber cuántas posiciones están alineadas antes de que empiece la secuencia real.
while sequences[temp_name][add_spaces_in]=='-': 
  add_spaces_in= add_spaces_in+1

while sequences[temp_name][len(sequences[temp_name])-1-add_spaces_fin]=='-': 
  add_spaces_fin= add_spaces_fin+1

print(len(sequences[temp_name]))
#print(sequences[temp_name][add_spaces_in])
#print(sequences[temp_name][len(sequences[temp_name])-add_spaces_fin])
print(str(add_spaces_in)+"    "+str(add_spaces_fin))

#Imprime los nombres de las secuencias encontradas.
print("Secuencias:")
for name in sequences.keys():
  print(name)
  
num_seqs = len(sequences)
dist_matrix = np.zeros((num_seqs, num_seqs))

for i in range(num_seqs):
  for j in range(i+1, num_seqs):
    # Usamos list(sequences.items()) para convertir el diccionario a una lista de tuplas (nombre, secuencia)
    name1, seq1 = list(sequences.items())[i]
    name2, seq2 = list(sequences.items())[j]
    diff_count = 0
    for k in range(0,len(seq1)):
      # Para cada par de secuencias, cuenta el número de diferencias (mutaciones) entre ellas en los rangos especificados.
      #if seq1[k] != seq2[k]:
      #if seq1[k] != seq2[k] and seq1[k]!= 'N' and seq2[k]!= 'N':
      if seq1[k] != seq2[k] and seq1[k]!= '-' and seq2[k]!= '-':
        diff_count += 1
      
      #Guarda las distancias en la matriz de distancia simétrica
      dist_matrix[i,j] = diff_count
      dist_matrix[j,i] = diff_count

print(dist_matrix)

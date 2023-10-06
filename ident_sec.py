#Jaime Salvador Lopez Viveros 06/10/23 
import re

with open('nombredelarchivo.txt') as f:
  lines = f.readlines()

#se crean los  vectores donde despues guardaremos las secuencias 
cured_lines = []
uncured_lines = [] 
uncured_names = [] 
cured_names = []

for line in lines:
  #si comienza con un > entonces es el nombre de la secuencia
  if line.startswith('>'):
    seq_name = line.strip()

  else:
    #si no es su contenido
    seq = line.strip()
    
    #verificamos que no haya cadenas del simbolo '-' cuyo tamaño sea multiplo de 3
    dash_groups = re.findall('-+', seq)
    invalid = any(len(g) % 3 != 0 for g in dash_groups)
    
    if invalid or re.search('[^ACTGN-]', seq):
      uncured_names.append(seq_name)
      uncured_lines.append(seq_name + '\n' + seq + '\n')
    else:
      cured_names.append(seq_name)
      cured_lines.append(seq_name + '\n' + seq + '\n')

with open('necesitacuracion.txt', 'w') as f:
  f.writelines(uncured_lines)
  
with open('nonecesitacuracion.txt', 'w') as f:
  f.writelines(cured_lines)

with open('nombresnecesitacuracion.txt', 'w') as f:
  f.writelines('\n'.join(uncured_names))
  
with open('nombresnonecesitacuracion.txt', 'w') as f:
  f.writelines('\n'.join(cured_names))

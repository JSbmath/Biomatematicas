import csv

# Abrir archivos
epicov_file = open('epicov.csv', 'r')
names_file = open('nombres.txt', 'r')
cur_file = open('cur.csv', 'w')

# Leer epicov.csv
epicov_reader = csv.reader(epicov_file)

# Leer nombres.txt y obtener subcadenas  
names = []
for line in names_file:
    parts = line.split('/')
    substring = parts[2] 
    names.append(substring)

# Escribir filas que coinciden en cur.csv
cur_writer = csv.writer(cur_file)
for row in epicov_reader:
    virus_name = row[0]
    if any(name in virus_name for name in names):
        cur_writer.writerow(row)
   
# Cerrar archivos        
epicov_file.close()
names_file.close()
cur_file.close()

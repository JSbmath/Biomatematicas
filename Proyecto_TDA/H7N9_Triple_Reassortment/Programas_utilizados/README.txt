Descripción de cada uno de los programas:

programa_matriz_dist_edit_genoma_new.py ---> Lee un archivo con terminación fasta llamado "archivo.fasta". Imprime las secuencias presentes en el archivo e imprime la matriz de distancia. 

programa_matriz_dist_edit_genoma_new_skip_part.py ---> Lo mismo que el anterior solo que no concidera las diferencias de las secuencias en el rango (6820-7636) de cada genoma 

concatenator.py ---> Toma los archivos llamados HA_alligned,NA_alligned y PB2_alligned. Encuentra las secuencias que provengan de la misma variante, las concatena y las imprime en un archivo llamado "final.txt"

concatena+espacios.sh ---> Guarda todas las secuencias que estén en la misma carpeta con terminación .fasta en un solo archivo llado "final.txt"

separa_partes_influenza.py ---> Separa las secuencias de los archivos presentes en la lista "file_list" y las guarda en "segment{i}.txt" dependiendo de la proteina a la cual pertenezca.

concatena_segmentos_formar_uno_solo.py ---> Concatena todos los segmentos presentes en los archivos "segment{i}.txt" de acuerdo a la variante a la que pertenezcan y las guarda en un archivo llamado "final.txt"

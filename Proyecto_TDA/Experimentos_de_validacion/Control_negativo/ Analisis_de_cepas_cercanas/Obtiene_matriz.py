import random
import os
import re
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import numpy as np

def normalize_sequence_strain(header):
    """
    Extrae y normaliza el nombre de la cepa de un encabezado
    """
    try:
        # Buscar el patrón (A/lugar/número/año(HxNx))
        match = re.search(r'\((A/[^)]+)\)', header)
        if match:
            return match.group(1)
        return None
    except Exception as e:
        print(f"Error normalizando cepa: {header}, Error: {str(e)}")
        return None

def check_required_sequences(file_path):
    """
    Verifica la presencia de las secuencias requeridas en el archivo
    """
    print("\nVerificando secuencias requeridas...")
    
    required_sequences = [
        "A/Nebraska/23/2009(H1N1)",
        "A/New York/0461/2009(H1N1)",
        "A/Wisconsin/629-D01038/2009(H1N1)",
        "A/California/VRDL133/2009(H1N1)",
        "A/Wisconsin/629-D00459/2009(H1N1)",  # Nueva secuencia agregada
        "A/New York/7236/2009(H1N1)",
        "A/California/VRDL121/2009(H1N1)"
    ]
    
    print("\nBuscando las siguientes secuencias:")
    for seq in required_sequences:
        print(f"- {seq}")
    
    sequences_by_strain = {}  # {strain_name: {segment_num: header}}
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    header = line.strip()[1:]
                    strain_name = normalize_sequence_strain(header)
                    if strain_name:
                        # Extraer número de segmento
                        segment_match = re.search(r'segment (\d+)', header.lower())
                        if segment_match:
                            segment_num = int(segment_match.group(1))
                            sequences_by_strain.setdefault(strain_name, {})[segment_num] = header
        
        # Verificar secuencias requeridas
        found_required = set()
        missing_sequences = set()
        
        for req_seq in required_sequences:
            found = False
            for strain in sequences_by_strain:
                if req_seq.replace('_', ' ') in strain or strain in req_seq.replace('_', ' '):
                    if len(sequences_by_strain[strain]) == 8:  # Verificar que tenga los 8 segmentos
                        found_required.add(strain)
                        found = True
                        break
            if not found:
                missing_sequences.add(req_seq)
        
        if missing_sequences:
            print("\nSecuencias requeridas no encontradas o incompletas:")
            for seq in missing_sequences:
                print(f"No está presente o incompleta: {seq}")
            return None
            
        print("\nTodas las secuencias requeridas fueron encontradas.")
        
        # Obtener solo las secuencias obligatorias
        selected_headers = []
        for strain in found_required:
            for segment_num in range(1, 9):
                if segment_num in sequences_by_strain[strain]:
                    selected_headers.append(sequences_by_strain[strain][segment_num])
        
        print(f"\nTotal de cepas seleccionadas: {len(found_required)}")
        print("\nCepas seleccionadas:")
        for strain in found_required:
            print(f"- {strain}")
        
        return selected_headers
        
    except FileNotFoundError:
        print(f"\nNo se encontró el archivo: {file_path}")
        return None
    except Exception as e:
        print(f"\nError al verificar secuencias: {str(e)}")
        print("Traza completa del error:")
        import traceback
        print(traceback.format_exc())
        return None

def count_initial_dashes(target_name):
    """
    Cuenta los guiones iniciales en la secuencia especificada
    """
    leading_gaps = 0
    sequence = None
    
    try:
        with open('segmento4final_aligned.fasta', 'r') as f:
            found = False
            current_sequence = []
            
            for line in f:
                if line.startswith('>'):
                    if found and current_sequence:
                        sequence = ''.join(current_sequence)
                        break
                    if target_name in line:
                        found = True
                    else:
                        found = False
                elif found:
                    current_sequence.append(line.strip())
            
            if found and current_sequence:
                sequence = ''.join(current_sequence)
        
        if sequence:
            for char in sequence:
                if char == '-':
                    leading_gaps += 1
                else:
                    break
            print(f"\nGaps encontrados en {target_name}: {leading_gaps}")
            return leading_gaps
        else:
            print(f"No se encontró la secuencia {target_name}")
            return 0
            
    except FileNotFoundError:
        print("No se encontró el archivo segmento4final_aligned.fasta")
        return 0
    except Exception as e:
        print(f"Error al procesar la secuencia: {str(e)}")
        return 0

def extract_segments(input_file, selected_sequences):
    """
    Procesa el archivo FASTA de secuencias de influenza para las secuencias seleccionadas
    """
    print("Paso 1: Extracción de segmentos")
    segment_identifiers = {
        'segment 1': 1,
        'PB2': 1,
        'segment 2': 2,
        'PB1': 2,
        'segment 3': 3,
        'PA': 3,
        'segment 4': 4,
        'HA': 4,
        'segment 5': 5,
        'NP': 5,
        'segment 6': 6,
        'NA': 6,
        'segment 7': 7,
        'matrix': 7,
        'segment 8': 8,
        'NS': 8
    }
    
    def get_segment_number(header):
        """Determina el número de segmento basado en el encabezado"""
        for i in range(1, 9):
            if f'segment {i}' in header.lower():
                return i
        
        for identifier, num in segment_identifiers.items():
            if identifier in header.lower() and 'segment' not in identifier:
                if identifier == 'matrix' and 'segment 8' in header.lower():
                    continue
                if identifier == 'NS' and 'segment 7' in header.lower():
                    continue
                return num
        
        return None

    # Procesar archivo
    segments = {i: [] for i in range(1, 9)}
    current_header = ''
    current_sequence = []
    sequence_found = False
    
    try:
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # Guardar secuencia anterior si fue seleccionada
                    if current_header and sequence_found and current_sequence:
                        segment_num = get_segment_number(current_header)
                        if segment_num:
                            segments[segment_num].append((current_header, ''.join(current_sequence)))
                    
                    # Preparar para nueva secuencia
                    current_header = line.strip()
                    current_sequence = []
                    sequence_found = any(seq in current_header for seq in selected_sequences)
                elif sequence_found:
                    current_sequence.append(line.strip())
            
            # Procesar última secuencia
            if current_header and sequence_found and current_sequence:
                segment_num = get_segment_number(current_header)
                if segment_num:
                    segments[segment_num].append((current_header, ''.join(current_sequence)))
        
        # Escribir segmentos en archivos
        for segment_num, sequences in segments.items():
            output_file = f'segmento{segment_num}.fasta'
            with open(output_file, 'w') as f:
                for header, sequence in sequences:
                    f.write(f'{header}\n')
                    f.write(f'{sequence}\n')
                    
    except Exception as e:
        print(f"Error procesando {input_file}: {str(e)}")
        return False
        
    return True

def concatenate_segments():
    """Concatena los archivos de segmentos en un solo archivo por segmento"""
    print("\nPaso 2: Concatenando segmentos...")
    for segment_num in range(1, 9):
        input_file = f'segmento{segment_num}.fasta'
        output_file = f'segmento{segment_num}final.fasta'
        
        try:
            if os.path.exists(input_file):
                with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                    outfile.write(infile.read())
        except Exception as e:
            print(f"Error concatenando segmento {segment_num}: {str(e)}")
            continue

def process_header_line(header):
    """Procesa una línea de encabezado FASTA reemplazando espacios por guiones bajos"""
    if not header.startswith('>'):
        return header
    
    paren_pos = header.rfind('(')
    
    if paren_pos == -1:
        return header.replace(' ', '_')
    else:
        before_paren = header[:paren_pos]
        after_paren = header[paren_pos:]
        return before_paren.replace(' ', '_') + after_paren

def replace_spaces_with_underscores():
    """Reemplaza espacios con guiones bajos en los encabezados"""
    print("\nPaso 3: Procesando encabezados...")
    for segment_num in range(1, 9):
        input_file = f'segmento{segment_num}final.fasta'
        output_file = f'segmento{segment_num}final_modified.fasta'
        
        if not os.path.exists(input_file):
            continue
            
        modified_lines = []
        with open(input_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    modified_lines.append(process_header_line(line.strip()) + '\n')
                else:
                    modified_lines.append(line)
        
        with open(output_file, 'w') as f:
            f.writelines(modified_lines)

def simplify_headers():
    """Simplifica los encabezados manteniendo solo el nombre de la cepa"""
    print("\nPaso 4: Simplificando encabezados...")
    for segment_num in range(1, 9):
        input_file = f'segmento{segment_num}final_modified.fasta'
        output_file = f'segmento{segment_num}final_modified_simple.fasta'
        
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        start_idx = line.find('(')
                        end_idx = line.find(')') + 1
                        if start_idx != -1 and end_idx != -1:
                            strain_name = line[start_idx-1:end_idx]
                            outfile.write(f'>{strain_name}\n')
                    else:
                        outfile.write(line)
                        
        except FileNotFoundError:
            continue

def remove_underscore_after_greater():
    """Elimina guiones bajos después del símbolo '>'"""
    print("\nPaso 5: Limpiando encabezados...")
    for segment_num in range(1, 9):
        input_file = f'segmento{segment_num}final_modified_simple.fasta'
        output_file = f'segmento{segment_num}final_modified_simple_no_underscore.fasta'
        sequence = ""
        
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                for line in infile:
                    if line.startswith('>'):
                        if sequence:
                            outfile.write(sequence)
                            sequence = ""
                        
                        if line.startswith('>_'):
                            line = '>' + line[2:]
                        
                        outfile.write(line)
                    else:
                        sequence += line
                
                if sequence:
                    outfile.write(sequence)
                    
        except FileNotFoundError:
            continue

def align_sequences():
    """Alinea las secuencias usando ClustalO"""
    print("\nPaso 6: Alineando secuencias...")
    for segment_num in range(1, 9):
        input_file = f'segmento{segment_num}final_modified_simple_no_underscore.fasta'
        output_file = f'segmento{segment_num}final_aligned.fasta'
        
        if not os.path.exists(input_file):
            continue
            
        try:
            import subprocess
            print(f"\nAlineando segmento {segment_num}...")
            cmd = ['clustalo', 
                   '-i', input_file,
                   '-o', output_file,
                   '--force',
                   '--auto']
            
            result = subprocess.run(cmd, 
                                 capture_output=True, 
                                 text=True)
            
            if result.returncode != 0:
                print(f"Error en el alineamiento del segmento {segment_num}")
                print(f"Error: {result.stderr}")
                continue
                
            print(f"Segmento {segment_num} alineado correctamente")
                
        except Exception as e:
            print(f"Error en el alineamiento del segmento {segment_num}: {str(e)}")
            continue

def concatenate_aligned_segments():
    """Concatena todos los segmentos alineados en un solo archivo"""
    print("\nPaso 7: Concatenando segmentos alineados...")
    archivos = [f'segmento{i}final_aligned.fasta' for i in range(1, 9)]
    secuencias = {}

    for archivo in archivos:
        try:
            current_sequence = ''
            current_header = ''
            
            with open(archivo, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        if current_header and current_sequence:
                            if current_header in secuencias:
                                secuencias[current_header] += current_sequence
                            else:
                                secuencias[current_header] = current_sequence
                        current_header = line.strip()[1:]
                        current_sequence = ''
                    else:
                        current_sequence += line.strip()
                
                if current_header and current_sequence:
                    if current_header in secuencias:
                        secuencias[current_header] += current_sequence
                    else:
                        secuencias[current_header] = current_sequence
            
        except FileNotFoundError:
            print(f"Archivo no encontrado: {archivo}")
            continue
    
    return secuencias

def calculate_distance_matrix(sequences):
    """Calcula la matriz de distancia entre las secuencias"""
    print("\nPaso 8: Calculando matriz de distancia...")
    
    # Obtener gaps iniciales para la secuencia de referencia
    leading_gaps = count_initial_dashes("A/wild_bird/Korea/A14/2011(H7N9)")
    
    # Calcular suma de los tres primeros segmentos
    range1_end = 2341 + 2358 + 2248  # Suma de segmentos 1, 2 y 3
    range2_start = range1_end + leading_gaps
    
    print(f"\nRangos de comparación:")
    print(f"Primera parte: 0-{range1_end}")
    print(f"Segunda parte: {range2_start}-final")
    
    print("\nSecuencias a analizar:")
    for name in sequences.keys():
        print(name)
    
    sequence_names = list(sequences.keys())
    num_seqs = len(sequences)
    dist_matrix = np.zeros((num_seqs, num_seqs))
    
    sequence_list = [sequences[name] for name in sequence_names]
    
    for i in range(num_seqs):
        for j in range(i+1, num_seqs):
            try:
                seq1 = sequence_list[i]
                seq2 = sequence_list[j]
                diff_count = 0
                
                # Contar diferencias en el primer rango (0 a suma de primeros 3 segmentos)
                for k in range(min(range1_end, len(seq1), len(seq2))):
                    if seq1[k] != seq2[k] and seq1[k] != '-' and seq2[k] != '-':
                        diff_count += 1
                
                # Contar diferencias en el segundo rango (después de los gaps)
                for k in range(range2_start, min(len(seq1), len(seq2))):
                    if seq1[k] != seq2[k] and seq1[k] != '-' and seq2[k] != '-':
                        diff_count += 1
                
                dist_matrix[i,j] = diff_count
                dist_matrix[j,i] = diff_count
                
            except Exception as e:
                print(f"Error en cálculo de matriz: {str(e)}")
                continue
    
    # Guardar solo la matriz numérica
    np.savetxt('matriz_resultante.txt', dist_matrix, fmt='%d')
    
    return dist_matrix

def cleanup_intermediate_files():
    """Elimina los archivos intermedios generados durante el proceso"""
    print("\nPaso 9: Limpiando archivos temporales...")
    patterns = [
        'segmento*.fasta',
        'segmento*final.fasta',
        'segmento*_modified.fasta',
        'segmento*_simple.fasta',
        'segmento*_no_underscore.fasta'
    ]
    
    import glob
    for pattern in patterns:
        for file_path in glob.glob(pattern):
            try:
                os.remove(file_path)
            except OSError:
                continue

def main():
    input_file = 'Influenza_A.fasta'

    try:
        print("\n=== Iniciando procesamiento de secuencias ===")
        
        # Verificar secuencias requeridas y seleccionar adicionales
        selected_sequences = check_required_sequences(input_file)
        if not selected_sequences:
            print("\nFinalizando programa debido a secuencias faltantes.")
            return
            
        print(f"\nProcesando {len(selected_sequences)} secuencias...")
        
        # Extraer segmentos del archivo
        if not extract_segments(input_file, selected_sequences):
            print("\nError en la extracción de segmentos. Finalizando programa.")
            cleanup_intermediate_files()
            return
        
        # Concatenar segmentos
        concatenate_segments()
        
        # Procesar encabezados
        replace_spaces_with_underscores()
        simplify_headers()
        remove_underscore_after_greater()
        
        # Alinear secuencias
        align_sequences()
        
        # Concatenar segmentos alineados y calcular matriz de distancia
        sequences = concatenate_aligned_segments()
        if sequences:
            calculate_distance_matrix(sequences)
            print("\nMatriz de distancia guardada en 'matriz_resultante.txt'")
        else:
            print("\nNo se encontraron secuencias para analizar")
            
        print("\n=== Procesamiento completado ===")
        
    except Exception as e:
        print(f"\nError durante la ejecución: {str(e)}")
        cleanup_intermediate_files()
        raise
        
    finally:
        # Limpieza final
        cleanup_intermediate_files()

if __name__ == "__main__":
    main()
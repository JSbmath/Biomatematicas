import random
import os
import re
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
import numpy as np

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

def extract_segments_multi(avian_files, human_files, individual_files, avian_samples=5, human_samples=2):  # Changed default values here
    """
    Procesa múltiples archivos FASTA de secuencias de influenza
    """
    print("Paso 1: Extracción de segmentos")
    # Mapeo de identificadores de segmentos a números
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
    
    def read_fasta(filename):
        """Lee un archivo FASTA y agrupa las secuencias por cepa"""
        sequences_by_strain = {}
        current_header = ''
        current_sequence = []
        
        with open(filename, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    if current_header:
                        sequences_by_strain.setdefault(strain_name, []).append(
                            (current_header, ''.join(current_sequence)))
                    current_header = line.strip()
                    current_sequence = []
                    strain_name = current_header[current_header.find('(')+1:current_header.find(')')]
                else:
                    current_sequence.append(line.strip())
            
            if current_header:
                sequences_by_strain.setdefault(strain_name, []).append(
                    (current_header, ''.join(current_sequence)))
        
        return sequences_by_strain

    def get_segment_number(header):
        """Determina el número de segmento basado en el encabezado"""
        for i in range(1, 9):
            if f'segment {i}' in header.lower():
                return i
        
        for identifier, num in segment_identifiers.items():
            if identifier in header and 'segment' not in identifier:
                if identifier == 'matrix' and 'segment 8' in header.lower():
                    continue
                if identifier == 'NS' and 'segment 7' in header.lower():
                    continue
                return num
        
        return None

    # Procesar archivos aviares
    print("\nProcesando archivos aviares...")
    for input_file, subtype in avian_files.items():
        segments = {i: [] for i in range(1, 9)}
        sequences_by_strain = read_fasta(input_file)
        
        num_samples = min(avian_samples, len(sequences_by_strain))  # Using new avian_samples parameter
        selected_strains = random.sample(list(sequences_by_strain.keys()), num_samples)
        
        for strain in selected_strains:
            for header, sequence in sequences_by_strain[strain]:
                segment_num = get_segment_number(header)
                if segment_num:
                    segments[segment_num].append((header, sequence))
        
        for segment_num, sequences in segments.items():
            output_file = f'segmento{segment_num}.fasta' if subtype == 'H9N2' else f'segmento{subtype}_{segment_num}.fasta'
            with open(output_file, 'w') as f:
                for header, sequence in sequences:
                    f.write(f'{header}\n')
                    f.write(f'{sequence}\n')

    # Procesar archivos humanos
    print("\nProcesando archivos humanos...")
    for input_file, subtype in human_files.items():
        segments = {i: [] for i in range(1, 9)}
        sequences_by_strain = read_fasta(input_file)
        
        num_samples = min(human_samples, len(sequences_by_strain))  # Using new human_samples parameter
        selected_strains = random.sample(list(sequences_by_strain.keys()), num_samples)
        
        for strain in selected_strains:
            for header, sequence in sequences_by_strain[strain]:
                segment_num = get_segment_number(header)
                if segment_num:
                    segments[segment_num].append((header, sequence))
        
        for segment_num, sequences in segments.items():
            output_file = f'segmento{segment_num}.fasta' if subtype == 'H9N2' else f'segmento{subtype}_{segment_num}.fasta'
            with open(output_file, 'a') as f:
                for header, sequence in sequences:
                    f.write(f'{header}\n')
                    f.write(f'{sequence}\n')
    
    # Procesar archivos individuales
    print("\nProcesando archivos individuales...")
    for input_file, subtype in individual_files.items():
        segments = {i: [] for i in range(1, 9)}
        with open(input_file, 'r') as f:
            current_header = ''
            current_sequence = []
            
            for line in f:
                if line.startswith('>'):
                    if current_header:
                        segment_num = get_segment_number(current_header)
                        if segment_num:
                            segments[segment_num].append((current_header, ''.join(current_sequence)))
                    current_header = line.strip()
                    current_sequence = []
                else:
                    current_sequence.append(line.strip())
            
            if current_header:
                segment_num = get_segment_number(current_header)
                if segment_num:
                    segments[segment_num].append((current_header, ''.join(current_sequence)))
        
        for segment_num, sequences in segments.items():
            output_file = f'segmento{segment_num}.fasta' if subtype == 'H9N2' else f'segmento{subtype}_{segment_num}.fasta'
            with open(output_file, 'a') as f:
                for header, sequence in sequences:
                    f.write(f'{header}\n')
                    f.write(f'{sequence}\n')

def concatenate_segments():
    """Concatena los archivos de segmentos de diferentes subtipos"""
    print("\nPaso 2: Concatenando segmentos...")
    for segment_num in range(1, 9):
        input_files = [
            f'segmento{segment_num}.fasta',         # H9N2
            f'segmentoH7N9_{segment_num}.fasta',    # H7N9
            f'segmentoH7N3_{segment_num}.fasta'     # H7N3
        ]
        
        output_file = f'segmento{segment_num}final.fasta'
        
        with open(output_file, 'w') as outfile:
            for input_file in input_files:
                try:
                    with open(input_file, 'r') as infile:
                        outfile.write(infile.read())
                        outfile.write('\n')
                except FileNotFoundError:
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
            continue
    
    return secuencias

def calculate_distance_matrix(sequences):
    """Calcula la matriz de distancia entre las secuencias"""
    print("\nPaso 8: Calculando matriz de distancia...")
    
    # Obtener gaps iniciales
    leading_gaps = count_initial_dashes("A/wild_bird/Korea/A14/2011(H7N9)")
    
    # Calcular suma de los tres primeros segmentos
    range1_end = 2341 + 2358 + 2248  # Suma de segmentos 1, 2 y 3
    range2_start = range1_end + leading_gaps
    
    print(f"\nRangos de comparación:")
    print(f"Primera parte: 0-{range1_end}")
    print(f"Segunda parte: {range2_start}-final")
    
    print("\nSecuencias:")
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
                for k in range(0, range1_end):
                    if k >= len(seq1) or k >= len(seq2):
                        raise IndexError(f"Índice {k} fuera de rango")
                    if seq1[k] != seq2[k] and seq1[k] != '-' and seq2[k] != '-':
                        diff_count += 1
                
                # Contar diferencias en el segundo rango (después de los gaps)
                for k in range(range2_start, min(len(seq1), len(seq2))):
                    if seq1[k] != seq2[k] and seq1[k] != '-' and seq2[k] != '-':
                        diff_count += 1
                
                dist_matrix[i,j] = diff_count
                dist_matrix[j,i] = diff_count
                
            except Exception as e:
                continue
    
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
        'segmento*_no_underscore.fasta',
        'final.txt'
    ]
    
    import glob
    for pattern in patterns:
        for file_path in glob.glob(pattern):
            try:
                os.remove(file_path)
            except OSError:
                continue

def main():
    # Definir los archivos de entrada
    avian_files = {
        '2010-2012_H9N2_257_avian.fasta': 'H9N2',
        '2010-2013_H7N3_13_avian.fasta': 'H7N3',
        '2011-2013_H7N9_32_avian.fasta': 'H7N9'
    }

    human_files = {
        'sequence_human(H9N2).fasta': 'H9N2',
        'sequence_human(H7N3).fasta': 'H7N3',
        'sequence_human(H7N9).fasta': 'H7N9'
    }

    individual_files = {
        'A-baikal_teal-Xianghai-421-2011(H9N2).fasta': 'H9N2',
        'A-brambling-Beijing-16-2012(H9N2).fasta': 'H9N2',
        'A-chicken-Zhejiang-329-2011_(H9N2).fasta': 'H9N2',
        'A-duck-Wuxi-7-2010(H9N2).fasta': 'H9N2',
        'A-duck-Zhejiang-2-2011(H7N3).fasta': 'H7N3',
        'A-duck-Zhejiang-10-2011(H7N3).fasta': 'H7N3',
        'A-quail-Lebanon-273-2010(H9N2).fasta': 'H9N2',
        'A-quail-wuxi-7-2010_(H9N2).fasta': 'H9N2',
        'A-Shanghai-02-2013(H7N9).fasta': 'H7N9',
        'A-wild_bird-Korea-A14-11(H7N9).fasta': 'H7N9'
    }

    try:
        print("\n=== Iniciando procesamiento de secuencias ===")
        # Using new sample sizes (5 avian, 2 human)
        extract_segments_multi(avian_files, human_files, individual_files, avian_samples=5, human_samples=2)
        concatenate_segments()
        replace_spaces_with_underscores()
        simplify_headers()
        remove_underscore_after_greater()
        align_sequences()
        sequences = concatenate_aligned_segments()
        
        if sequences:
            calculate_distance_matrix(sequences)
            
        print("\n=== Procesamiento completado ===")
        
    except Exception as e:
        print(f"\nError durante la ejecución: {str(e)}")
        cleanup_intermediate_files()
        raise
        
    finally:
        cleanup_intermediate_files()

if __name__ == "__main__":
    main()
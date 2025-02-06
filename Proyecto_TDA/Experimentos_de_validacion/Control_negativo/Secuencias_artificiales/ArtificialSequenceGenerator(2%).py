import random
from Bio import SeqIO
import os
from datetime import datetime
import re

class ArtificialSequenceGenerator:
    def __init__(self, key_sequences, file_mapping, mutation_rate=0.02):
        """
        Inicializa el generador de secuencias artificiales.
        
        Args:
            key_sequences: Lista de las secuencias que conforman el agujero bidimensional
            file_mapping: Diccionario que mapea nombres de secuencias a nombres de archivos
            mutation_rate: Tasa de mutación para generar variantes
        """
        self.key_sequences = key_sequences
        self.file_mapping = file_mapping
        self.mutation_rate = mutation_rate
        self.segment_data = {}
        self.nucleotides = ['A', 'T', 'G', 'C']

    def get_segment_number(self, header):
        """Determina el número de segmento del encabezado FASTA"""
        for i in range(1, 9):
            if f'segment {i}' in header.lower():
                return i
        return None

    def read_segments(self):
        """Lee y almacena todos los segmentos de las secuencias clave"""
        for sequence in self.key_sequences:
            self.segment_data[sequence] = {}
            filename = self.file_mapping[sequence]
            
            try:
                with open(filename, 'r') as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        segment_num = self.get_segment_number(record.description)
                        if segment_num:
                            self.segment_data[sequence][segment_num] = str(record.seq)
            except FileNotFoundError:
                print(f"No se encontró el archivo: {filename}")

    def generate_mutated_sequence(self, original_sequence):
        """Genera una secuencia mutada basada en la secuencia original"""
        mutated_seq = list(original_sequence)
        num_mutations = int(len(original_sequence) * self.mutation_rate)
        
        positions = random.sample(range(len(original_sequence)), num_mutations)
        for pos in positions:
            current_base = mutated_seq[pos]
            new_base = random.choice([b for b in self.nucleotides if b != current_base])
            mutated_seq[pos] = new_base
            
        return ''.join(mutated_seq)

    def modify_sequence_name(self, sequence_name, variant_num):
        """Modifica el nombre de la secuencia según el formato especificado"""
        # Separar la secuencia en partes
        parts = sequence_name.split('/')
        
        # Modificar la segunda parte (nombre del huésped)
        if len(parts) >= 2:
            parts[1] = f"{parts[1]}v{variant_num}"
        
        # Reconstruir el nombre
        return '/'.join(parts)

    def generate_artificial_variants(self):
        """Genera 5 variantes artificiales para cada secuencia clave"""
        self.read_segments()
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        output_file = f"artificial_sequences_{timestamp}.fasta"
        modifications_file = f"modifications_log_{timestamp}.txt"
        
        with open(output_file, 'w') as out_f, open(modifications_file, 'w') as log_f:
            for sequence in self.key_sequences:
                # Generar 5 variantes
                for variant_num in range(1, 6):
                    # Decidir aleatoriamente qué segmentos modificar
                    num_segments_to_modify = random.randint(1, 4)
                    segments_to_modify = random.sample(range(1, 9), num_segments_to_modify)
                    
                    # Generar el nuevo nombre de la secuencia
                    variant_name = self.modify_sequence_name(sequence, variant_num)
                    variant_name = f"{variant_name}_variant_{variant_num}"
                    
                    modification_info = f"\nVariante {variant_num} de {sequence}:\n"
                    modification_info += f"Segmentos modificados: {segments_to_modify}\n"
                    modification_info += f"Nuevo nombre: {variant_name}\n"
                    
                    # Generar y escribir cada segmento
                    for segment_num in range(1, 9):
                        if segment_num in self.segment_data[sequence]:
                            out_f.write(f">{variant_name} segment {segment_num}\n")
                            
                            if segment_num in segments_to_modify:
                                mutated_seq = self.generate_mutated_sequence(
                                    self.segment_data[sequence][segment_num]
                                )
                                out_f.write(f"{mutated_seq}\n")
                            else:
                                out_f.write(f"{self.segment_data[sequence][segment_num]}\n")
                    
                    log_f.write(modification_info)

def main():
    # Lista de secuencias que conforman el agujero bidimensional
    key_sequences = [
        'A/brambling/Beijing/16/2012(H9N2)',
        'A/quail/Wuxi/7/2010(H9N2)',
        'A/quail/Lebanon/273/2010(H9N2)',
        'A/Shanghai/02/2013(H7N9)',
        'A/wild bird/Korea/A14/2011(H7N9)',
        'A/duck/Zhejiang/10/2011(H7N3)'
    ]
    
    # Mapeo de nombres de secuencias a nombres de archivos
    file_mapping = {
        'A/brambling/Beijing/16/2012(H9N2)': 'A-brambling-Beijing-16-2012(H9N2).fasta',
        'A/quail/Wuxi/7/2010(H9N2)': 'A-quail-wuxi-7-2010_(H9N2).fasta',
        'A/quail/Lebanon/273/2010(H9N2)': 'A-quail-Lebanon-273-2010(H9N2).fasta',
        'A/Shanghai/02/2013(H7N9)': 'A-Shanghai-02-2013(H7N9).fasta',
        'A/wild bird/Korea/A14/2011(H7N9)': 'A-wild_bird-Korea-A14-11(H7N9).fasta',
        'A/duck/Zhejiang/10/2011(H7N3)': 'A-duck-Zhejiang-10-2011(H7N3).fasta'
    }
    
    # Crear el generador y generar las secuencias artificiales
    generator = ArtificialSequenceGenerator(key_sequences, file_mapping, mutation_rate=0.02)
    generator.generate_artificial_variants()

if __name__ == "__main__":
    main()
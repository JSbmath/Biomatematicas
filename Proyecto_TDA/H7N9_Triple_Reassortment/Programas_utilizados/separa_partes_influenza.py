import os
import re

def process_files(file_list):
    segment_data = {f"segment{i}.txt": [] for i in range(1, 9)}

    for file_path in file_list:
        with open(file_path, "r") as file:
            content = file.read()
            sequences = content.split(">")[1:]

            for seq in sequences:
                header, sequence = seq.split("\n", 1)
                name_match = re.search(r"\(.*?\)", header)
                segment_match = re.search(r"segment \d", header)

                if name_match and segment_match:
                    name = name_match.group()
                    segment = segment_match.group()[-1]

                    output_file = f"segment{segment}.txt"
                    segment_data[output_file].append((name, sequence.replace("\n", "")))
                else:
                    print(f"Skipping sequence: {header}")

    for output_file, data in segment_data.items():
        with open(output_file, "w") as file:
            for name, sequence in data:
                file.write(f">{name}\n{sequence}\n\n")

# Lista de archivos a procesar (puedes modificar esto según tus necesidades)
file_list = ["A-baikal_teal-Xianghai-421-2011(H9N2).fasta", "A-brambling-Beijing-16-2012.fasta", "A-chicken-Zhejiang-329-2011_(H9N2).fasta","A-duck-Wuxi-7-2010(H9N2).fasta","A-duck-Zhejiang-2-2011(H7N3).fasta","A-duck-Zhejiang-10-2011(H7N3).fasta","A-quail-Lebanon-273-2010(H9N2).fasta","A-quail-wuxi-7-2010_(H9N2).fasta", "A-Shanghai-02-2013(H7N9).fasta","A-wild_bird-Korea-A14-11(H7N9).fasta"]

process_files(file_list)
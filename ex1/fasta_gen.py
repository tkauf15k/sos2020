import random

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

min_length = 1000
genome_count = 500

input_path = "human_data_new.txt"
output_path = "human_data_new.fasta"

sequences = []
with open(input_path, 'r') as file:
    idx = 0
    for line in file:
        if len(line) < min_length:
            continue
        line = line.replace('\n', '').split()[0]
        sequences.append(SeqRecord(Seq(line), id="id:{}".format(idx), name="name:{}".format(idx), description="description:{}".format(idx)))
        idx += 1

        if idx >= genome_count:
            break

random.shuffle(sequences)

with open(output_path, "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")

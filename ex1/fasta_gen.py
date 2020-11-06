from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

min_length = 100
genome_count = 10000

input_path = "human_data.txt"
output_path = "human_data.fasta"

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

with open(output_path, "w") as output_handle:
    SeqIO.write(sequences, output_handle, "fasta")

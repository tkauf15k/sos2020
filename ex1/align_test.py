from Bio import AlignIO
from Bio import SeqIO
from Bio.Align import MultipleSeqAlignment

input_path = "human_data.fasta"
max_count = 10

with open("human_data_500.aln", "r") as input_handle:
    test = AlignIO.parse(input_handle,"clustal")
    a = list(test)

alignment = MultipleSeqAlignment([])

with open(input_path, "r") as input_handle:
    sequences = SeqIO.parse(input_handle, 'fasta')
    # sequences = AlignIO.read(input_handle, "fasta")
    idx = 0
    for sequence in sequences.records:
        alignment.add_sequence(sequence.description, str(sequence.seq))
        idx += 1
        if idx >= max_count:
            break

for record in alignment:
    print("%s %i" % (record.id, len(record)))

print('alignment-length:{}'.format(alignment.get_alignment_length()))  # equals the maximum

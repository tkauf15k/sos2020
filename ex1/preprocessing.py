import pickle
import sys

from lib.helperfunctions import preprocess

num_sequences = int(sys.argv[1])
input_path = "data/human_data_{}.fasta".format(num_sequences)

if __name__ == "__main__":
    cost_matrix = preprocess(input_path, num_sequences)
    pickle.dump(cost_matrix, open('data/human_data_{}.cm'.format(num_sequences), "wb"))

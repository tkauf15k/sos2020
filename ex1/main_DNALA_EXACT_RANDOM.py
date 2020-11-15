from Bio import SeqIO
import numpy as np
import timeit
import sys
from functools import lru_cache
from operator import itemgetter
from typing import List, Tuple
import random

from data.ExactWeightedMatching import ExactWeightedMatching
from lib.helperfunctions import preprocess, DNALA, fitness


def random_pairing(num_sequences: int):
    indizes = list(range(num_sequences))
    random.shuffle(indizes)
    pairing = []
    for i in range(num_sequences // 2):
        pairing.append((indizes[2 * i], indizes[2 * i + 1]))
    return pairing


if __name__ == "__main__":
    print("just some tests for reference algorithm and preprocessing...")

    num_sequences = int(sys.argv[1])
    input_path = "data/human_data_{}.fasta".format(num_sequences)

    cost_matrix = preprocess(input_path, num_sequences)
    print()

    print("---------- USING DNALA ----------")
    start = timeit.default_timer()
    pairing_dnala = DNALA(num_sequences, cost_matrix)
    runtime = timeit.default_timer() - start
    print("Fitness for DNALA: {0}, Runtime {1}s".format(fitness(pairing_dnala, cost_matrix), runtime))
    print()
    
    print("---------- USING EXACT ----------")
    start = timeit.default_timer()
    exact_pairing = ExactWeightedMatching(num_sequences, cost_matrix)
    runtime = timeit.default_timer() - start
    print("Fitness for EXACT: {0}, Runtime {1}s".format(fitness(exact_pairing, cost_matrix), runtime))
    print()

    print("---------- USING RANDOM ----------")
    print("Fitness for RANDOM: {0}".format(fitness(random_pairing(num_sequences), cost_matrix)))

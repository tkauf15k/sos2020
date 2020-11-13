from Bio import SeqIO
import numpy as np
import IUBAmbiguities
import timeit
from functools import lru_cache
from operator import itemgetter
from typing import List, Tuple
import random
import DNALA

from data.ExactWeightedMatching import ExactWeightedMatching


def find_gen(x, y, levels, hierarchy):
    if x == y:
        return x

    if x == ord('-') or y == ord('-'):  # just a shortcut...
        return ord('N')

    set_x = {x}
    set_y = {y}
    level_x = levels[x]
    level_y = levels[y]

    current_level = min(level_x, level_y)

    while current_level < levels[ord('N')]:

        if level_x == current_level:
            set_x = {b for a in set_x for b in hierarchy[a]}
            level_x += 1

        if level_y == current_level:
            set_y = {b for a in set_y for b in hierarchy[a]}
            level_y += 1

        xy = set_x.intersection(set_y)
        if len(xy) > 0:
            return list(xy)[0]

        current_level += 1

    return ord('N')


@lru_cache(maxsize=None)
def gen_cost(x, y):
    levels = IUBAmbiguities._levels
    hierarchy = IUBAmbiguities._hierarchy
    first = min(x, y)
    second = max(x, y)

    if first == second:
        return first, 0

    generalizer = find_gen(first, second, levels, hierarchy)
    cost = 2 * levels[generalizer] - levels[first] - levels[second]

    return generalizer, cost


def preprocess(input_path: str, num_sequences: int):
    with open(input_path, "r") as input_handle:
        sequences: List = list(SeqIO.parse(input_handle, 'fasta'))

    sequences: List[List[int]] = [[ord(c) for c in list(str(a.seq))] for a in sequences]

    data: np.array = np.array(sequences, dtype=np.int)

    snvr_indices = np.where((np.max(data, axis=0) - np.min(data, axis=0)) != 0)[0]

    snvrs = data[:, snvr_indices]

    sequence_length = snvrs.shape[1]

    # fill caches of all pairs
    for v, w in [(a, b) for a in IUBAmbiguities._alphabet for b in IUBAmbiguities._alphabet if ord(a) < ord(b)]:
        gen, cost = gen_cost(ord(v), ord(w))

    assert num_sequences == snvrs.shape[0]

    cost_matrix = np.zeros(shape=(num_sequences, num_sequences), dtype=int)

    start = timeit.default_timer()

    for a, b in [(a, b) for a in range(0, num_sequences) for b in range(0, num_sequences) if a < b]:
        gen_cost_vect_func = np.vectorize(lambda a, b: itemgetter(1)(gen_cost(a, b)))
        sum_cost = np.sum(gen_cost_vect_func(snvrs[a, :], snvrs[b, :]))
        cost_matrix[a, b] = cost_matrix[b, a] = sum_cost

    end = timeit.default_timer()

    print("preparation of costmatrix required: {} s".format(end - start))

    return cost_matrix


def fitness(pairing: List[Tuple[int, int]], cost_matrix: np.ndarray):
    return sum(cost_matrix[i][j] for i, j in pairing)


def fitness_tuple(pairing: List[Tuple[int, int]], cost_matrix: np.ndarray):
    '''
    ea lib requires tuple as return value...
    :param pairing:
    :param cost_matrix:
    :return:
    '''
    return (fitness(pairing, cost_matrix),)


def random_pairing(num_sequences: int):
    indizes = list(range(num_sequences))
    random.shuffle(indizes)
    pairing = []
    for i in range(num_sequences // 2):
        pairing.append((indizes[2 * i], indizes[2 * i + 1]))
    return pairing


def validate(individual):
    nodes = [i for a in individual for i in a]
    assert len(set(nodes)) == len(nodes)


if __name__ == "__main__":
    print("just some tests for reference algorithm and preprocessing...")

    num_sequences = 20  # 20 50 100
    input_path = "data/human_data_{}.fasta".format(num_sequences)

    cost_matrix = preprocess(input_path, num_sequences)
    pairing = DNALA(num_sequences, cost_matrix)
    exact_pairing = ExactWeightedMatching(num_sequences, cost_matrix)

    assert len(pairing) == len(exact_pairing)

    print(pairing)

    print(fitness(pairing, cost_matrix))
    print(fitness(exact_pairing, cost_matrix))
    print(fitness(random_pairing(num_sequences), cost_matrix))

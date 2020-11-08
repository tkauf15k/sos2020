from Bio import SeqIO
import numpy as np
from itertools import tee
import IUBAmbiguities
import timeit
from functools import lru_cache
from operator import itemgetter
from typing import Set, List, Tuple
import random

num_sequences = 20  # 20 50 100
input_path = "data/human_data_{}.fasta".format(num_sequences)


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

def preprocess(input_path: str):
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
        # print('{} : {} => {}, {}'.format(v, w, chr(gen), cost))

    assert num_sequences == snvrs.shape[0]

    cost_matrix = np.zeros(shape=(num_sequences, num_sequences), dtype=int)

    start = timeit.default_timer()

    for a, b in [(a, b) for a in range(0, num_sequences) for b in range(0, num_sequences) if a < b]:

        gen_cost_vect_func = np.vectorize(lambda a,b: itemgetter(1)(gen_cost(a, b)))
        sum_cost = np.sum(gen_cost_vect_func(snvrs[a, :], snvrs[b, :]))
        cost_matrix[a, b] = cost_matrix[b, a] = sum_cost

    end = timeit.default_timer()
    
    # print("just some arbitrary output: {}".format(np.average(cost_matrix)))
    print("preparation of costmatrix required: {} s".format(end - start))

    return cost_matrix

def DNALA(num_sequences: int, cost_matrix: np.ndarray):
    def get_min_indizes(cost_matrix: np.ndarray, i: int, S: Set[int]) -> Set[int]:
        row = cost_matrix[i, :]
        possible_ind = S.difference({i})
        
        mn = np.amin(np.take(row, list(possible_ind)))
        ind = np.where(row == mn)[0]
        return set(ind).intersection(possible_ind)

    start = timeit.default_timer()
    P = []
    S = set(range(num_sequences))

    while len(S) > 0:
        for s in np.random.permutation(list(S)):
            if s not in S:
                continue

            min_indizes_s = get_min_indizes(cost_matrix, s, S)

            for c in np.random.permutation(list(min_indizes_s)):
                min_indizes_c = get_min_indizes(cost_matrix, c, S)
                if s in min_indizes_c:
                    P.append((s,c))
                    S.remove(s)
                    S.remove(c)
                    break
    return P

def fitness(pairing: List[Tuple[int,int]], cost_matrix: np.ndarray):
    return sum(cost_matrix[i][j] for i,j in pairing)

def random_pairing(num_sequences: int):
    indizes = list(range(num_sequences))
    random.shuffle(indizes)
    pairing = []
    for i in range(num_sequences//2):
        pairing.append((indizes[2*i], indizes[2*i+1]))
    return pairing

if __name__ == "__main__":
    cost_matrix = preprocess(input_path)
    pairing = DNALA(num_sequences, cost_matrix)

    print(fitness(pairing, cost_matrix))
    print(fitness(random_pairing(num_sequences), cost_matrix))


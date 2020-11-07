from Bio import SeqIO
import numpy as np
from itertools import tee
import IUBAmbiguities

num_sequences = 10  # 20 50 100
input_path = "data/human_data_{}.fasta".format(num_sequences)

with open(input_path, "r") as input_handle:
    sequences = list(SeqIO.parse(input_handle, 'fasta'))

sequences = [[ord(c) for c in list(str(a.seq))] for a in sequences]

data = np.array(sequences, dtype=np.int)

snvr_indices = np.where((np.max(data, axis=0) - np.min(data, axis=0)) != 0)[0]
snvrs = data[:, snvr_indices]

_cache = {}


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


def gen_cost(x, y, _cache, levels, hierarchy):
    first = min(x, y)
    second = max(x, y)

    if (first, second) in _cache:
        return _cache[(first, second)]

    if first == second:
        return 0

    generalizer = find_gen(first, second, levels, hierarchy)
    cost = 2 * levels[generalizer] - levels[first] - levels[second]

    _cache[(first, second)] = (generalizer, cost)

    return generalizer, cost


alphabet = ['A', 'C', 'G', 'T', 'R', 'Y', 'S', 'W', 'M', 'K', 'B', 'D', 'H', 'V', '-', 'N']

# fill caches of all pairs
for v, w in [(a, b) for a in alphabet for b in alphabet if ord(a) < ord(b)]:
    gen, cost = gen_cost(ord(v), ord(w), _cache, IUBAmbiguities._levels, IUBAmbiguities._hierarchy)
    print('{} : {} => {}, {}'.format(v, w, chr(gen), cost))

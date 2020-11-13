import numpy as np
from typing import Set


def DNALA(num_sequences: int, cost_matrix: np.ndarray):
    def get_min_indizes(cost_matrix: np.ndarray, i: int, S: Set[int]) -> Set[int]:
        row = cost_matrix[i, :]
        possible_ind = S.difference({i})

        mn = np.amin(np.take(row, list(possible_ind)))
        ind = np.where(row == mn)[0]
        return set(ind).intersection(possible_ind)

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
                    P.append((s, c))
                    S.remove(s)
                    S.remove(c)
                    break
    return P
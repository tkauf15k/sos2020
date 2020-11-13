import numpy as np


def random_matching_initializer(container, func, n):
    perm = np.random.permutation(range(0, n))

    num_pairs = int(np.ceil(n / 2))

    parts = [(x[0], x[1]) for x in np.split(perm, num_pairs)]

    return func(parts)
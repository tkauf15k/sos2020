import sys

import numpy as np
import random

import DNAnon


def prepare_pairs(individual, indpb, toolbox):
    probabilities = np.random.random(size=len(individual))

    condition = probabilities <= 2 * indpb

    pairs = [idx for idx, elem in enumerate(condition) if elem == True]
    random.shuffle(pairs)

    if pairs is None or len(pairs) < 2:  # for small instances, we sometimes cannot mutate properly
        return (False, toolbox.clone(individual),)

    return (True, pairs)

def select_random_edge(p1, p2, cost_matrix):
    indices = np.random.randint(0, 2, 2)

    al1 = (p1[indices[0]], p2[indices[1]])
    al2 = (p1[1 - indices[0]], p2[1 - indices[1]])
    return al1, al2

def select_best_edge(p1,p2, cost_matrix):

    best = sys.maxsize
    best_combination = None
    for i,j in [(i,j) for i in range(2) for j in range(2)]:
        al1 = (p1[i], p2[j])
        al2 = (p1[1 - i], p2[1 - j])

        cost = cost_matrix[al1[0], al1[1]] + cost_matrix[al2[0], al2[1]]
        if cost < best:
            best = cost
            best_combination = (al1,al2)


    return best_combination



def two_opt_mutation(individual, indpb, toolbox, edge_selector, cost_matrix):
    success, pairs = prepare_pairs(individual, indpb, toolbox)
    if not success:
        return (pairs,)

    truncated_length = int(len(pairs) / 2) * 2

    cloned_individual = toolbox.clone(individual)
    for i, j in [(i, (i + 1)) for i in range(truncated_length - 1)]:

        al1, al2 = edge_selector(cloned_individual[pairs[i]], cloned_individual[pairs[j]], cost_matrix)

        cloned_individual[pairs[i]] = al1
        cloned_individual[pairs[j]] = al2

    DNAnon.validate(cloned_individual)
    return cloned_individual,

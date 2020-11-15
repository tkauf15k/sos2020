from collections import Set

import numpy as np
import random
from .helperfunctions import validate

import sys



def arbitrary_matched_subgraph(ind1, ind2, toolbox):
    c1 = create_arbitrary_matched_subgraph_offsprint(ind1, ind2, toolbox)
    c2 = create_arbitrary_matched_subgraph_offsprint(ind2, ind1, toolbox)

    validate(c1)
    validate(c2)
    return c1, c2


def create_arbitrary_matched_subgraph_offsprint(ind1, ind2, toolbox):
    ind1_indices = [idx for idx, a in enumerate(np.random.randint(0, 2, len(ind1))) if a == 1]
    pairs = {ind1[i] for i in ind1_indices}

    remainder_edges = {a for a in ind1 if a not in pairs}
    remainder_nodes = {i for a in remainder_edges for i in a}

    ind2_matching_arcs = [idx for idx, (a, b) in enumerate(ind2) if a in remainder_nodes and b in remainder_nodes]
    matching_arcs = [ind2[i] for i in ind2_matching_arcs]
    unmatched_nodes = [(a if a in remainder_nodes else b) for idx, (a, b) in enumerate(ind2) if
                       idx not in ind2_matching_arcs and (a in remainder_nodes or b in remainder_nodes)]

    assert len(unmatched_nodes) % 2 == 0
    random.shuffle(unmatched_nodes)
    new_edges = [(unmatched_nodes[i], unmatched_nodes[i + 1]) for i in range(0, len(unmatched_nodes) - 1, 2)] + matching_arcs

    remaining_indices = list(set(range(0, len(ind1))) - set(ind1_indices))
    freeIndices = 0
    new_arc_indices = 0

    c1 = toolbox.clone(ind1)

    for i in range(0, len(ind1)):
        if ind1[i] in pairs:
            continue

        next_free_index = remaining_indices[freeIndices]
        freeIndices += 1
        c1[next_free_index] = new_edges[new_arc_indices]
        new_arc_indices+=1

    return c1


def ordered_matching_best_match(ind1, ind2, toolbox):
    c1 = create_offspring(ind1, ind2, toolbox)
    c2 = create_offspring(ind2, ind1, toolbox)

    validate(c1)
    validate(c2)

    return c1, c2


def create_offspring(ind1, ind2, toolbox):
    split_index = np.random.randint(0, len(ind1))
    offspr = ind1[0:split_index]
    remainder_edges = ind1[split_index:]
    remainder_nodes = {i for a in ind1[split_index:] for i in a}
    ind2_matching_arcs = [idx for idx, (a, b) in enumerate(ind2) if a in remainder_nodes and b in remainder_nodes]
    matching_arcs = [ind2[i] for i in ind2_matching_arcs]
    unmatched_nodes = [(a if a in remainder_nodes else b) for idx, (a, b) in enumerate(ind2) if
                       idx not in ind2_matching_arcs and (a in remainder_nodes or b in remainder_nodes)]
    assert len(unmatched_nodes) % 2 == 0
    random.shuffle(unmatched_nodes)
    new_arcs = [(unmatched_nodes[i], unmatched_nodes[i + 1]) for i in range(0, len(unmatched_nodes) - 1, 2)]
    ind = offspr + matching_arcs + new_arcs
    assert len(offspr) + len(remainder_edges) == len(ind1)
    assert len(ind) == len(ind1)
    c1 = toolbox.clone(ind1)
    for i in range(len(c1)):
        c1[i] = ind[i]
    return c1


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

    validate(cloned_individual)
    return cloned_individual,

def random_matching_initializer(container, func, n):
    perm = np.random.permutation(range(0, n))

    num_pairs = int(np.ceil(n / 2))

    parts = [(x[0], x[1]) for x in np.split(perm, num_pairs)]

    return func(parts)

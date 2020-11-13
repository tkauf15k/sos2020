import numpy as np
import random
import DNAnon


def ordered_matching_best_match(ind1, ind2, toolbox):
    c1 = create_offspring(ind1, ind2, toolbox)
    c2 = create_offspring(ind2, ind1, toolbox)

    DNAnon.validate(c1)
    DNAnon.validate(c2)

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

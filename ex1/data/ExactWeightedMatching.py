import numpy as np
import networkx.algorithms
import networkx as nx


def ExactWeightedMatching(num_sequences: int, cost_matrix: np.ndarray):
    g = nx.Graph()

    for i in range(num_sequences):
        g.add_node(i)

    max_value = np.max(cost_matrix) + 1

    edges = [(a, b, {'weight': max_value - cost_matrix[a, b]}) for a in range(num_sequences) for b in range(num_sequences) if a < b]

    g.add_edges_from(edges)

    return list(networkx.algorithms.max_weight_matching(g))

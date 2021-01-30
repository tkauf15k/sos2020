import sklearn
import sklearn.metrics
import numpy as np

num_nodes = 1000
feature_dim = 4
radius = 0.75

np.random.seed(4)

input = np.random.uniform(low=-0.5, high=0.5, size=(num_nodes, feature_dim))

assignment_x = np.random.randint(0, 10, size=(num_nodes,))
assignment_y = np.random.randint(0, 10, size=(num_nodes,))

# print(assignment_x)
# print(assignment_y)

distances = sklearn.metrics.pairwise_distances(input)

tmp = distances < radius
np.fill_diagonal(tmp, False)

tmp2 = np.tril(tmp)
tmp3 = tmp2.astype(np.int)
index_matrix = np.array([list(range(0, num_nodes)), ] * num_nodes).transpose()
tmp4 = np.multiply(tmp3, index_matrix)

tmp5 = np.where(tmp4 > 0, tmp4, -1)

lines = set()

for i in range(0, num_nodes):
    my_coords = (assignment_x[i], assignment_y[i])

    my_partners = tmp5[:, i]
    my_partners_filtered = np.where(tmp5[:, i] > -1)

    if len(my_partners_filtered[0]) == 0:
        continue

    partner_x_coords = np.vectorize(lambda x: assignment_x[x])(
        my_partners_filtered)
    partner_y_coords = np.vectorize(lambda y: assignment_y[y])(
        my_partners_filtered)

    coords = np.concatenate([partner_x_coords, partner_y_coords],
                            axis=0).transpose()

    array_of_tuples = list(map(list, coords))

    neighbors = {tuple(val) for val in array_of_tuples}

    neighbors = {t for t in neighbors if t != my_coords}

    for n in neighbors:
        lines.add((my_coords, n))

    # print(f'{i} has {len(neighbors)} neighbors')

# print(lines)

print(f'number of lines: {len(lines)}')
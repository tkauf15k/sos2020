from typing import List, Tuple
import numpy as np

from lib import helperfunctions


def localsearch(solution: List[Tuple[int, int]], cost_matrix: np.ndarray, num_iterations=20, first_improvement = False):
    incumbent_solution = solution.copy()
    incumbent_solution_fitness = helperfunctions.fitness(incumbent_solution, cost_matrix)

    for i in range(0, num_iterations):

        delta, m1, m2 = improving_neighbor(incumbent_solution, cost_matrix, first_improvement)

        if delta is None or delta >= 0:
            break
        assert m1 != None and m2 != None

        apply(m1, m2, incumbent_solution)
        incumbent_solution_fitness += delta

        assert helperfunctions.fitness(incumbent_solution, cost_matrix) == incumbent_solution_fitness

    return (incumbent_solution, incumbent_solution_fitness)


def apply(m1: Tuple[int, Tuple[int, int]], m2: Tuple[int, Tuple[int, int]], solution: List[Tuple[int, int]]):
    m1_idx = m1[0]
    m1_matching = m1[1]

    m2_idx = m2[0]
    m2_matching = m2[1]

    solution[m1_idx] = m1_matching
    solution[m2_idx] = m2_matching

    helperfunctions.validate(solution)


def improving_neighbor(solution: List[Tuple[int, int]], cost_matrix: np.ndarray, first_improvement):
    best_improvement = 0
    m1 = None
    m2 = None

    for i, j in [(i, j) for i in range(0, len(solution)) for j in range(0, len(solution)) if i < j]:
        m11, m12 = solution[i]
        m21, m22 = solution[j]

        base_cost = cost_matrix[m11, m12] + cost_matrix[m21, m22]

        improved_cost = cost_matrix[m11, m21] + cost_matrix[m12, m22]
        if improved_cost < base_cost and improved_cost - base_cost < best_improvement:
            best_improvement = improved_cost - base_cost
            m1 = (i, (m11, m21))
            m2 = (j, (m12, m22))

        improved_cost = cost_matrix[m11, m22] + cost_matrix[m12, m21]
        if improved_cost < base_cost and improved_cost - base_cost < best_improvement:
            best_improvement = improved_cost - base_cost
            m1 = (i, (m11, m22))
            m2 = (j, (m12, m21))

        if best_improvement < 0 and first_improvement:
            break


    return best_improvement, m1, m2

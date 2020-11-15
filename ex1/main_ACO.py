# Contains the ACO solution

import sys
import timeit

import pants_modified as pants
import math
import random

from lib.helperfunctions import opt_gap, DNALA, fitness, preprocess
from data.ExactWeightedMatching import ExactWeightedMatching

num_sequences = int(sys.argv[1])
max_time = int(sys.argv[2])
input_path = "data/human_data_{}.fasta".format(num_sequences)

cost_matrix = preprocess(input_path, num_sequences)
exact_pairing = ExactWeightedMatching(num_sequences, cost_matrix)
optimum = fitness(exact_pairing, cost_matrix)
greedy_heuristic_pairing = DNALA(num_sequences, cost_matrix)
greedy_heuristic = fitness(greedy_heuristic_pairing, cost_matrix)

def arc_cost(a, b):
    return cost_matrix[a, b]


best_time = None
best_iteration = 0
start = 0

def improvement_callback(fitness, iteration):
    gap = opt_gap(fitness, optimum)
    global best_time
    global best_iteration
    global start
    best_iteration = iteration
    best_time = timeit.default_timer()
    print("Runtime = {}s, Fitness = {}, Delta={}, Gap={}".format(timeit.default_timer() - start, fitness, (fitness - optimum), gap))

if __name__ == "__main__":
    world = pants.World(range(0, num_sequences), arc_cost)

    solver = pants.Solver(ant_count=100, limit=max_time, optimum=optimum, improvement_callback=improvement_callback)

    start = timeit.default_timer()

    solution, iterations = solver.solve(world)

    best_fitness = solution.distance
    print("Tour: {0}".format(solution.tour))
    print()

    overall_runtime = timeit.default_timer() - start
    print("Terminated after: {} s / {} iterations, Best found after: {} s".format(overall_runtime, iterations, best_time))
    print(
        "Fitness = {}, Delta={}, Gap={}, Optimum={}, Greedy={}".format(best_fitness, (best_fitness - optimum), opt_gap(best_fitness, optimum),
                                                                    optimum, greedy_heuristic))

    # print(solution.path)

import sys
import timeit

import pants
import math
import random

import DNALA
import DNAnon
import evaluation

num_sequences = 40  # 20 50 100
input_path = "data/human_data_{}.fasta".format(num_sequences)

cost_matrix = DNAnon.preprocess(input_path, num_sequences)
exact_pairing = DNAnon.ExactWeightedMatching(num_sequences, cost_matrix)
optimum = DNAnon.fitness(exact_pairing, cost_matrix)
greedy_heuristic_pairing = DNALA.DNALA(num_sequences, cost_matrix)
greedy_heuristic = DNAnon.fitness(greedy_heuristic_pairing, cost_matrix)


def arc_cost(a, b):
    return cost_matrix[a, b]


best_time = None
best_iteration = 0


def improvement_callback(fitness, iteration):
    gap = evaluation.opt_gap(fitness, optimum)
    global best_time
    global best_iteration
    best_iteration = iteration
    best_time = timeit.default_timer()
    print("Fitness = {}, Delta={}, Gap={}".format(fitness, (fitness - optimum), gap))


world = pants.World(range(0, num_sequences), arc_cost)

solver = pants.Solver(ant_count=100, limit=10000, optimum=optimum, improvement_callback=improvement_callback)

start = timeit.default_timer()

solution, iterations = solver.solve(world)

best_fitness = solution.distance
print(solution.tour)

overall_runtime = timeit.default_timer() - start
print("Terminated after: {} s / {} iterations, Best found after: {} s".format(overall_runtime, iterations, best_time))
print(
    "Fitness = {}, Delta={}, Gap={}, Optimum={}, Greedy={}".format(best_fitness, (best_fitness - optimum), evaluation.opt_gap(best_fitness, optimum),
                                                                   optimum, greedy_heuristic))

# print(solution.path)

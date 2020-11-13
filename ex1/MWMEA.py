import random
import sys
import timeit

import numpy as np
from deap import creator, base, tools, algorithms
import DNAnon
import DNALA
import DNAnonMutations
import DNAnonRecombinations
import DNAnonInitializer
import evaluation

num_sequences = 70  # 20 50 100
input_path = "data/human_data_{}.fasta".format(num_sequences)

if __name__ == "__main__":
    cost_matrix = DNAnon.preprocess(input_path, num_sequences)
    greedy_heuristic_pairing = DNALA.DNALA(num_sequences, cost_matrix)
    exact_pairing = DNAnon.ExactWeightedMatching(num_sequences, cost_matrix)

    optimum = DNAnon.fitness(exact_pairing, cost_matrix)
    greedy_heuristic = DNAnon.fitness(greedy_heuristic_pairing, cost_matrix)

    creator.create("FitnessFunction", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessFunction)

    toolbox = base.Toolbox()
    toolbox.register("individual", DNAnonInitializer.random_matching_initializer, None, creator.Individual, n=num_sequences)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    mutation_rate = 1 / (num_sequences / 2)

    toolbox.register("evaluate", DNAnon.fitness_tuple, cost_matrix=cost_matrix)
    toolbox.register("mate", DNAnonRecombinations.arbitrary_matched_subgraph, toolbox=toolbox)
    toolbox.register("mutate", DNAnonMutations.two_opt_mutation, indpb=mutation_rate, toolbox=toolbox, cost_matrix=cost_matrix,
                     edge_selector=DNAnonMutations.select_best_edge)
    toolbox.register("select", tools.selTournament, tournsize=3)

    population = toolbox.population(n=300)

    best_solution = None
    best_fitness = sys.maxsize
    best_time = 0

    start = timeit.default_timer()
    NGEN = 500
    for gen in range(NGEN):
        offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.1)
        fits = toolbox.map(toolbox.evaluate, offspring)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit
        population = toolbox.select(offspring, k=len(population))
        best = tools.selBest(population, k=1)[0]

        if best.fitness.values[0] < best_fitness:
            best_fitness = best.fitness.values[0]
            best_solution = toolbox.clone(best)
            gap = evaluation.opt_gap(best_fitness, optimum)
            best_time = timeit.default_timer() - start
            print("Fitness = {}, Delta={}, Gap={}".format(best_fitness, (best_fitness - optimum), gap))
            if gap <= 0:
                break

    overall_runtime = timeit.default_timer() - start
    print("Terminated after: {} s / {} iterations, Best found after: {} s".format(overall_runtime, gen, best_time))
    print("Fitness = {}, Delta={}, Gap={}, Optimum={}, Greedy={}".format(best_fitness, (best_fitness - optimum), evaluation.opt_gap(best_fitness, optimum),
                                                                         optimum,
                                                                         greedy_heuristic))

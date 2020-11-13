import random
import sys
import timeit

import numpy as np
from deap import creator, base, tools, algorithms
import DNAnon
import DNALA

num_sequences = 50  # 20 50 100
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


    def mwm_individual_initializer(container, func, n):
        perm = np.random.permutation(range(0, n))

        num_pairs = int(np.ceil(n / 2))

        parts = [(x[0], x[1]) for x in np.split(perm, num_pairs)]

        return func(parts)


    toolbox.register("individual", mwm_individual_initializer, None, creator.Individual, n=num_sequences)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)


    def evaluate_individual(individual):
        return (DNAnon.fitness(individual, cost_matrix),)


    def validate(individual):
        nodes = [i for a in individual for i in a]
        assert len(set(nodes)) == len(nodes)


    def crossover(ind1, ind2):

        c1 = create_offspring(ind1, ind2)
        c2 = create_offspring(ind2, ind1)

        validate(c1)
        validate(c2)

        return c1, c2


    def create_offspring(ind1, ind2):
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


    def mutation(individual, indpb):

        probabilities = np.random.random(size=len(individual))

        condition = probabilities <= 2 * indpb

        pairs = [idx for idx, elem in enumerate(condition) if elem == True]
        random.shuffle(pairs)

        if pairs is None or len(pairs) < 2:  # for small instances, we sometimes cannot mutate properly
            return toolbox.clone(individual),

        truncated_length = int(len(pairs) / 2) * 2

        cloned_individual = toolbox.clone(individual)
        for i, j in [(i, (i + 1)) for i in range(truncated_length - 1)]:
            p1 = cloned_individual[pairs[i]]
            p2 = cloned_individual[pairs[j]]

            indices = np.random.randint(0, 2, 2)

            al1 = (p1[indices[0]], p2[indices[1]])
            al2 = (p1[1 - indices[0]], p2[1 - indices[1]])

            cloned_individual[pairs[i]] = al1
            cloned_individual[pairs[j]] = al2

        validate(cloned_individual)
        return cloned_individual,


    def opt_gap(fitness):
        return ((1.0 * fitness - optimum) / optimum) * 100  # in percent...


    mutation_rate = 1 / (num_sequences / 2)

    toolbox.register("evaluate", evaluate_individual)
    toolbox.register("mate", crossover)
    toolbox.register("mutate", mutation, indpb=mutation_rate)
    toolbox.register("select", tools.selTournament, tournsize=3)

    population = toolbox.population(n=300)

    best_solution = None
    best_fitness = sys.maxsize
    best_time = 0

    start = timeit.default_timer()
    NGEN = 1000
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
            gap = opt_gap(best_fitness)
            best_time = timeit.default_timer() - start
            print("Fitness = {}, Delta={}, Gap={}".format(best_fitness, (best_fitness - optimum), gap))
            if gap <= 0:
                break

    overall_runtime = timeit.default_timer() - start
    print("Terminated after: {} s, Best found after: {} s".format(overall_runtime, best_time))
    print("Fitness = {}, Delta={}, Gap={}, Optimum={}, Greedy={}".format(best_fitness, (best_fitness - optimum), opt_gap(best_fitness), optimum,
                                                                         greedy_heuristic))

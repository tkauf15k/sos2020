import timeit
from deap import creator, base, tools, algorithms
from lib import localsearch
from lib.helperfunctions import preprocess, DNALA, fitness, fitness_tuple, opt_gap
from lib.GA_operations import *
from data.ExactWeightedMatching import ExactWeightedMatching

num_sequences = int(sys.argv[1])
max_time = int(sys.argv[2])
assert max_time > 0
input_path = "data/human_data_{}.fasta".format(num_sequences)

if __name__ == "__main__":
    cost_matrix = preprocess(input_path, num_sequences)
    greedy_heuristic_pairing = DNALA(num_sequences, cost_matrix)
    exact_pairing = ExactWeightedMatching(num_sequences, cost_matrix)
    optimum = fitness(exact_pairing, cost_matrix)
    greedy_heuristic = fitness(greedy_heuristic_pairing, cost_matrix)

    creator.create("FitnessFunction", base.Fitness, weights=(-1.0,))
    creator.create("Individual", list, fitness=creator.FitnessFunction)

    toolbox = base.Toolbox()
    toolbox.register("individual", random_matching_initializer, None, creator.Individual, n=num_sequences)
    toolbox.register("population", tools.initRepeat, list, toolbox.individual)

    mutation_rate = 1 / (num_sequences / 2)

    toolbox.register("evaluate", fitness_tuple, cost_matrix=cost_matrix)
    toolbox.register("mate", arbitrary_matched_subgraph, toolbox=toolbox)
    toolbox.register("mutate", two_opt_mutation, indpb=mutation_rate, toolbox=toolbox, cost_matrix=cost_matrix,
                     edge_selector=select_best_edge)
    toolbox.register("select", tools.selTournament, tournsize=3)

    population = toolbox.population(n=300)

    best_solution = None
    best_fitness = sys.maxsize
    best_time = 0

    start = timeit.default_timer()

    gen = 0
    while True:
        offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.1)
        fits = toolbox.map(toolbox.evaluate, offspring)
        for fit, ind in zip(fits, offspring):
            ind.fitness.values = fit
        population = toolbox.select(offspring, k=len(population))

        # ls improvement
        # population = [ls_convert(localsearch(f, cost_matrix, num_iterations=100, first_improvement=True)) for f in population]

        best = tools.selBest(population, k=1)[0]

        runtime = timeit.default_timer() - start
        if best.fitness.values[0] < best_fitness:
            best_fitness = best.fitness.values[0]
            best_solution = toolbox.clone(best)
            gap = opt_gap(best_fitness, optimum)
            best_time = timeit.default_timer() - start
            print("Runtime = {}s, Fitness = {}, Delta={}, Gap={}".format(runtime, best_fitness, (best_fitness - optimum), gap))
            if gap <= 0:
                break
        if runtime > max_time:
            break
        gen += 1

    overall_runtime = timeit.default_timer() - start
    print("Terminated after: {} s / {} iterations, Best found after: {} s".format(overall_runtime, gen, best_time))
    print("Fitness = {}, Delta={}, Gap={}, Optimum={}, Greedy={}".format(best_fitness, (best_fitness - optimum), opt_gap(best_fitness, optimum),
                                                                         optimum,
                                                                         greedy_heuristic))

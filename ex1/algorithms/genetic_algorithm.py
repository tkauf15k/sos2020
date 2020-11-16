from deap import base, tools, algorithms
from algorithms.algorithm_base import algorithm_base
from lib import localsearch
from lib.GA_operations import *
from lib.helperfunctions import fitness_tuple


class genetic_algorithm(algorithm_base):

    def __init__(self, population_size, improving_mutation, mutation_rate_factor, selection_strategy, selection_k, localsearch):
        super().__init__()
        self.population_size = population_size
        self.improving_mutation = improving_mutation
        self.mutation_rate_factor = mutation_rate_factor
        self.selection_strategy = selection_strategy
        self.selection_k = selection_k
        self.localsearch = localsearch

    def register_selection_function(self, toolbox):
        if self.selection_strategy == 'Tournament':
            toolbox.register("select", tools.selTournament, tournsize=self.selection_k)
        elif self.selection_strategy == 'SUS':
            toolbox.register("select", tools.selStochasticUniversalSampling, k=self.selection_k)
        elif self.selection_strategy == 'Roulette':
            toolbox.register("select", tools.selRoulette, k=self.selection_k)
        elif self.selection_strategy == 'Best':
            toolbox.register("select", tools.selBest, k=self.selection_k)
        else:
            assert False

    def solve(self, num_sequences, time_limit):
        cost_matrix, optimum = self.prepare(num_sequences)
        creator.create("FitnessFunction", base.Fitness, weights=(-1.0,))
        creator.create("Individual", list, fitness=creator.FitnessFunction)

        toolbox = base.Toolbox()
        toolbox.register("individual", random_matching_initializer, None, creator.Individual, n=num_sequences)
        toolbox.register("population", tools.initRepeat, list, toolbox.individual)

        mutation_rate = (1 / (num_sequences / 2)) * self.mutation_rate_factor

        toolbox.register("evaluate", fitness_tuple, cost_matrix=cost_matrix)
        toolbox.register("mate", arbitrary_matched_subgraph, toolbox=toolbox)
        self.register_selection_function(toolbox)
        toolbox.register("mutate", two_opt_mutation, indpb=mutation_rate, toolbox=toolbox, cost_matrix=cost_matrix,
                         edge_selector=select_best_edge if self.improving_mutation else select_random_edge)

        population = toolbox.population(n=self.population_size)

        best_fitness = sys.maxsize
        self.start_run()

        gen = 0
        while True:
            offspring = algorithms.varAnd(population, toolbox, cxpb=0.5, mutpb=0.1)
            fits = toolbox.map(toolbox.evaluate, offspring)
            for fit, ind in zip(fits, offspring):
                ind.fitness.values = fit
            population = toolbox.select(offspring, k=len(population))

            if self.localsearch:
                population = [ls_convert(localsearch(f, cost_matrix, num_iterations=100, first_improvement=True)) for f in population]

            best = tools.selBest(population, k=1)[0]

            if best.fitness.values[0] < best_fitness:
                best_fitness = best.fitness.values[0]
                best_solution = toolbox.clone(best)
                if self.report_improvement(best_fitness, gen, best_solution):
                    break

            if self.terminated(time_limit):
                break
            gen += 1

        self.stop_run(gen)

        # print("Terminated after: {} s / {} iterations, Best found after: {} s".format(overall_runtime, gen, best_time))
        # print("Fitness = {}, Delta={}, Gap={}, Optimum={}, Greedy={}".format(best_fitness, (best_fitness - optimum), opt_gap(best_fitness, optimum),
        #                                                                      optimum,
        #                                                                      greedy_heuristic))

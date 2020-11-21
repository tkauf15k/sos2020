import pickle
import timeit
from abc import abstractmethod

from data.ExactWeightedMatching import ExactWeightedMatching
from lib.evaluation import opt_gap
from lib.helperfunctions import preprocess, fitness


class algorithm_base:

    def __init__(self):
        self.optimum = 0
        self.trace = []
        self.start_time = 0
        self.stop_time = 0
        self.best_solution = None
        self.overall_iterations = 0
        self.best_iteration = 0

    @abstractmethod
    def solve(self, num_sequences, time_limit):
        pass

    def report(self, output):
        with open(output, 'w') as out_file:
            out_file.write("optimum={}\n".format(self.optimum))
            out_file.write("execution_time={},start_time={},stop_time={}\n".format(self.stop_time - self.start_time, self.start_time, self.stop_time))
            out_file.write("best_iteration={},overall_iterations={}\n".format(self.best_iteration, self.overall_iterations))
            out_file.write(str(self.best_solution) + "\n")

            for t in self.trace:
                (fitness, iteration, gap, delta, timestamp) = t
                out_file.write('fitness={},iteration={},gap={},delta={},timestamp={}\n'.format(fitness, iteration, gap, delta, timestamp))
        # if output.endswith('.txt'):
        #     pickle_path = output.replace('.txt', '.pickle')
        #     pickle.dump(self, open(pickle_path, "wb"))

    def prepare(self, num_sequences):
        input_path = "data/human_data_{}.fasta".format(num_sequences)
        cost_matrix = preprocess(input_path, num_sequences)

        exact_pairing = ExactWeightedMatching(num_sequences, cost_matrix)
        optimum = fitness(exact_pairing, cost_matrix)
        self.optimum = optimum

        return cost_matrix, optimum

    def report_improvement(self, fitness, iteration, solution):
        gap = opt_gap(fitness, self.optimum)
        delta = fitness - self.optimum
        timestamp = timeit.default_timer()
        self.trace.append((fitness, iteration, gap, delta, timestamp))
        self.best_solution = solution
        self.best_iteration = iteration

        print("fitness={}, iteration={}, gap={}, delta={}, execution-time={}".format(fitness, iteration, gap, delta, timestamp))
        return gap <= 0

    def start_run(self):
        self.start_time = timeit.default_timer()

    def stop_run(self, iterations):
        self.stop_time = timeit.default_timer()
        self.overall_iterations = iterations

    def terminated(self, time_limit):
        if timeit.default_timer() - self.start_time >= time_limit:
            print('terminated!')
            return True
        return False

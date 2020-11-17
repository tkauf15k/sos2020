from algorithms.algorithm_base import algorithm_base
from lib.helperfunctions import DNALA, fitness
from lib.localsearch import localsearch


class dnala(algorithm_base):

    def __init__(self, localsearch):
        super().__init__()
        self.localsearch = localsearch

    def solve(self, num_sequences, time_limit):
        cost_matrix, optimum = self.prepare(num_sequences)

        self.start_run()
        iterations = 1
        pairing_dnala = DNALA(num_sequences, cost_matrix)
        cost = fitness(pairing_dnala, cost_matrix)
        self.report_improvement(cost, 0, pairing_dnala)

        if self.localsearch:
            iterations += 1
            pairing_dnala, cost = localsearch(pairing_dnala, cost_matrix, num_iterations=50, first_improvement=True, termination=(lambda: self.terminated(time_limit)))
            self.report_improvement(cost, 1, pairing_dnala)

        self.stop_run(iterations)





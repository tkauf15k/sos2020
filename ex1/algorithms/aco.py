import pants

from algorithms.algorithm_base import algorithm_base


class aco(algorithm_base):

    def __init__(self, populationsize, alpha, beta, rho, localsearch):
        super().__init__()
        self.populationsize = populationsize
        self.alpha = alpha
        self.beta = beta
        self.rho = rho
        self.localsearch = localsearch

    def solve(self, num_sequences, time_limit):
        cost_matrix, optimum = self.prepare(num_sequences)

        def arc_cost(a, b):
            return cost_matrix[a, b]

        def improvement_callback(fitness, iteration, solution):
            matching = [(solution.tour[i], solution.tour[i + 1]) for i in range(0, len(solution.tour) - 1, 2)]
            self.report_improvement(fitness, iteration, matching)

        def local_search_callback(tour):
            matching = [(tour[i], tour[i + 1]) for i in range(0, len(tour) - 1, 2)]
            assert len(matching) == num_sequences / 2
            from lib.localsearch import localsearch
            improved_matching, fitness = localsearch(matching, cost_matrix, num_iterations=50, first_improvement=True,
                                                     termination=(lambda: self.terminated(time_limit)))
            improved_tour = [j for i in improved_matching for j in i]
            assert len(improved_tour) == len(tour)
            assert len(set(improved_tour)) == len(improved_tour)
            return improved_tour, fitness

        world = pants.World(range(0, num_sequences), arc_cost)

        lscb = local_search_callback if self.localsearch else None
        solver = pants.Solver(ant_count=self.populationsize, alpha=self.alpha, beta=self.beta, rho=self.rho,
                              limit=time_limit, optimum=optimum, improvement_callback=improvement_callback,
                              local_search_callback=lscb)

        self.start_run()

        solution, iterations = solver.solve(world)

        self.stop_run(iterations)

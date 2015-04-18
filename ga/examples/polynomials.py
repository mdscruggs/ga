try:
    import pylab as py
    py.style.use('ggplot')
except ImportError:
    py = None

from ..algorithms import BaseGeneticAlgorithm
from ..chromosomes import Chromosome
from ..translators import BinaryFloatTranslator


class PolyModelGA(BaseGeneticAlgorithm):
    """
    A GA that attempts to find optimal coefficients of a polynomial.
    """
    def __init__(self, coefficients, num_x, significand_length, *args, **kwargs):
        """
        Construct a new ``PolyModelGA`` instance.

        coefficients:  list of polynomial coefficients
        num_x:  number of x-values to compute, starting at 1
        significand_length:  length of significand in genes
        *args, **kwargs forwarded to ``BaseGeneticAlgorithm`` constructor
        """
        super().__init__(*args, **kwargs)
        self.translator = BinaryFloatTranslator(significand_length)

        self.coefficients = coefficients
        self.num_x = num_x
        self.expected_values = self.compute_y(coefficients, num_x)

        self.fitness_cache = {}

    def compute_y(self, coefficients, num_x):
        """ Return calculated y-values for the domain of x-values in [1, num_x]. """
        y_vals = []

        for x in range(1, num_x + 1):
            y = sum([c * x ** i for i, c in enumerate(coefficients[::-1])])
            y_vals.append(y)

        return y_vals

    def compute_err(self, solution_y, coefficients):
        """
        Return an error value by finding the absolute difference for each
        element in a list of solution-generated y-values versus expected values.

        Compounds error by 50% for each negative coefficient in the solution.

        solution_y:  list of y-values produced by a solution
        coefficients:  list of polynomial coefficients represented by the solution

        return:  error value
        """
        error = 0
        for modeled, expected in zip(solution_y, self.expected_values):
            error += abs(modeled - expected)

        if any([c < 0 for c in coefficients]):
            error *= 1.5

        return error

    def eval_fitness(self, chromosome):
        """
        Evaluate the polynomial equation using coefficients represented by a
        solution/chromosome, returning its error as the solution's fitness.

        return:  fitness value
        """
        if chromosome.dna in self.fitness_cache:
            return self.fitness_cache[chromosome.dna]

        coefficients = self.translator.translate_chromosome(chromosome)
        solution_y = self.compute_y(coefficients, self.num_x)
        fitness = -1 * self.compute_err(solution_y, coefficients)

        self.fitness_cache[chromosome.dna] = fitness
        return fitness


def run(coefficients=(0.001, 0.01, 0.1, 1), num_x=10, generations=5000, plot=True):
    # fit a polynomial equation to expected values
    
    poly_str = ''
    for i, c in enumerate(coefficients[::-1]):
        if i > 1:
            poly_str = '{0:+} x^{1} '.format(c, i) + poly_str
        elif i == 1:
            poly_str = '{0:+} x '.format(c) + poly_str
        else:
            poly_str = '{0:+}'.format(c) + poly_str
            
    poly_str = '$' + poly_str + '$'
    
    # gene lengths are a little complex:
    #   1 bit for significand sign
    #   8 bits for significand body
    #   1 bit for exponent sign
    #   2 bits for exponent body (represent 0, 1, 2)
    significand_length = 8
    gene_length = (1 + significand_length + 1 + 2,) * len(coefficients)  # 1 gene per polynomial coefficient
    chromosomes = Chromosome.create_random(gene_length, n=20)
    
    poly_ga = PolyModelGA(coefficients, num_x, significand_length, chromosomes, abs_fit_weight=1, rel_fit_weight=0)
    
    p_mutate = 0.15
    p_cross = 0.50
    best = poly_ga.run(generations, p_mutate, p_cross, elitist=True)
    
    best_coeff = poly_ga.translator.translate_chromosome(best)
    best_y = poly_ga.compute_y(best_coeff, num_x)

    print("run took", poly_ga.run_time_s, "seconds")
    print("best solution coefficients =", best_coeff, "error =", poly_ga.compute_err(best_y, best_coeff))

    if plot:
        if py:
            # plot a curve for every solution that caused an upset
            # older solutions have higher transparency
            x = range(1, num_x + 1)

            for new_best_gen in poly_ga.new_fittest_generations:
                gen_fittest = poly_ga.generation_fittest[new_best_gen]
                sol_coefficients = poly_ga.translator.translate_chromosome(gen_fittest)
                modeled_y = poly_ga.compute_y(sol_coefficients, num_x)
                alpha = min(0.5, max(0.2, new_best_gen/max(poly_ga.new_fittest_generations)))

                if new_best_gen == poly_ga.new_fittest_generations[-1]:
                    py.plot(x, modeled_y, color='green', label='best', marker='x', linewidth=4, markersize=12, alpha=alpha)
                else:
                    py.plot(x, modeled_y, color='red', linestyle='--', alpha=alpha)

            py.plot(x, poly_ga.expected_values, color='blue', marker='o', label='expected')
            py.ylim(min(poly_ga.expected_values) - min(poly_ga.expected_values), max(poly_ga.expected_values) * 2)
            py.legend(loc='best')
            py.title('PolyModelGA ({} gen.)\n{}'.format(generations, poly_str))

            py.show()
        else:
            print("Did not plot example results because matplotlib not installed")
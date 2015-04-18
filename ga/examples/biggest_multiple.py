import functools

try:
    import pylab as py
    py.style.use('ggplot')
except ImportError:
    py = None

from ..algorithms import BaseGeneticAlgorithm
from ..chromosomes import Chromosome
from ..translators import BinaryIntTranslator


class BiggestMultipleGA(BaseGeneticAlgorithm):
    """
    A GA that tries to find the largest multiple of a set of factors,
    within the encodable space of the given chromosomes.

    Forces use of a ``translators.BinaryIntTranslator``.
    """
    def __init__(self, factors, *args, **kwargs):
        """
        Construct a new ``BiggestMultipleGA`` instance.

        Only allows chromosomes with 1 gene each.

        factors:  sequence of integer factors >= 1 that solutions must have.
        """
        kwargs['translator'] = BinaryIntTranslator()
        super().__init__(*args, **kwargs)

        # only 1-gene chromosomes allowed
        for c in self.chromosomes:
            assert len(c.genes) == 1

        assert factors and all(isinstance(i, int) and i >= 1 for i in factors)
        self.factors = factors

        self.max_encoded_val = int('1' * self.chromosomes[0].length, base=2)
        self.best_possible_fit = len(self.factors)
        self.found_best = False

    def eval_fitness(self, chromosome):
        """
        Convert a 1-gene chromosome into an integer and calculate its fitness
        by checking it against each required factor.

        return:  fitness value
        """
        score = 0

        number = self.translator.translate_gene(chromosome.genes[0])

        for factor in self.factors:
            if number % factor == 0:
                score += 1
            else:
                score -= 1

        # check for optimal solution
        if score == len(self.factors):
            min_product = functools.reduce(lambda a,b: a * b, self.factors)

            if number + min_product > self.max_encoded_val:
                print("Found best solution:", number)
                self.found_best = True

        # scale number of factors achieved/missed by ratio of
        # the solution number to the maximum possible integer
        # represented by binary strings of the given length
        return score * number / self.max_encoded_val

    def should_terminate(self, overall_fittest):
        return self.found_best


def run(factors=(2, 3, 7, 11), gene_length=16, generations=2**32, plot=True):
    # find largest encodable integer that has all factors
    chromosomes = Chromosome.create_random(gene_length, n=10)

    bm_ga = BiggestMultipleGA(factors, chromosomes)
    
    # run until we (hopefully) find the largest encodable 
    # value that is a product of all the factors
    p_mutate = 0.15
    p_cross = 0.25
    best = bm_ga.run(generations, p_mutate, p_cross, elitist=True)
    best_num = bm_ga.translator.translate_gene(best.genes[0])
    
    print("run took", bm_ga.run_time_s, "seconds")
    print("Best solution is:", best, best_num)
    
    for factor in factors:
        assert best_num % factor == 0

    if plot:
        if py:
            # plot fitness progression
            py.plot([v for k, v in sorted(bm_ga.overall_fittest_fit.items())], label='run')

            py.xlabel('generation')
            py.ylabel('fitness')
            py.legend(loc='best')

            py.show()
        else:
            print("Did not plot example results because matplotlib not installed")

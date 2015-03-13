try:
    import pylab as py
    py.style.use('ggplot')
except ImportError:
    py = None

from ..algorithms import BiggestMultipleGA
from .. import util


def run(factors=(2, 3, 7, 11), gene_length=16, generations=2**32):
    # find largest encodable integer that has all factors
    chromosomes = util.random_chromosomes(10, gene_length)

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

    if py:
        # plot fitness progression
        py.plot([v for k, v in sorted(bm_ga.overall_fittest_fit.items())], label='run')

        py.xlabel('generation')
        py.ylabel('fitness')
        py.legend(loc='best')

        py.show()
    else:
        print("Did not plot example results because matplotlib not installed")

import pylab as py
py.style.use('ggplot')

from ..algorithms import BiggestMultipleGA
from .. import util

def run_example():
    # find largest encodable integer that has all factors
    gene_length = 16
    chromosomes = util.random_chromosomes(10, gene_length)

    factors = (2, 3, 7, 11)
    bm_ga = BiggestMultipleGA(factors, chromosomes)
    
    # run until we (hopefully) find the largest encodable 
    # value that is a product of all the factors
    generations = 2**32 
    p_mutate = 0.15
    p_cross = 0.25
    best = bm_ga.run(generations, p_mutate, p_cross, elitist=True)
    best_num = bm_ga.translator.translate_gene(best.genes[0])
    
    print("run took", bm_ga.run_time_s, "seconds")
    print("Best solution is:", best, best_num)
    
    for factor in factors:
        assert best_num % factor == 0
    
    # plot fitness progression
    py.plot(bm_ga.overall_fittest_fit, label='run')
    
    py.xlabel('generation')
    py.ylabel('fitness')
    py.legend(loc='best')
    
    py.show()

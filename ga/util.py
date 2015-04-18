import random


def compute_fitness_cdf(chromosomes, ga):
    """
    Return a list of fitness-weighted cumulative probabilities for a set of chromosomes.
    
    chromosomes:  chromosomes to use for fitness-based calculations
    ga:  ``algorithms.BaseGeneticAlgorithm`` used to obtain fitness values using its ``eval_fitness`` method
    
    return:  list of fitness-weighted cumulative probabilities in [0, 1]
    """
    ga.sort(chromosomes)
    
    fitness = [ga.eval_fitness(c) for c in chromosomes]
    min_fit = min(fitness)
    fit_range = max(fitness) - min_fit
    
    if fit_range == 0:
        # all chromosomes have equal chance of being chosen
        n = len(chromosomes)
        return [i / n for i in range(1, n + 1)]
        
    return [(fit - min_fit) / fit_range for fit in fitness]
    
    
def weighted_choice(seq, cdf):
    """
    Select a random element from a sequence, given cumulative probabilities of selection.
    
    See ``compute_fitness_cdf`` function for obtaining cumulative probabilities.
    
    seq:  sequence to select from
    cdf:  sequence with 1 cumulative probability value in [0, 1] for each element in ``seq``
    
    return:  randomly selected element
    """
    assert len(seq) == len(cdf)
    rand = random.random()
    
    for i, e in enumerate(seq):
        cp = cdf[i]
        assert 0 <= cp <= 1
        
        if rand < cp:
            return e
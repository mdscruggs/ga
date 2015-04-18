import abc
import random
import time

from .chromosomes import Chromosome
from .util import weighted_choice, compute_fitness_cdf


class BaseGeneticAlgorithm(abc.ABC):
    """
    Base genetic algorithm.
    
    Subclasses must override the ``eval_fitness`` method.
    """
    def __init__(self, chromosomes, translator=None, abs_fit_weight=0.25, rel_fit_weight=0.75):
        """
        Construct a new ``BaseGeneticAlgorithm`` instance.
        
        chromosomes:  collection of initial chromosomes
        
        translator (default=None):  ``translators.BaseTranslator`` instance that may be needed 
                                    by ``eval_fitness`` method
        
        abs_fit_weight (default=0.25, must be in [0, 1]):  
          The fraction of survival probabilities apportioned to a run's all-time fitness range.
          
          Increasing this weight places greater emphasis on a solution's fitness versus the
          overall observed fitness range within a run.
          
          Think of this as "environmental" pressure where bad solutions cannot survive even if
          they are much better than all other solutions (e.g. a powerful shark on dry land).
          
        rel_fit_weight (default=0.75, must be in [0, 1]):
          The fraction of survival probabilities apportioned to relative fitness each generation.
          
          Increasing this weight places greater emphasis on a solution's fitness versus the
          fitness range within the current generation.
          
          Think of this as "competitive" pressure where even bad solutions which are relatively
          better get an relative advantage.
          
        Asserts that (abs_fit_weight + rel_fit_weight) equals 1.
        """
        assert all(isinstance(c, Chromosome) for c in chromosomes)
        assert 0 <= abs_fit_weight <= 1
        assert 0 <= rel_fit_weight <= 1
        assert abs_fit_weight + rel_fit_weight == 1
        
        self.chromosomes = chromosomes
        self.translator = translator
        self.abs_fit_weight = abs_fit_weight
        self.rel_fit_weight = rel_fit_weight
        
        self.orig_pop_size = len(self.chromosomes)
        self.min_fit_ever = None
        self.max_fit_ever = None
        
        # run results
        self.generation_fittest = {}
        self.generation_fittest_fit = {}
        self.overall_fittest_fit = {}
        self.new_fittest_generations = []
        self.run_time_s = None

    @abc.abstractmethod
    def eval_fitness(self, chromosome):
        """ Return a fitness score for a chromosome. """
        pass
            
    def get_fittest(self):
        """ Get the chromosome with the highest fitness score. """
        return max(self.chromosomes, key=self.eval_fitness)
        
    def get_weakest(self):
        """ Get the chromosome with the lowest fitness score. """
        return min(self.chromosomes, key=self.eval_fitness)
        
    def sort(self, chromosomes):
        """ 
        Sort a list of chromosomes into ascending order based on fitness score.
        
        chromosomes:  list of chromosomes to sort in-place
        """
        chromosomes.sort(key=self.eval_fitness)
        
    def compete(self, chromosomes):
        """
        Simulate competition/survival of the fittest.
        
        The fitness of each chromosome is used to calculate a survival probability
        based on how it compares to the overall fitness range of the run and the
        fitness range of the current generation. The ``abs_fit_weight`` and 
        ``rel_fit_weight`` attributes determine the degree to which overall and
        current fitness ranges affect the final survival probability.
        
        return:  list of surviving chromosomes
        """
        # update overall fitness for this run
        self.sort(chromosomes)
        min_fit = self.eval_fitness(chromosomes[0])
        max_fit = self.eval_fitness(chromosomes[-1])
        
        if min_fit < self.min_fit_ever:
            self.min_fit_ever = min_fit
        if max_fit > self.max_fit_ever:
            self.max_fit_ever = max_fit
        
        overall_fit_range = self.max_fit_ever - self.min_fit_ever  # "absolute" fitness range
        current_fit_range = max_fit - min_fit                      # "relative" fitness range
        
        # choose survivors based on relative fitness within overall fitness range
        survivors = []
        for chromosome in chromosomes:
            fit = self.eval_fitness(chromosome)
            p_survival_absolute = (fit - self.min_fit_ever) / overall_fit_range if overall_fit_range != 0 else 1
            p_survival_relative = (fit - min_fit) / current_fit_range if current_fit_range != 0 else 1
            
            # compute weighted average survival probability
            # a portion accounts for absolute overall fitness for all chromosomes ever encountered (environment-driven)
            # the other portion accounts for relative fitness within the current population (competition-driven)
            p_survival = p_survival_absolute * self.abs_fit_weight + p_survival_relative * self.rel_fit_weight

            if random.random() < p_survival:
                survivors.append(chromosome)

        if not survivors:
            # rarely, nothing survives -- allow everyone to live
            return chromosomes
                
        return survivors
        
    def reproduce(self, survivors, p_crossover, target_size=None):
        """
        Reproduces the population from a pool of surviving chromosomes
        until a target population size is met. Offspring are created
        by selecting a survivor. Survivors with higher fitness have a
        greater chance to be selected for reproduction.
        
        Genetic crossover events may occur for each offspring created.
        Crossover mates are randomly selected from the pool of survivors.
        Crossover points are randomly selected from the length of the crossed chromosomes.
        If crossover does not occur, an offspring is an exact copy of the selected survivor.
        Crossover only affects the DNA of the offspring, not the survivors/parents.
        
        survivors:  pool of parent chromosomes to reproduce from
        p_crossover:  probability in [0, 1] that a crossover event will
                      occur for each offspring
        target_size (default=original population size): target population size
        
        return:  list of survivors plus any offspring
        """
        assert 0 <= p_crossover <= 1
        
        if not target_size:
            target_size = self.orig_pop_size
            
        num_survivors = len(survivors)
    
        # compute reproduction cumulative probabilities
        # weakest member gets p=0 but can be crossed-over with
        cdf = compute_fitness_cdf(survivors, self)
        
        offspring = []
        while num_survivors + len(offspring) < target_size:
            # pick a survivor to reproduce
            c1 = weighted_choice(survivors, cdf).copy()
            
            # crossover
            if random.random() < p_crossover:
                # randomly pick a crossover mate from survivors
                # same chromosome can be c1 and c2
                c2 = random.choice(survivors).copy()
                crosspoint = random.randrange(0, c1.length)
                c1.crossover(c2, crosspoint)
                
            offspring.append(c1)
            
        return survivors + offspring
        
    def mutate(self, chromosomes, p_mutate):
        """ 
        Call every chromosome's ``mutate`` method. 
        
        p_mutate:  probability of mutation in [0, 1]
        """
        assert 0 <= p_mutate <= 1
        
        for chromosome in chromosomes:
            chromosome.mutate(p_mutate)
            
    def refresh(self, chromosomes, p_mutate=0.5):
        """
        Refresh chromosomes with a high mutation rate.
        
        chromosomes:  chromosomes to mutate
        p_mutate (default=0.5):  mutation rate in [0, 1]
        """
        self.mutate(chromosomes, p_mutate)
        
    def run(self, generations, p_mutate, p_crossover, elitist=True, refresh_after=None, quit_after=None):
        """
        Run a standard genetic algorithm simulation for a set number
        of generations (iterations), each consisting of the following
        ordered steps:
          
          1. competition/survival of the fittest (``compete`` method)
          2. reproduction (``reproduce`` method)
          3. mutation (``mutate`` method)
          4. check if the new population's fittest is fitter than the overall fittest
          4a.  if not and the ``elitist`` option is active, replace the weakest solution
               with the overall fittest
        
        generations:  how many generations to run
        p_mutate:  probability of mutation in [0, 1]
        p_crossover:  probability in [0, 1] that a crossover event will occur for each offspring
        elitist (default=True):  option to replace the weakest solution with the 
                                 strongest if a new one is not found each generation
                                 
        return:  the overall fittest solution (chromosome)
        """
        start_time = time.time()
        
        assert 0 <= p_mutate <= 1
        assert 0 <= p_crossover <= 1
        
        # these values guaranteed to be replaced in first generation
        self.min_fit_ever =  1e999999999
        self.max_fit_ever = -1e999999999
        
        self.generation_fittest.clear()
        self.generation_fittest_fit.clear()
        self.overall_fittest_fit.clear()
        self.new_fittest_generations.clear()
        
        overall_fittest = self.get_fittest()
        overall_fittest_fit = self.eval_fitness(overall_fittest)
        gens_since_upset = 0
    
        for gen in range(1, generations + 1):
            survivors = self.compete(self.chromosomes)
            self.chromosomes = self.reproduce(survivors, p_crossover)
            self.mutate(self.chromosomes, p_mutate)
                    
            # check for new fittest
            gen_fittest = self.get_fittest().copy()
            gen_fittest_fit = self.eval_fitness(gen_fittest)
            
            if gen_fittest_fit > overall_fittest_fit:
                overall_fittest = gen_fittest
                overall_fittest_fit = gen_fittest_fit
                self.new_fittest_generations.append(gen)
                gens_since_upset = 0
            else:
                gens_since_upset += 1
                
                if elitist:
                    # no new fittest found, replace least fit with overall fittest
                    self.sort(self.chromosomes)
                    self.chromosomes[0].dna = overall_fittest.dna
            
            if quit_after and gens_since_upset >= quit_after:
                print("quitting on generation", gen, "after", quit_after, "generations with no upset")
                break
            
            if refresh_after and gens_since_upset >= refresh_after:
                # been a very long time since a new best solution -- mix things up
                print("refreshing on generation", gen)
                self.mutate(self.chromosomes, 0.5)
                gens_since_upset = 0
                
            self.generation_fittest[gen] = gen_fittest
            self.generation_fittest_fit[gen] = gen_fittest_fit
            self.overall_fittest_fit[gen] = overall_fittest_fit
            
            if self.should_terminate(overall_fittest):
                break
            
        self.run_time_s = time.time() - start_time
        
        return overall_fittest
        
    def should_terminate(self, overall_fittest):
        """ 
        Return whether the current run should terminate, called at the end of each generation.
        
        overall_fittest:  best solution/chromosome encountered thus far in the run
        
        return:  bool
        """
        return False

# ga
A lightweight genetic algorithm library written in pure Python.

## API structure

The API consists of 4 main concepts.

  1. genes (genes.py)
  2. chromosomes (chromosomes.py)
  3. translators (translators.py)
  4. genetic algorithms (algorithms.py)
  
Base classes (some abstract) are provided for each of these concepts in the files listed above.
You could just skip down to the genetic algorithms section, but to get the most out of this package, I recommend the full read.

### Genes (genes.py)

* The "abstract" `BaseGene` class represents a single, uninterrupted sequence of genetic material (DNA).
    * NB: `BaseGene` is not technically abstract (it does not use the `abc` package) but should be treated as such.
* You can create genes with DNA randomly generated from the `GENETIC_MATERIAL_OPTIONS` class variable:

        length = 4
        g = BinaryGene.create_random(length)  # BinaryGene.GENETIC_MATERIAL_OPTIONS = '10'
        print(g.dna)
        > 4 random characters such as '0110'
        
* Genes have a `mutate()` method that randomly mutates DNA elements based on a probability:
        
        g = BinaryGene('1111')
        g.mutate(1.0)  # 100% chance that each element will mutate
        print(g.dna)
        > '0000'

* A few concrete gene classes are provided, with the following `GENETIC_MATERIAL_OPTIONS`:
    * `BinaryGene` - 01
    * `Base10Gene` - 0123456789
    * `AlphabetGene` - ABCDEFGHIJKLMNOPQRSTUVWXYZ
    * `DNAGene` - ATCG
    
* To subclass `BaseGene`, override the `GENETIC_MATERIAL_OPTIONS` class variable with a string of supported characters:
  
        class VowelGene(BaseGene):
            GENETIC_MATERIAL_OPTIONS = 'AEIOU'
            
### Chromosomes (chromosomes.py)

* The concrete `Chromosome` class represents an ordered sequence of genes, and facilitates evolutionary mechanisms such as mutation and crossover.
* To create a chromosome with randomly generated genes (`BinaryGene` by default):

        # 1 chromosome with 1 gene
        gene_length = 4
        c = Chromosome.create_random(gene_length)
        print(c.dna)
        > something like Chromosome<1100>
        
* To create a list of randomly generated chromosomes, with multiple genes of different lengths:

        # 30 chromosomes with 5 genes each
        gene_length = (2, 4, 6, 8, 10)
        chromosomes = Chromosome.create_random(gene_length, n=30, gene_class=Base10Gene)
        print(chromosomes[0])
        > something like Chromosome<01,9221,429183,23491832,1937261837>
        
* Chromosomes have a `crossover()` method to exchange DNA with each other:

        c1 = Chromosome([BinaryGene('11')])
        c2 = Chromosome([BinaryGene('00')])
        c1.crossover(c2, 1)  # 1=crossover index
        print(c1, c2)
        > Chromosome<10> Chromosome<01>
        
* They also have a `mutate()` method which calls each gene's `mutate()` method with the same mutation probability.
* NB: if you want to save a snapshot of a chromosome that may change, be sure to call its `copy()` method:

        c1 = Chromosome([BinaryGene('1111')])
        c1_copy = c1.copy()
        c1.mutate(1.0)
        print(c1_copy, c1)  # the original has changed
        > Chromosome<1111> Chromosome<0000>
        
### Translators (translators.py)

* A "translator" is responsible for translating DNA into an object (such as a number or custom class instance) that a genetic algorithm can use to measure the fitness of a solution/chromosome.
* The abstract `BaseTranslator` class provides 2 methods:
    * `translate_gene` (abstract) - translates a single gene into a useful value
    * `translate_chromosome` (concrete) - returns a list of `translate_gene` results for a chromosome

* A few concrete translators are provided:
    * `BinaryIntTranslator` - translates binary DNA into base-10 integers
    * `BinaryFloatTranslator` - translates binary DNA into base-10 floating point real numbers (see its docstrings)
    * `Base10IntTranslator` - translates base-10 DNA into positive base-10 integers
    
* Custom genetic algorithms may require new ways of translating DNA into values/objects. For instance, you might need
to create a new subclass that runs an external program and retrieves results:

        class ModelTranslator(BaseTranslator):
            def translate_gene(self, gene):
                <run a subprocess>
                return <some results>
                
### Genetic algorithms (algorithms.py)

* The abstract `BaseGeneticAlgorithm` class facilitates the execution and evolution of chromosomes toward an optimal state.
* Genetic algorithm instances are given an initial population of chromosomes when created--this population is evolved over successive iterations (aka "generations"):

        generations = 100
        chromosomes = Chromosome.create_random(gene_length=20, n=10)
        ga = MyGA(chromosomes)
        ga.run(generations, ...)

* To create a new GA subclass, one must at least override the `eval_fitness()` method:

        class MostOnesGA(BaseGeneticAlgorithm):
            """ A simple GA to maximize the number of 1's in a solution. """
            
            def eval_fitness(self, chromosome):
                return chromosome.dna.count('1')
                
* To run this GA for 100 generations and retrieve the best solution it encountered:

        p_mutate = 0.15
        p_crossover = 0.60
        chromosomes = Chromosome.create_random(gene_length=20, n=10) 
        ga = MostOnesGA(chromosomes)
        best = ga.run(100, p_mutate, p_crossover)
        print(ga.eval_fitness(best))
        > 20
        
* Most GAs will need a translator to help the `eval_fitness()` method:

        class BiggestIntGA(BaseGeneticAlgorithm):
            """ A simple GA to find the largest encodable integer. """
            
            def eval_fitness(self, chromosome):
                ints = self.translator.translate_chromosome(chromosome)
                return sum(ints)
                
        translator = BinaryIntTranslator()
        ga = BiggestNumberGA(chromosomes, translator=translator)
        
* The `run()` method accepts several arguments that can help improve the success of your search:
    * `p_mutate` - the probability given to chromosome and gene `mutate()` methods
    * `p_crossover` - the probability given to the GA's `reproduce()` method
    * `elitist` (default=True) - at the end of a generation that did not find a new best chromosome, replace the current weakest with the overall strongest from the current run
    * `refresh_after` (default=None) - after N generations of not finding a new best chromosome, call the GA's `refresh()` method to force exploring the solution space
    * `quit_after` (default=None) - after N generations of not finding a new best chromosome, stop the search
* `BaseGeneticAlgorithm` has several methods you could override or otherwise play with as needed. 
* General default behavior for these methods (in the order they are called each generation):
    * `compete()` - simulate survival of the fittest (relative and absolute)
    * `reproduce()` - survivors of competition reproduce (with the potential for genetic cross-over) to restore the population to its original size
    * `mutate()` - check all chromosomes (and all genes) for mutation events
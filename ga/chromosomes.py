import random

from .genes import BaseGene, BinaryGene


class Chromosome:
    """
    Represents a chromosome, a single strand of DNA with
    at least 1 gene. Genes are ordered along the chromosome.
    """
    @classmethod
    def create_random(cls, gene_length, n=1, gene_class=BinaryGene):
        """
        Create 1 or more chromosomes with randomly generated DNA.

        gene_length:  int (or sequence of ints) describing gene DNA length
        n:  number of chromosomes to create (default=1); returns a list if n>1, else a single chromosome
        gene_class:  subclass of ``ga.chromosomes.BaseGene`` to use for genes

        return:  new chromosome
        """
        assert issubclass(gene_class, BaseGene)
        chromosomes = []

        # when gene_length is scalar, convert to a list to keep subsequent code simple
        if not hasattr(gene_length, '__iter__'):
                gene_length = [gene_length]

        for _ in range(n):
            genes = [gene_class.create_random(length) for length in gene_length]
            chromosomes.append(cls(genes))

        if n == 1:
            return chromosomes[0]
        else:
            return chromosomes

    def __init__(self, genes):
        """
        Construct a new ``Chromosome`` instance.
        
        genes:  a collection of ``genes.BaseGene`` instances
        """
        assert all(isinstance(g, BaseGene) for g in genes)
        self.genes = genes
        
    @property
    def dna(self):
        """ Return the full DNA string for all genes in this chromosome. """
        return ''.join(g.dna for g in self.genes)
        
    @dna.setter
    def dna(self, dna):
        """ 
        Replace this chromosome's DNA with new DNA of equal length,
        assigning the new DNA to the chromosome's genes sequentially.
        
        For example, if a chromosome contains these genes...
          1. 100100
          2. 011011
          
          ...and the new DNA is 111111000000, the genes become:
          1. 111111
          2. 000000
        """
        assert self.length == len(dna)
        i = 0
        
        for gene in self.genes:
            gene.dna = dna[i:i + gene.length]
            
            i += gene.length
        
    @property
    def length(self):
        """ Return the length of this chromosome's full DNA string. """
        return len(self.dna)
        
    def crossover(self, chromosome, point):
        """
        Exchange DNA with another chromosome of equal length at a common point.
        
        For example, consider chromosomes:
          1. 11110000
          2. 00001111
          
        If the crossover point is 4, the exchange results in a new DNA arrangement:
          1. 11111111
          2. 00000000
          
        chromosome:  other ``Chromosome`` to exchange DNA with
        point:  zero-based index used for the crossover point
        """
        assert self.length == chromosome.length
        new_dna = self.dna[:point] + chromosome.dna[point:]
        other_new_dna = chromosome.dna[:point] + self.dna[point:]
        
        self.dna = new_dna
        chromosome.dna = other_new_dna
        
    def mutate(self, p_mutate):
        """ 
        Check all genes in this chromosome for mutation. 
        
        p_mutate:  probability for mutation to occur
        """
        assert 0 <= p_mutate <= 1
        
        for gene in self.genes:
            gene.mutate(p_mutate)
            
    def copy(self):
        """ Return a new instance of this chromosome by copying its genes. """
        genes = [g.copy() for g in self.genes]
        return type(self)(genes)
        
    def __iter__(self):
        for g in self.genes:
            yield g
            
    def __str__(self):
        return 'Chromosome<{}>'.format(','.join(g.dna for g in self.genes))
        

class ReorderingSetChromosome(Chromosome):
    """
    A chromosome for completely representing a set of values.
    Each gene represents a distinct value in the set.
    The chromosome must represent each value in the set of choices exactly once.
    Evolutionary mechanisms are limited to reordering the original genes.
    """
    def __init__(self, genes, dna_choices):
        """
        Construct a new ``ReorderingSetChromosome`` instance.

        genes:  genes in this chromosome, which must pass the ``self.check_genes()`` test.
        dna_choices:  collection of unique DNA choices that genes in this chromosome must represent
        """
        super().__init__(genes)
        self.dna_choices = tuple(dna_choices)
        self.dna_choices_set = set(dna_choices)
        assert len(self.dna_choices) == len(self.dna_choices_set)
        self.check_genes()
        
    def check_genes(self):
        """ Assert that every DNA choice is represented by exactly one gene. """
        gene_dna_set = set([g.dna for g in self.genes])
        assert gene_dna_set == self.dna_choices_set
        
    def mutate(self, p_mutate):
        # gene-swapping mutation
        for g1_idx, gene in enumerate(self.genes):
            if random.random() < p_mutate:
                num_genes = len(self.genes)
                g2_idx = g1_idx

                while g1_idx == g2_idx:
                    g2_idx = random.randrange(num_genes)

                tmp = self.genes[g2_idx].copy()
                self.genes[g2_idx] = gene.copy()
                self.genes[g1_idx] = tmp
            
        self.check_genes()
            
    def crossover(self, chromosome, point):
        # find gene on other chromosome at point
        i = 0
        other_gene_idx = 0
        for g in chromosome.genes:
            i += g.length
            if point < i:
                break
                
            other_gene_idx += 1
        other_gene = chromosome.genes[other_gene_idx]
            
        # find idx of gene on this chromosome
        for i, g in enumerate(self.genes):
            if g.dna == other_gene.dna:
                # perform swap
                tmp = self.genes[other_gene_idx].copy()
                self.genes[other_gene_idx] = g.copy()
                self.genes[i] = tmp
                break
                
        self.check_genes()
        
    def copy(self):
        genes = [g.copy() for g in self.genes]
        return ReorderingSetChromosome(genes, self.dna_choices)

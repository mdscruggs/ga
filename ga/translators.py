import abc

from .chromosomes import Chromosome


class BaseTranslator(abc.ABC):
    """
    A "translator" is responsible for translating DNA into an object 
    (such as a number or custom class instance) that a genetic algorithm 
    can use to measure the fitness of a solution/chromosome.
    
    Consider how in real cells, DNA serves as the "raw data" that gets
    transcribed and translated by cellular mechanisms (RNA, ribosomes, etc.) 
    into proteins, which have their own properties, capabilities, and functions.
    
    Translators perform these cellular mechanisms in this library, separating 
    raw chromosomes/genes/DNA from the logic needed to express them.
    """
    @abc.abstractmethod
    def translate_gene(self, gene):
        """
        Translate a gene into the object it represents.
        
        gene:  a ``genes.BaseGene`` instance to translate
        
        return:  the object that the gene's DNA represents
        """
        raise NotImplementedError
        
    def translate_chromosome(self, chromosome):
        """
        Translate all the genes in a chromosome.
        
        chromosome:  a ``chromosomes.Chromosome`` instance to translate
        
        return:  a list of translation products for the genes in the chromosome
        """
        assert isinstance(chromosome, Chromosome)
        return [self.translate_gene(g) for g in chromosome]
        
        
class BinaryIntTranslator(BaseTranslator):
    """
    A translator that translates binary DNA into base-10 integers.
    """
    def translate_gene(self, gene):
        """ Return the base-10 integer represented by a gene's binary DNA. """
        return int(gene.dna, base=2)
        
        
class BinaryFloatTranslator(BaseTranslator):
    """
    A translator that translates binary DNA into base-10 floating point real numbers.
    """
    def __init__(self, significand_length, signed=True):
        """
        Construct a new ``BinaryFloatTranslator``.
        
        significand_length:  the length of the binary DNA that contains an integer.
        signed:  whether genes start with a single sign bit (0=negative, 1=positive)
        """
        self.significand_length = significand_length
        self.signed = signed

    def translate_gene(self, gene):
        """
        Translate a gene with binary DNA into a base-10 floating point real number.
        
        Parses the DNA in this manner:
          1. The first bit determines the sign of the integer portion of the result (0=positive, 1=negative)
          2. The next ``significand_length`` bits of the DNA are converted into a base-10 integer,
             and given a positive/negative sign based on step (1).
          3. The next bit determines the sign of the exponent portion of the result (0=positive, 1=negative)
          4. The remaining bits in the DNA are converted into a base-10 integer, and given a
             positive/negative sign based on step (3).
          5. The result of step (2) is multiplied by 10 raised to the power of the result of step (4).
          
        Example:  let DNA="001111", significand_length=3
          1. "0" indicates a positive sign for the integer portion
          2. "011" is converted into the base-10 integer 3, its sign stays positive due to step (1)
          3. "1" indicates a negative sign for the exponent portion
          4. The remaining "1" bit is converted into the base-10 integer 1 and becomes -1 due to step (3)
          5. The final result becomes:  3 * 10^-1 = 0.3
        """
        if self.signed:
            sign = 1 if gene.dna[0] == '0' else -1
            base_start_idx = 1
        else:
            sign = 1
            base_start_idx = 0

        base = sign * int(gene.dna[base_start_idx:base_start_idx + self.significand_length], base=2)
        
        exponent_sign = 1 if gene.dna[1 + self.significand_length] == '0' else -1
        exponent = exponent_sign * int(gene.dna[self.significand_length + 2:], base=2)
        
        return float(base * 10 ** exponent)
        
        
class Base10IntTranslator(BaseTranslator):
    """
    A translator that translates base-10 DNA into a positive base-10 integer.
    """
    def translate_gene(self, gene):
        return int(gene.dna)
import random
import string


class BaseGene:
    """
    A gene which DNA that represents a single feature or attribute, such as hair color.
    
    The ``GENETIC_MATERIAL_OPTIONS`` static attribute is a string
    containing the set of characters a single DNA element can have.
    Subclasses must override this with a string of characters, usually unique.
    """
    GENETIC_MATERIAL_OPTIONS = ''
    
    @classmethod
    def create_random(cls, length, **kwargs):
        """
        Return a new instance of this gene class with random DNA,
        with characters chosen from ``GENETIC_MATERIAL_OPTIONS``.
        
        length:  the number of characters in the randomized DNA
        **kwargs:  forwarded to the ``cls`` constructor
        """
        dna = ''.join([random.choice(cls.GENETIC_MATERIAL_OPTIONS) for _ in range(length)])
        return cls(dna, **kwargs)

    def __init__(self, dna, suppressed=False, name=None):
        """
        Construct a new BaseGene.
        
        dna:  DNA string, which must only contain characters in ``GENETIC_MATERIAL_OPTIONS``
        suppressed (default=False):  whether expression of this gene is suppressed somehow
                                     (has no effect in this base class; may be useful in subclasses)
        name (default=None):  optional name for this gene
        """
        self._check_dna(dna)
        
        self._dna = dna
        self.suppressed = suppressed
        self.name = name
        
    @property
    def length(self):
        """ Return the length of this gene's DNA string. """
        return len(self.dna)
        
    @property
    def dna(self):
        """ Return this gene's DNA string. """
        return self._dna
        
    @dna.setter
    def dna(self, dna):
        """ 
        Set this gene's DNA string.
        Checks that it only contains characters in ``GENETIC_MATERIAL_OPTIONS``.
        """
        self._check_dna(dna)
        self._dna = dna
        
    def mutate(self, p_mutate):
        """
        Simulate mutation against a probability.
        
        p_mutate:  probability for mutation to occur
        """
        new_dna = []

        for bit in self.dna:
            if random.random() < p_mutate:
                new_bit = bit

                while new_bit == bit:
                    new_bit = random.choice(self.GENETIC_MATERIAL_OPTIONS)
                bit = new_bit

            new_dna.append(bit)

        self.dna = ''.join(new_dna)
        
    def copy(self):
        """ Return a new instance of this gene with the same DNA. """
        return type(self)(self.dna, suppressed=self.suppressed, name=self.name)
        
    def _check_dna(self, dna):
        """ Check that a DNA string only contains characters in ``GENETIC_MATERIAL_OPTIONS``. """
        valid_chars = set(self.GENETIC_MATERIAL_OPTIONS)
        assert all(char in valid_chars for char in dna)
        
    def __str__(self):
        s = 'Gene'
        if self.name:
            s += '[' + self.name + ']'
        return s + '<' + self.dna + '>'


class BinaryGene(BaseGene):
    """
    The most common gene type for genetic algorithms,
    this gene class uses binary DNA string values ("0" and "1").
    """
    GENETIC_MATERIAL_OPTIONS = '01'
        
    def mutate(self, p_mutate):
        """
        Check each element for mutation, swapping "0" for "1" and vice-versa.
        """
        new_dna = []
        
        for bit in self.dna:
            if random.random() < p_mutate:
                bit = '1' if bit == '0' else '0'
                
            new_dna.append(bit)
            
        self.dna = ''.join(new_dna)
        
        
class Base10Gene(BaseGene):
    """
    A gene that uses base-10 numbers for DNA.
    Only supports positive valued integers.
    """
    GENETIC_MATERIAL_OPTIONS = "0123456789"
        
       
class AlphabetGene(BaseGene):
    """
    A gene that uses uppercase alphabet characters for DNA.
    """
    GENETIC_MATERIAL_OPTIONS = string.ascii_uppercase


class DNAGene(BaseGene):
    GENETIC_MATERIAL_OPTIONS = 'ATCG'
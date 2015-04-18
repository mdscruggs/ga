"""
GA solution to reddit daily programmer from 3/18/15:
http://www.reddit.com/r/dailyprogrammer/comments/2zezvf/20150318_challenge_206_intermediate_maximizing/
"""
from io import StringIO
import math
import time

from ..algorithms import BaseGeneticAlgorithm
from ..chromosomes import Chromosome
from ..translators import BinaryIntTranslator

MAP = (  # "x" denotes a crop cell
"""......x...x....x............x............x.................................x...............
.........x...........x...................x.....x...........xx.............x................
...........x.................x.x............x..........................x................x..
......x...x.....................x.....x....x.........x......x.......x...x..................
.x...x.....x................xx...........................x.....xx.....x............x.......
.....xx.......x..x........x.............xx........x.......x.....................x.......x..
...x..x.x..x......x..............................................................x...x.....
........x..x......x......x...x.x....x.......x........x..x...........x.x...x..........xx....
...............x.x....x...........x......x.............x..........................x........
...................x........x..............................................................
..x.x.....................................x..x.x......x......x.............................
......x.............................................................................x..x...
......x....x...............x...............................................................
............x.............x.............................x...............x................x.
..xx........xx............x...x......................x.....................................
........x........xx..............x.....................x.x.......x........................x
.......x....................xx.............................................................
............x...x.........x...xx...............x...........................................
.............................x...............xx..x...........x....x........x...x.......x.x.
..........x.......................x.....................................x..................
...xx..x.x..................x........................x.....................x..x.......x....
.............xx..........x...............x......................x.........x.........x....x.
...............................x.....................x.x...................................
...................x....x............................x...x.......x.............x....x.x....
.x.xx........................x...................................x.....x.......xx..........
.......x...................................................................................
.........x.....x.................x.................x...x.......x..................x........
.......x................x.x...................................x...xx....x.....x...x........
..............................................x..................x.........................
............................x........x.......x............................................x
..x.............x.....x...............x............x...x....x...x..........................
.......................xx.................x...................x...................x.......x
.x.x.............x....x.................................x...........x..x..........x.....x..
...x..x..x......................x...........x..........x.............xxx....x..........x...
...........................................................x...............................
x......x.....x................x...............x....................................x.......
..x...........................x............x..........x....x..............................x
.......................x.......xx...............x...x.x.................x..x............x..
x................x.......x........x.............................x.x.x...................x.x
.......................x...x.......................................................x.......
.x..................x.....x..........................................x...........x.........
.x...................x........x.................x..........xx..................x..x........
.x..........x...x...........................x.x....................x..x.......x............
.............x...x..................x................x..x.x.....xxx..x...xx..x.............
.x...................x.x....x...x.................x.............................x.....x....
......................x.x........x...........x...................................x......x..
................x....................................x....x....x......x..............x..x..
......x.........................................x..x......x.x.......x......................
.x..............................x..........x.x....x.................x......................
x..x...........x..x.x...x..........................................x..............xx.......
..xx......x.......x.x.................x......................................x.............""")


def mapstr_to_list(mapstr):
    """ Convert an ASCII map string with rows to a list of strings, 1 string per row. """
    maplist = []

    with StringIO(mapstr) as infile:
        for row in infile:
            maplist.append(row.strip())

    return maplist


def sprinkler_reaches_cell(x, y, sx, sy, r):
    """
    Return whether a cell is within the radius of the sprinkler.

    x:  column index of cell
    y:  row index of cell
    sx:  column index of sprinkler
    sy:  row index of sprinkler
    r:  sprinkler radius
    """
    dx = sx - x
    dy = sy - y
    return math.sqrt(dx ** 2 + dy ** 2) <= r


class IrrigationGA(BaseGeneticAlgorithm):
    def __init__(self, mapstr, h, w, r, *args, **kwargs):
        """
        Construct a new ``IrrigationGA`` instance.

        Assumes chromosomes have 2 genes with binary DNA.

        map:  ASCII map as a single string, 1 line per row
        h:  height of ASCII crop map (number of rows)
        w:  width of ASCII crop map (number of columns)
        r:  how many cells away the sprinkler can reach
        *args:  forwarded to BaseGeneticAlgorithm constructor
        **kwargs:  forwarded to BaseGeneticAlgorithm constructor
        """
        super().__init__(*args, **kwargs)
        self.mapstr = mapstr
        self.h = h
        self.w = w
        self.r = r

        self.translator = BinaryIntTranslator()

        # parse map string into a list of strings, 1 string per row
        self.maplist = mapstr_to_list(mapstr)

        self.fitness_cache = {}  # DNA -> fitness value

    def eval_fitness(self, chromosome):
        """
        Return the number of plants reached by the sprinkler.

        Returns a large penalty for sprinkler locations outside the map.
        """
        if chromosome.dna in self.fitness_cache:
            return self.fitness_cache[chromosome.dna]

        # convert DNA to represented sprinkler coordinates
        sx, sy = self.translator.translate_chromosome(chromosome)

        # check for invalid points
        penalty = 0
        if sx >= self.w:
            penalty += self.w * self.h
        if sy >= self.h:
            penalty += self.w * self.h

        if penalty > 0:
            self.fitness_cache[chromosome.dna] = -penalty
            return -penalty

        # calculate number of crop cells watered by sprinkler
        crops_watered = 0

        # check bounding box around sprinkler for crop cells
        # circle guaranteed to be within square with side length 2*r around sprinkler
        row_start_idx = max(0, sy - self.r)
        row_end_idx = min(sy + self.r + 1, self.h)
        col_start_idx = max(0, sx - self.r)
        col_end_idx = min(sx + self.r + 1, self.w)

        for y, row in enumerate(self.maplist[row_start_idx:row_end_idx], row_start_idx):
            for x, cell in enumerate(row[col_start_idx:col_end_idx], col_start_idx):
                if cell == 'x' and sprinkler_reaches_cell(x, y, sx, sy, self.r):
                    crops_watered += 1

        # reduce score by 1 if sprinkler placed on a crop cell
        if self.maplist[sy][sx] == 'x':
            crops_watered -= 1

        self.fitness_cache[chromosome.dna] = crops_watered
        return crops_watered

    def map_sprinkler(self, sx, sy, watered_crop='^', watered_field='_', dry_field=' ', dry_crop='x'):
        """
        Return a version of the ASCII map showing reached crop cells.
        """
        # convert strings (rows) to lists of characters for easier map editing
        maplist = [list(s) for s in self.maplist]

        for y, row in enumerate(maplist):
            for x, cell in enumerate(row):
                if sprinkler_reaches_cell(x, y, sx, sy, self.r):
                    if cell == 'x':
                        cell = watered_crop
                    else:
                        cell = watered_field
                else:
                    cell = dry_crop if cell == 'x' else dry_field

                maplist[y][x] = cell

        maplist[sy][sx] = 'O'  # sprinkler

        return '\n'.join([''.join(row) for row in maplist])


def run(generations=500, p_mutate=0.10, p_crossover=0.65):
    # create GA instance
    gene_length = (7, 6)  # 2^6 = 64 > 51; 2^7 = 128 > 91
    chromosomes = Chromosome.create_random(gene_length, n=20)
    irrig_ga = IrrigationGA(MAP, 51, 91, 9, chromosomes)

    # run GA
    t0 = time.time()
    best = irrig_ga.run(generations, p_mutate, p_crossover)
    t1 = time.time()

    # report results
    best_x, best_y = irrig_ga.translator.translate_chromosome(best)
    best_plant_count = irrig_ga.eval_fitness(best)

    print("Best solution is (row, col) ({}, {}) with {} plants reached".format(best_y, best_x, best_plant_count))
    print("took", round(t1-t0, 2), "sec")
    print(irrig_ga.map_sprinkler(best_x, best_y))
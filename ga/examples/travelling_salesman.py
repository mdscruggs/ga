import math
import random

try:
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    plt.style.use('ggplot')
except ImportError:
    plt = None

from ..algorithms import BaseGeneticAlgorithm
from ..chromosomes import ReorderingSetChromosome
from ..genes import BinaryGene
from ..translators import BinaryIntTranslator


class TravellingSalesmanGA(BaseGeneticAlgorithm):
    def __init__(self, city_distances, *args, **kwargs):
        """
        city_distances:  2-deep mapping of city_id -> city_id -> distance
        """
        super().__init__(*args, **kwargs)
        self.city_distances = city_distances

        self.fitness_cache = {}
        self.max_distance = max([max(subdict.values()) for subdict in city_distances.values()])
        self.num_cities = len(self.city_distances)

    def calc_distance(self, chromosome, pow=1):
        # get list of city IDs
        city_ids = self.translator.translate_chromosome(chromosome)

        # compute distance travelled
        tot_dist = 0
        for i, start_city_id in enumerate(city_ids[:-1]):
            end_city_id = city_ids[i + 1]
            tot_dist += self.city_distances[start_city_id][end_city_id] ** pow

        tot_dist += self.city_distances[city_ids[-1]][city_ids[0]] ** pow
        return tot_dist

    def eval_fitness(self, chromosome):
        """
        Calculate the distance travelled by the salesman by converting
        the solution/chromosome into a sequence of visited city IDs.

        Penalty distance is added each time any of these conditions occur:
          1. cities are visited multiple times
          2. not all cities are visited
          3. an invalid city ID is encountered

        return:  fitness value
        """
        if chromosome.dna in self.fitness_cache:
            return self.fitness_cache[chromosome.dna]

        fitness = -1 * self.calc_distance(chromosome, pow=2)
        self.fitness_cache[chromosome.dna] = fitness
        return fitness


def run(num_cities=20, num_chromosomes=20, generations=2500, plot=True):
    # solve a simple travelling salesman problem
    rs = random.randint(1, 1000000)
    random.seed(100)
    
    gene_length = -1
    for i in range(1, num_cities + 1):
        if 2 ** i >= num_cities:
            gene_length = i + 1
            break
    assert gene_length >= 1
    
    city_ids = list(range(1, num_cities + 1))
    city_points = {city_id: (random.random() * 100, random.random() * 100) for city_id in city_ids}
    
    city_distances = {}
    for start_city_id in city_ids:
        city_distances[start_city_id] = {}
        x1, y1 = city_points[start_city_id]
        
        for end_city_id in city_ids:
            x2, y2 = city_points[end_city_id]
            dist = math.sqrt((x2 - x1)**2 + (y2 - y1)**2)
            city_distances[start_city_id][end_city_id] = dist
    
    print("city distances:")
    for start_city_id in city_ids:
        for end_city_id in city_ids:
            print("distance from", start_city_id, "to", end_city_id, "=", city_distances[start_city_id][end_city_id])
    
    random.seed(rs)
    
    chromosomes = []
    for x in range(num_chromosomes):
        genes = []
        
        for city_id in city_ids:
            dna = ('{:0' + str(gene_length) + 'd}').format(int(bin(city_id)[2:]))
            g = BinaryGene(dna, name='city ' + str(x))
            genes.append(g)
            
        choices = [g.dna for g in genes]
        c = ReorderingSetChromosome(genes, choices)
        chromosomes.append(c)
        
    ts_ga = TravellingSalesmanGA(city_distances, chromosomes, 
                                 translator=BinaryIntTranslator(), 
                                 abs_fit_weight=0, rel_fit_weight=1)
    
    p_mutate = 0.10
    p_cross = 0.50
    
    best = ts_ga.run(generations, p_mutate, p_cross, elitist=True, refresh_after=generations/2)
    best_city_ids = ts_ga.translator.translate_chromosome(best)
    best_dist = ts_ga.calc_distance(best)
    
    print("run took", ts_ga.run_time_s, "seconds")
    print("best solution =", best_city_ids)
    print("best distance =", best_dist)

    if plot:
        if plt:
            # plot fitness progression
            plt.plot([v for k, v in sorted(ts_ga.overall_fittest_fit.items())], label='run best')
            plt.plot([v for k, v in sorted(ts_ga.generation_fittest_fit.items())], label='gen best')
            plt.legend(loc='best')
            plt.show()

            fig, ax = plt.subplots()

            def iter_generations():
                for gen in ts_ga.new_fittest_generations:
                    yield gen

            def animate(generation):
                chromosome = ts_ga.generation_fittest[generation]
                ax.clear()

                x, y = [], []
                for city_id, point in city_points.items():
                    x.append(point[0])
                    y.append(point[1])
                ax.plot(x, y, marker='s', linestyle='', label='cities', alpha=0.6)

                # plot optimal route
                chrom_city_ids = ts_ga.translator.translate_chromosome(chromosome)
                dist = round(ts_ga.calc_distance(chromosome), 2)

                ax.set_title("generation " + str(generation) + "\ndistance = " + str(dist))

                for i, start_city_id in enumerate(chrom_city_ids):
                    end_city_idx = i + 1

                    if end_city_idx == num_cities:
                        # distance from last city to first
                        end_city_idx = 0

                    end_city_id = chrom_city_ids[end_city_idx]

                    x1, y1 = city_points[start_city_id]
                    x2, y2 = city_points[end_city_id]
                    mid_x = (x2 - x1) / 2 + x1
                    mid_y = (y2 - y1) / 2 + y1

                    plt.arrow(x1, y1, x2 - x1, y2 - y1, head_width=1.5, fc='k', ec='k', alpha=0.7, linestyle='dotted', length_includes_head=True)
                    plt.text(mid_x, mid_y, str(i + 1))

            ani = animation.FuncAnimation(fig, animate, iter_generations,
                                          repeat=True, interval=1000, repeat_delay=12000)
            plt.show()
        else:
            print("Did not plot example results because matplotlib not installed")

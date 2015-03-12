import math
import random

import pylab as py
py.style.use('ggplot')

from ..genes import BinaryGene
from ..chromosomes import ReorderingSetChromosome
from ..algorithms import TravelingSalesmanGA
from ..translators import BinaryIntTranslator


def run(num_cities=20, generations=2500):
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
    for x in range(20):
        genes = []
        
        for city_id in city_ids:
            dna = ('{:0' + str(gene_length) + 'd}').format(int(bin(city_id)[2:]))
            g = BinaryGene(dna, name='city ' + str(x))
            genes.append(g)
            
        choices = [g.dna for g in genes]
        c = ReorderingSetChromosome(genes, choices)
        chromosomes.append(c)
        
    ts_ga = TravelingSalesmanGA(city_distances, chromosomes, translator=BinaryIntTranslator())
    
    p_mutate = 0.10
    p_cross = 0.50
    
    best = ts_ga.run(generations, p_mutate, p_cross, elitist=True)
    best_city_ids = ts_ga.translator.translate_chromosome(best)
    best_fit = ts_ga.eval_fitness(best)
    
    best_dist = 0
    for i, start_city_id in enumerate(best_city_ids[:-1]):
        end_city_id = best_city_ids[i + 1]
        
        best_dist += city_distances[start_city_id][end_city_id]
    
    print("run took", ts_ga.run_time_s, "seconds")
    print("best solution =", best_city_ids)
    print("best fitness =", best_fit, "dist =", best_dist)
    
    # plot fitness progression
    py.plot(ts_ga.overall_fittest_fit, label='run best')
    py.plot(ts_ga.generation_fittest_fit, label='gen best')
    py.legend(loc='best')
    py.show()
    
    # plot optimal route
    x, y = [], []
    for city_id, point in city_points.items():
        x.append(point[0])
        y.append(point[1])
    py.plot(x, y, marker='s', linestyle='', label='cities', alpha=0.6)
        
    for i, start_city_id in enumerate(best_city_ids):
        end_city_idx = i + 1
        
        if end_city_idx == num_cities:
            # distance from last city to first
            end_city_idx = 0
            
        end_city_id = best_city_ids[end_city_idx]
        
        x1, y1 = city_points[start_city_id]
        x2, y2 = city_points[end_city_id]
        mid_x = (x2 - x1) / 2 + x1
        mid_y = (y2 - y1) / 2 + y1
        
        py.arrow(x1, y1, x2 - x1, y2 - y1, head_width=1.5, fc='k', ec='k', alpha=0.7, linestyle='dotted', length_includes_head=True)
        py.text(mid_x, mid_y, str(i + 1))
        
    py.title("Traveling salesman solution\ndistance = {}".format(round(best_dist, 2)))
    py.xlabel('x ("longitude")')
    py.ylabel('y ("latitude")')
    py.legend(loc='best')
    py.show()

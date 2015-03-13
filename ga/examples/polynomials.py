try:
    import pylab as py
    py.style.use('ggplot')
except ImportError:
    py = None

from ..algorithms import PolyModelGA
from .. import util


def run(coefficients=(0.001, 0.01, 0.1, 1), num_x=10, generations=5000):
    # fit a polynomial equation to expected values
    
    poly_str = ''
    for i, c in enumerate(coefficients[::-1]):
        if i > 1:
            poly_str = '{0:+} x^{1} '.format(c, i) + poly_str
        elif i == 1:
            poly_str = '{0:+} x '.format(c) + poly_str
        else:
            poly_str = '{0:+}'.format(c) + poly_str
            
    poly_str = '$' + poly_str + '$'
    
    # gene lengths are a little complex:
    #   1 bit for significand sign
    #   8 bits for significand body
    #   1 bit for exponent sign
    #   2 bits for exponent body (represent 0, 1, 2)
    significand_length = 8
    gene_length = (1 + significand_length + 1 + 2,) * len(coefficients)  # 1 gene per polynomial coefficient
    chromosomes = util.random_chromosomes(20, gene_length)
    
    poly_ga = PolyModelGA(coefficients, num_x, significand_length, chromosomes, abs_fit_weight=1, rel_fit_weight=0)
    
    p_mutate = 0.15
    p_cross = 0.50
    best = poly_ga.run(generations, p_mutate, p_cross, elitist=True)
    
    best_coeff = poly_ga.translator.translate_chromosome(best)
    best_y = poly_ga.compute_y(best_coeff, num_x)

    print("run took", poly_ga.run_time_s, "seconds")
    print("best solution coefficients =", best_coeff, "error =", poly_ga.compute_err(best_y, best_coeff))

    if py:
        # plot a curve for every solution that caused an upset
        # older solutions have higher transparency
        x = range(1, num_x + 1)

        for new_best_gen in poly_ga.new_fittest_generations:
            gen_fittest = poly_ga.generation_fittest[new_best_gen-1]
            sol_coefficients = poly_ga.translator.translate_chromosome(gen_fittest)
            modeled_y = poly_ga.compute_y(sol_coefficients, num_x)
            alpha = min(0.5, max(0.2, new_best_gen/max(poly_ga.new_fittest_generations)))

            if new_best_gen == poly_ga.new_fittest_generations[-1]:
                py.plot(x, modeled_y, color='green', label='best', marker='x', linewidth=4, markersize=12, alpha=alpha)
            else:
                py.plot(x, modeled_y, color='red', linestyle='--', alpha=alpha)

        py.plot(x, poly_ga.expected_values, color='blue', marker='o', label='expected')
        py.ylim(min(poly_ga.expected_values) - min(poly_ga.expected_values), max(poly_ga.expected_values) * 2)
        py.legend(loc='best')
        py.title('PolyModelGA ({} gen.)\n{}'.format(generations, poly_str))

        py.show()
    else:
        print("Did not plot example results because matplotlib not installed")
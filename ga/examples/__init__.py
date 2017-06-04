__all__ = ['biggest_multiple', 'polynomials', 'travelling_salesman', 'irrigation']

from . import biggest_multiple
from . import polynomials
from . import travelling_salesman
from . import irrigation


def run_all(plot=True, seed=None):
    """ Run all examples. """
    if seed is not None:
        import random
        random.seed(seed)

    print("Running biggest_multiple.py")
    biggest_multiple.run(plot=plot)

    print("Running polynomials.py")
    polynomials.run(plot=plot)

    print("Running travelling_salesman.py")
    travelling_salesman.run(plot=plot)

    print("Running irrigation.py")
    irrigation.run()
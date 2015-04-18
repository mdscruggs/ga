__all__ = ['biggest_multiple', 'polynomials', 'travelling_salesman', 'irrigation']

from . import biggest_multiple
from . import polynomials
from . import travelling_salesman
from . import irrigation


def run_all(plot=True):
    """ Run all examples. """
    print("Running biggest_multiple.py")
    biggest_multiple.run(plot=plot)

    print("Running polynomials.py")
    polynomials.run(plot=plot)

    print("Running travelling_salesman.py")
    travelling_salesman.run(plot=plot)

    print("Running irrigation.py")
    irrigation.run()
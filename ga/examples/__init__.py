__all__ = ['biggest_multiple', 'polynomials', 'travelling_salesman', 'irrigation']

from . import biggest_multiple
from . import polynomials
from . import travelling_salesman
from . import irrigation


def run_all(plot=True):
    """ Run all examples. """
    biggest_multiple.run(plot=plot)
    polynomials.run(plot=plot)
    travelling_salesman.run(plot=plot)
    irrigation.run()
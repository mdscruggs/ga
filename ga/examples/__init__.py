__all__ = ['biggest_multiple', 'polynomials', 'travelling_salesman', 'irrigation']

from . import biggest_multiple
from . import polynomials
from . import travelling_salesman
from . import irrigation


def run_all():
    """ Run all examples listed in __all__. Assumes they all got imported globally. """
    globals_copy = dict(globals())
    
    for modname in __all__:
        print("Running", modname)
        globals_copy[modname].run()
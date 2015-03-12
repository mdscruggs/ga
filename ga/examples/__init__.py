__all__ = ['biggest_multiple', 'polynomials', 'travelling_salesman.py']

from . import biggest_multiple
from . import polynomials
from . import travelling_salesman


def run_all():
    globals_copy = dict(globals())
    
    for modname in __all__:
        print("Running", modname)
        globals_copy[modname].run()
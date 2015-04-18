from argparse import ArgumentParser

from . import run_all


if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-p', '--plot', help='Show matplotlib plots for examples that have them', action='store_true')
    args = parser.parse_args()

    run_all(plot=args.plot)
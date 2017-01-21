from numpy import *

from consts import (HM_FILE, FLY_FILE, WORM_FILE,
                    SPECIES, GRAPH_FILE)

def save_cut_file():
    for fn in [HM_FILE, FLY_FILE, WORM_FILE]:
        arr = loadtxt(fn, delimiter=',', skiprows=1)
        idx = fn.rfind('.')
        save(fn[:idx], arr)


def save_graph_file():
    for sp in SPECIES:
        fn = GRAPH_FILE[sp]
        arr = loadtxt(fn)
        idx = fn.rfind('.')
        save(fn[:idx], arr)


save_graph_file()

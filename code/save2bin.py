from numpy import *

from consts import HM_FILE, FLY_FILE, WORM_FILE

for fn in [HM_FILE, FLY_FILE, WORM_FILE]:
    arr = loadtxt(fn, delimiter=',', skiprows=1)
    idx = fn.rfind('.')
    save(fn[:idx], arr)

# visualize a chosen quantity for each co-latitudes theta
# pass model file name and quantity as parameters

from ester import *
import matplotlib.pyplot as plt
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "model",
    type=str,
    help="the path to the folder containing ESTER model",
)
parser.add_argument(
    "--plot",
    default="rho",
    type=str,
    help=f"which parameter to plot"
)
ARGS = parser.parse_args()

a=star2d(ARGS.model)
plt.plot(a.r[:],np.log10(getattr(a, ARGS.plot)[:]))
plt.xlabel("r")
plt.ylabel(ARGS.plot)
plt.show()

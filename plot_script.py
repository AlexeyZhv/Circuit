import numpy as np
import scipy.optimize as opt
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
from colorsys import hsv_to_rgb as hsv
import pathlib

def line(x, a, b):
    return a + b * x

font = {'weight' : 'normal',
        'size'   : 14}

matplotlib.rc('font', **font)

amp_path = pathlib.Path("./in_out/output/ampermeters")
vol_path = pathlib.Path("./in_out/output/voltmeters")

for file in amp_path.iterdir():
    plt.figure(figsize=[10, 7], dpi=300)
    df = pd.read_csv(str(file))
    plt.plot(df["t"][2:-1], df["v"][2:-1], color="red", lw=1)
    plt.xlabel("$Time$, sec")
    plt.ylabel("$Current$, A")
    plt.title("Ampermeter " + str(file)[-5])
    plt.savefig(str(file)[:-26] + "plots\\" + str(file)[-26:-3] + "png")

for file in vol_path.iterdir():
    plt.figure(figsize=[10, 7], dpi=300)
    df = pd.read_csv(str(file))
    plt.plot(df["t"][2:-1], df["v"][2:-1], color="blue", lw=1)
    plt.xlabel("$Time$, sec")
    plt.ylabel("$Current$, V")
    plt.title("Voltmeter " + str(file)[-5])
    plt.savefig(str(file)[:-25] + "plots\\" + str(file)[-25:-3] + "png")

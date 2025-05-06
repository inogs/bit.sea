import argparse

import matplotlib

matplotlib.use("Agg")
import numpy as np
import matplotlib.pyplot as plt
from bitsea.commons.utils import addsep


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Generates plots comparison of 4PFTs with HPLC dataset
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument("--inputdir", "-i", type=str, required=True, help="")
    parser.add_argument("--outdir", "-o", type=str, required=True, help="")

    parser.add_argument(
        "--ptype",
        "-p",
        type=str,
        required=True,
        help="Phytoplankton Functional Type",
    )
    return parser.parse_args()


args = argument()

INDIR = addsep(args.inputdir)
OUTDIR = addsep(args.outdir)
P_type = args.ptype


if P_type == "P1l":
    P_name = "DIATO"
if P_type == "P2l":
    P_name = "NANO"
if P_type == "P3l":
    P_name = "PICO"
if P_type == "P4l":
    P_name = "DINO"

infile_MEAN = INDIR + P_type + "_OPTION3_newSeasons3_means.csv"
infile_STD = INDIR + P_type + "_OPTION3_newSeasons3_std.csv"


headers = ["#", "meanWref", "meanWmod", "meanSref", "meanSmod"]
headers_STD = ["#", "meanWref", "meanWmod", "meanSref", "meanSmod"]
layers = ["0-10", "10-30", "30-60", "60-100 ", "100-150", "150-300", "0-300"]
layer_mean_depth = ["5", "20", "50", "80", "125", "225", "150"]
depths = [5, 20, 50, 80, 125, 225, 150]


Pl = np.loadtxt(infile_MEAN, skiprows=1, usecols=range(1, 5))
Pl_std = np.loadtxt(infile_STD, skiprows=1, usecols=range(1, 5))


fig, (ax1, ax2) = plt.subplots(1, 2)

ax1.plot(Pl[:-1, 0], layers[:-1], ":", color="blueviolet")
ax1.plot(Pl[:-1, 1], layers[:-1], "-", color="cyan")
ax2.plot(Pl[:-1, 2], layers[:-1], ":", color="coral")
ax2.plot(Pl[:-1, 3], layers[:-1], "-", color="orange")


ax1.fill_betweenx(
    layers[:-1],
    Pl[:-1, 0] - Pl_std[:-1, 0],
    Pl[:-1, 0] + Pl_std[:-1, 0],
    color="lavender",
)
ax1.fill_betweenx(
    layers[:-1],
    Pl[:-1, 1] - Pl_std[:-1, 1],
    Pl[:-1, 1] + Pl_std[:-1, 1],
    color="lightcyan",
)
ax2.fill_betweenx(
    layers[:-1],
    Pl[:-1, 2] - Pl_std[:-1, 2],
    Pl[:-1, 2] + Pl_std[:-1, 2],
    color="mistyrose",
)
ax2.fill_betweenx(
    layers[:-1],
    Pl[:-1, 3] - Pl_std[:-1, 3],
    Pl[:-1, 3] + Pl_std[:-1, 3],
    color="moccasin",
)

ax1.set_xlim([0, 0.25])
ax2.set_xlim([0, 0.25])

# Set the tick positions
ax1.set_yticks(layers[:-1])
# Set the tick labels
ax1.set_yticklabels(layer_mean_depth[:-1])

ax1.tick_params(axis="both", which="major", labelsize=10)
ax2.tick_params(axis="both", which="major", labelsize=10)

ax2.set_yticklabels([])
ax1.set_xlabel("$[ mg {\,} m^{-3} ]$").set_fontsize(12)
ax2.set_xlabel("$[ mg {\,} m^{-3} ]$").set_fontsize(12)

ax1.set_ylabel("depth [m]").set_fontsize(12)


ax1.invert_yaxis()
ax2.invert_yaxis()
ax1.legend(["HPLC", "mod"])
ax2.legend(["HPLC", "mod"])
ax1.set_title("WINTER")
ax2.set_title("SUMMER")
fig.suptitle(P_name)
fig.savefig(OUTDIR + "HPLC_profiles_comparison_STD_" + P_type + "_seasons.png")


plt.close(fig)

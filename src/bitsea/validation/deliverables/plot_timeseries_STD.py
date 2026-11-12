import argparse

import matplotlib

from bitsea.utilities.argparse_types import date_from_str
from bitsea.utilities.argparse_types import existing_dir_path

matplotlib.use("Agg")
import matplotlib.pyplot as pl
import matplotlib.dates as mdates
from bitsea.validation.deliverables import netcdf_validation_file
from matplotlib.ticker import FormatStrFormatter
from bitsea.instruments.var_conversions import SAT_VARS
from bitsea.commons.Timelist import TimeInterval
from datetime import timedelta


def argument():
    parser = argparse.ArgumentParser(
        description="""
    Plot timeseries for QUID.
    Reference doc: CMEMS-Med-QUID-006-008-V2-V1.0.docx
    Usually Figure IV.2
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--outdir",
        "-o",
        type=existing_dir_path,
        required=True,
        help=""" Output image dir""",
    )
    parser.add_argument(
        "--inputdir", "-i", type=existing_dir_path, required=True
    )
    parser.add_argument(
        "--var",
        "-v",
        type=str,
        required=True,
        choices=[
            "P_l",
            "kd490",
            "P1l",
            "P2l",
            "P3l",
            "P4l",
            "RRS412",
            "RRS443",
            "RRS490",
            "RRS510",
            "RRS555",
            "RRS670",
        ],
        help=""" model var name""",
    )
    parser.add_argument(
        "--datestart",
        "-s",
        type=date_from_str,
        required=True,
        help=""" Date start for time interval to consider for validation, format yyyymmdd""",
    )
    parser.add_argument(
        "--dateend",
        "-e",
        type=date_from_str,
        required=True,
        help=""" Date end for time interval to consider for validation,format yyyymmdd""",
    )
    parser.add_argument(
        "--coastness",
        "-c",
        type=str,
        required=True,
        choices=["coast", "open_sea", "everywhere"],
    )
    parser.add_argument(
        "--zone",
        "-z",
        type=str,
        required=False,
        default="Med",
        help="Areas to generate the STATISTICS mean. std, bias and RMSD with respect satellite",
        choices=["Med", "rivers", "coastal","mer"],
    )

    return parser.parse_args()


args = argument()


TI = TimeInterval.fromdatetimes(args.datestart, args.dateend)
dr = netcdf_validation_file.dir_reader(
    TI, args.inputdir, args.var, args.coastness
)

if args.zone == "Med":
    from bitsea.basins import V2 as OGS
if args.zone == "rivers":
    from bitsea.basins import RiverBoxes as OGS
if args.zone == "coastal":
    from bitsea.basins import COASTAL12nm as OGS
if args.zone == "mer":
    domain = "ION"
    from bitsea.basins.basin import ComposedBasin, SimplePolygonalBasin
    from bitsea.basins.region import Rectangle
    mer = Rectangle(lonmin=4.0, lonmax=20, latmin=32.0, latmax=48.0)
    mer = SimplePolygonalBasin(domain, mer, f"{domain} Sea")
    OGS = ComposedBasin("P", [mer], f"{domain} Sea")
    OGS.P = ComposedBasin("P", [mer], f"{domain} Sea")


if args.var == "kd490":
    units = "[m$^{-1}$]"

elif args.var.startswith("RRS"):
    units = "[st$^{-1}$]"

else:
    units = "[mg/m$^3$]"

var_label = SAT_VARS[args.var] + " " + units


def get_vmin(var: str) -> float:
    if var == "kd490":
        vmin = 0.02
    else:
        vmin = 0.0
    return vmin


def get_vmax(var: str, sub: str, zone: str) -> float:
    if var.startswith("RRS"):
        vmax = 0.0025
        if var.startswith(("RRS5", "RRS6")):
            vmax = 0.01
    if var == "kd490":
        vmax = 0.12
    if var in ["P1l", "P2l", "P3l", "P4l", "P_l"]:
        if zone == "rivers":
            vmax = 1.2
            if sub.startswith("Po"):
                vmax = 4.0
            if sub.startswith("Piave"):
                vmax = 2.0
        else:
            if sub in ["alb", "swm1", "swm2", "nwm", "tyr1", "tyr2"]:
                vmax = 0.9
            else:
                vmax = 0.3
    return vmax


vmin = get_vmin(args.var)

color = "tab:green"
lightcolor = "palegreen"


if args.var == "kd490":
    color = "tab:blue"
    lightcolor = "lightsteelblue"


for isub, sub in enumerate(OGS.P):
    if sub.name == "atl":
        continue
    outfilename = args.outdir / f"{args.var}_{sub.name}_STD.png"
    print(outfilename)
    fig, ax = pl.subplots()
    fig.set_size_inches(12, 4)
    fig.set_dpi(220)
    ax.plot(dr.TIMES, dr.SAT___MEAN[:, isub], "-", label=" Sat", color=color)
    ax.plot(dr.TIMES, dr.SAT___MEAN[:, isub], "o", color=color)
    yfillbottom = dr.SAT___MEAN[:, isub] - dr.SAT____STD[:, isub]
    yfilltop = dr.SAT___MEAN[:, isub] + dr.SAT____STD[:, isub]
    ax.fill_between(dr.TIMES, yfillbottom, yfilltop, color=lightcolor, alpha=0.7)
    ax.plot(dr.TIMES, yfillbottom, "-", color=color, linewidth=0.6)
    ax.plot(dr.TIMES, yfilltop, "-", color=color, linewidth=0.6)
    ax.plot(dr.TIMES, dr.MODEL_MEAN[:, isub], "-k", label=" Model")
    ax.plot(
        dr.TIMES,
        dr.MODEL_MEAN[:, isub] - dr.MODEL__STD[:, isub],
        "-k",
        linewidth=0.6,
    )
    ax.plot(
        dr.TIMES,
        dr.MODEL_MEAN[:, isub] + dr.MODEL__STD[:, isub],
        "-k",
        linewidth=0.6,
    )

    ax.set_ylabel("%s - %s" % (sub.name.upper(), var_label)).set_fontsize(11)
    ax.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))
    ax.legend(
        loc="best",
        labelspacing=0,
        handletextpad=0,
        borderpad=0.1,
        fontsize=10,
        frameon=True,
    )
    leg = pl.gca().get_legend()
    ltext = leg.get_texts()
    pl.setp(ltext, fontsize=10)
    pl.rc("xtick", labelsize=10)
    pl.rc("ytick", labelsize=10)
    ax.tick_params(axis="both", labelsize=10)

    vmax = get_vmax(args.var, sub.name, args.zone)
    pl.ylim(vmin, vmax)
    if dr.TIMES[-1] - dr.TIMES[0] > timedelta(days=90):
        ax.xaxis.set_major_locator(mdates.MonthLocator())
        ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%Y"))
    ax.grid(True)
    xlabels = ax.get_xticklabels()
    pl.setp(xlabels, rotation=30, fontsize=8)
    ylabels = ax.get_yticklabels()
    pl.setp(ylabels, fontsize=8)
    pl.tight_layout()
    pl.savefig(outfilename, dpi=220)
    pl.close(fig)

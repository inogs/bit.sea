import argparse

import numpy as np

from bitsea.basins import V2
from bitsea.commons.dataextractor import DataExtractor
from bitsea.commons.mask import Mask
from bitsea.commons.submask import SubMask
from bitsea.commons.Timelist import TimeList
from bitsea.commons.utils import addsep


def argument():
    parser = argparse.ArgumentParser(
        description="""
    post of check and mis for BGC-Argo DA
    """,
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "--indadir",
        "-d",
        type=str,
        default=None,
        required=True,
        help="DA dir DA__FREQ_1",
    )

    parser.add_argument(
        "--outdir", "-o", type=str, default=None, required=True, help=""
    )

    parser.add_argument(
        "--maskfile", "-m", type=str, default=None, required=True, help=""
    )

    return parser.parse_args()


args = argument()


INDADIR = addsep(args.indadir)
OUTDIR = addsep(args.outdir)
TheMask = Mask.from_file(args.maskfile)

_, jpj, jpi = TheMask.shape
mask200 = TheMask.mask_at_level(200)


dtype = [(sub.name, bool) for sub in V2.P]
npSub = {}
SUB = np.zeros((jpj, jpi), dtype=dtype)
Nsub = 1
for sub in V2.Pred:
    Nsub += 1
    sbmask = SubMask(sub, TheMask)
    SUB[sub.name] = sbmask[0, :, :]
    SUB["med"] = SUB["med"] | SUB[sub.name]
    npSub[sub.name] = np.sum(SUB[sub.name])

npSub["med"] = np.sum(SUB["med"])


TLmis = TimeList.fromfilenames(
    None, INDADIR, "????????.chl_mis.nc", prefix="", dateformat="%Y%m%d"
)


for misfile, datemis in zip(TLmis.filelist, TLmis.Timelist):
    date8 = datemis.strftime("%Y%m%d")
    misfmean = np.zeros((4, Nsub))
    misfmean[:] = np.nan
    De = DataExtractor(TheMask, filename=misfile, varname="misfchl", dimvar=2)
    chlmisf = De.values
    chlmisf[chlmisf > 1.0e19] = np.nan

    for isub, sub in enumerate(V2.P):
        chlsub = chlmisf[SUB[sub.name]]
        if not np.all(np.isnan(chlsub)):
            misfmean[0, isub] = (np.nanmean(chlsub**2)) ** 0.5
        misfmean[3, isub] = float(np.sum(np.isfinite(chlsub))) / npSub[sub.name]
        chlsub = chlmisf[SUB[sub.name] & mask200]
        if not np.all(np.isnan(chlsub)):
            misfmean[1, isub] = (np.nanmean(chlsub**2)) ** 0.5
        chlsub = chlmisf[SUB[sub.name] & ~mask200]
        if not np.all(np.isnan(chlsub)):
            misfmean[2, isub] = (np.nanmean(chlsub**2)) ** 0.5

    fileout = OUTDIR + date8 + "_meanmisf.npy"
    np.save(fileout, misfmean)

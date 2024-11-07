import numpy as np
import netCDF4 as NC
def side_tmask(side, mask):
    if side == "N":
        tmask = mask.mask[:, -1, :]
    if side == "S":
        tmask = mask.mask[:, 0, :]
    if side == "E":
        tmask = mask.mask[:, :, -1]
    if side == "W":
        tmask = mask.mask[:, :, 0]
    return tmask


def zeroPadding(side, mask):
    if side in ["E", "W"]:
        return np.zeros((mask.jpk, mask.jpj), dtype=np.float32)
    if side in ["S", "N"]:
        return np.zeros((mask.jpk, mask.jpi), dtype=np.float32)


def writeSideCheckFile(OUTPUTDIR, M, Mask2, side, t, var, interpdir):
    tmask2 = side_tmask(side, Mask2)
    checkfile = OUTPUTDIR / ("CHECK/OBC_" + side + "." + t + "." + var + ".nc")
    Mcheck = M.copy()
    if interpdir is not None:
        missing_value = 1.0e20
        Mcheck[~tmask2] = missing_value
    else:
        missing_value = 0

    NCout = NC.Dataset(checkfile, "w")
    NCout.createDimension("Lon", Mask2.jpi)
    NCout.createDimension("Lat", Mask2.jpj)
    NCout.createDimension("Depth", Mask2.jpk)
    if side in ["E", "W"]:
        ncvar = NCout.createVariable(var, "f", ("Depth", "Lat"))
    if side in ["N", "S"]:
        ncvar = NCout.createVariable(var, "f", ("Depth", "Lon"))
    setattr(ncvar, "missing_value", missing_value)
    ncvar[:] = Mcheck
    NCout.close()
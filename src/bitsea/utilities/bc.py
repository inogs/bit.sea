from enum import Enum
from typing import Literal
from typing import Union

import netCDF4
import numpy as np
from bitsea.commons.mask import Mask


class Side(Enum):
    """A Side object is typically used to specify which side of the boundary
    conditions is being addressed."""

    NORTH = "N"
    EAST = "E"
    SOUTH = "S"
    WEST = "W"


def side_tmask(
    side: Union[Side, Literal["N", "S", "W", "E"]], mask: Mask
) -> np.ndarray:
    """
    Extracts a boundary slice from the specified side of a mask.

    Args:
        side (`Side`, Literal["N", "S", "W", "E"]): The cardinal direction
          (North, South, West, or East) representing the boundary to retrieve
          from the mask.
        mask (`Mask`): The mask from which to extract the boundary slice.

    Returns:
        `np.ndarray`: A 2D NumPy array containing the values from the specified
        boundary of the mask.
    """
    if isinstance(side, str):
        side = Side(side)

    if side == Side.NORTH:
        return mask.mask[:, -1, :]
    if side == Side.SOUTH:
        return mask.mask[:, 0, :]
    if side == Side.EAST:
        return mask.mask[:, :, -1]
    if side == Side.WEST:
        return mask.mask[:, :, 0]

    raise ValueError(f"Invalid side: {side}")


def zeroPadding(
    side: Union[Side, Literal["N", "S", "W", "E"]], mask: Mask
) -> np.ndarray:
    """
    Creates a zero-filled array matching the shape of a specified boundary of
    the mask.

    Args:
        side (Side, Literal["N", "S", "W", "E"]): The cardinal direction
          representing the boundary of interest.
        mask (Mask): The mask whose boundary shape will be used.

    Returns:
        np.ndarray: A 2D NumPy array filled with zeros, matching the shape of
        the specified boundary of the mask.
    """
    if isinstance(side, str):
        side = Side(side)

    if side in (Side.EAST, Side.WEST):
        return np.zeros((mask.jpk, mask.jpj), dtype=np.float32)
    if side in (Side.SOUTH, Side.NORTH):
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

    NCout = netCDF4.Dataset(checkfile, "w")
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

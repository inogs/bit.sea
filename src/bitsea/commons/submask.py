from numbers import Real
from typing import Optional
from typing import Union

import numpy as np

from bitsea.basins.basin import Basin
from bitsea.basins.region import Rectangle
from bitsea.basins.region import Region
from bitsea.commons.mask import Mask


class SubMask(Mask):
    def __init__(self, basin: Union[Basin, Region], mask: Mask):
        basin_mask = basin.is_inside(lon=mask.xlevels, lat=mask.ylevels)

        super().__init__(
            grid=mask.grid,
            zlevels=mask.zlevels,
            mask_array=np.logical_and(basin_mask, mask.get_sea_cells()),
            e3t=mask.e3t,
        )
        self._original_mask = mask
        self._basin = basin

    @property
    def basin(self):
        return self._basin

    def get_sea_cells(self):
        return self._original_mask.get_sea_cells()

    def get_water_cells(self):
        return self._original_mask.get_water_cells()

    @staticmethod
    def from_square_cutting(
        mask: Mask,
        degrees: Real,
        start_lon: Optional[Real],
        start_lat: Optional[Real],
    ):
        """
        Creates a list of SubMask objects from cutting a Mask object into
        square sections.  The cutting starts from a point and proceeds to cut
        along longitude then along latitude.

        Args:
            - *mask*: a Mask object.
            - *degrees*: number of degrees for a sigle section side.
            - *start_lon*: starting longitude point.
            - *start_lat*: starting latitude point.

        Returns: a list of basin objecs mapping a single section each.
        """
        # Get mask dimensions
        min_lon = mask.xlevels.min()
        max_lon = mask.xlevels.max()
        min_lat = mask.ylevels.min()
        max_lat = mask.ylevels.max()

        # Compute maximum value for degrees
        max_deg = min([(max_lon - min_lon), (max_lat - min_lat)])
        if degrees <= 0:
            raise ValueError("degrees must be greater than 0.")

        if degrees > max_deg:
            raise ValueError(
                f"The degrees value of {degrees} is too big for this mask (maximum: {max_deg})"
            )

        if start_lon is None:
            start_lon = min_lat - min(degrees / 100.0, 1e-8)
        if start_lon is None:
            start_lon = min_lon - min(degrees / 100.0, 1e-8)

        output = list()

        # Bottom Left point
        bl_point = [start_lon, start_lat]
        # Top Right point
        tr_point = [start_lon + degrees, start_lat + degrees]

        # Section indices
        lon_in = 0
        lat_in = 0

        # Latitude cycle
        while tr_point[1] <= max_lat:
            # Longitude cycle
            while tr_point[0] <= max_lon:
                # Create the Rectangle
                rect = Rectangle(
                    bl_point[0], tr_point[0], bl_point[1], tr_point[1]
                )

                # Create the SubMask and append it to output
                output.append(rect)
                # Increment longitude index
                lon_in += 1
                # Increment longitude
                bl_point[0] += degrees
                tr_point[0] += degrees
            # Increment latitude index
            lat_in += 1
            # Increment latitude
            bl_point[1] += degrees
            tr_point[1] += degrees
            # Reset Longitude index
            lon_in = 0
            # Reset Longitude
            bl_point[0] = start_lon
            tr_point[0] = start_lon + degrees
        return output


if __name__ == "__main__":
    # submask.nc generator using no-repetition technique
    # each cell can belong to an only subbasin

    TheMask = Mask(
        "/Users/gbolzon/Documents/workspace/ogstm_boundary_conditions/masks/meshmask_872.nc"
    )
    from bitsea.basins import V2

    already_assigned = np.zeros(TheMask.shape, dtype=bool)

    for sub in V2.Pred.basin_list:
        S = SubMask(sub, maskobject=TheMask)

        #         fig,ax = mapplot({'data':S.mask[0,:,:], 'clim':[0,1]}, mask=TheMask)
        #         ax.set_xlim([-9, 36])
        #         ax.set_ylim([30, 46])
        #         outfile=sub.name + ".png"
        #         fig.savefig(outfile)
        #         continue
        S.mask[already_assigned] = False
        S.save_as_netcdf("submask.nc", maskvarname=sub.name)
        already_assigned[S.mask] = True

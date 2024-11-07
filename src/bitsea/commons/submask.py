from itertools import product as cart_prod
from numbers import Real
from typing import Iterable
from typing import List
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
    ) -> List[Rectangle]:
        """
        Divides a `Mask` object into square sections, returning a list of
        `Rectangle` objects representing each section.
        The division starts from a specified starting point and proceeds along
        the longitude, followed by the latitude.

        Args:
            mask (Mask): The `Mask` object to be divided.
            degrees (Real): The side length of each square section,
              specified in degrees.
            start_lon (Optional[Real]): The longitude coordinate to start the
              division from. If not provided, the minimum longitude of the
              `Mask` object is used.
            start_lat (Optional[Real]): The latitude coordinate to start the
              division from. If not provided, the minimum latitude of the
              `Mask` object is used.

        Returns:
            List[Rectangle]: A list of `Rectangle` objects representing the
            square sections that collectively cover the entire area of the
            original `Mask` object.
        """
        if degrees <= 0:
            raise ValueError("degrees must be greater than 0.")

        # Get mask dimensions
        min_lon = mask.xlevels.min()
        max_lon = mask.xlevels.max()
        min_lat = mask.ylevels.min()
        max_lat = mask.ylevels.max()

        # Compute maximum value for degrees
        max_deg = min([(max_lon - min_lon), (max_lat - min_lat)])

        if degrees > max_deg:
            raise ValueError(
                f"The degrees value of {degrees} is too big for this mask "
                f"(maximum: {max_deg})"
            )

        if start_lon is None:
            start_lon = min_lon - min(degrees / 100.0, 1e-8)
        if start_lat is None:
            start_lat = min_lat - min(degrees / 100.0, 1e-8)

        lon_points = []
        lat_points = []

        current_lon = start_lon
        while current_lon < max_lon:
            lon_points.append(current_lon)
            current_lon += degrees

        current_lat = start_lat
        while current_lat < max_lat:
            lat_points.append(current_lat)
            current_lat += degrees

        output = []
        for square_lon, square_lat in cart_prod(lon_points, lat_points):
            rect = Rectangle(
                square_lon,
                square_lon + degrees,
                square_lat,
                square_lat + degrees,
            )
            output.append(rect)

        return output

    @staticmethod
    def build_partition(
        mask: Mask, regions: Iterable[Union[Basin, Region]]
    ) -> List:
        """
        Creates a list of submasks from a given `Mask` object, with each
        submask corresponding to a specified `Basin` or `Region`.
        The submasks are mutually exclusive, meaning that each point in the
        `Mask` belongs to only one submask. If a point lies
        within multiple basins or regions, it is assigned to the first one in
        the input order.

        Args:
            mask (Mask): The `Mask` object to be subdivided.
            regions (Iterable[Union[Basin, Region]]): An iterable of `Basin`
            or `Region` objects defining the areas for the submasks.

        Returns:
            List[SubMask]: A list of `SubMask` objects, one for each `Basin`
            or `Region` in the input iterable, with no overlapping areas.
        """
        not_assignable = ~mask[:]
        submasks = []

        for region in regions:
            submask = SubMask(region, mask)
            submask[not_assignable] = False
            not_assignable[submask] = True

            submasks.append(submask)

        return submasks

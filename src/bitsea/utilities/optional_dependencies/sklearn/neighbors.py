from logging import getLogger
from warnings import warn

import numpy as np

from bitsea.commons.geodistances import compute_great_circle_distance
from bitsea.utilities.optional_dependencies import OptionalDependencyMissing

LOGGER = getLogger(__name__)


try:
    import sklearn.neighbors

    _SKLEARN_AVAILABLE = True
except ModuleNotFoundError as e:
    LOGGER.debug('Module "sklearn.neighbors" not found.', exc_info=e)
    _SKLEARN_AVAILABLE = False


class BallTreeFallBack:
    WARNING_RAISED = False

    def __init__(self, data, metric="haversine"):
        if not self.WARNING_RAISED:
            warn(
                "sklearn.neighbors.BallTree is not available. "
                "This code will use a fallback implementation that is way "
                "slower. If you want to suppress this warning, install"
                "sklearn",
                OptionalDependencyMissing,
            )
        self.WARNING_RAISED = True

        self._data = np.rad2deg(data)
        self._metric = metric
        if self._metric != "haversine":
            raise ValueError(
                "Only haversine metric is supported by the fallback"
                "implementation of the BallTree. Install sklearn to use a "
                "different kind of metric."
            )

    def query(self, query_data, k: int = 1):
        if k != 1:
            raise ValueError(
                "Only k=1 is supported by the fallback implementation"
            )
        if query_data.ndim != 2:
            raise ValueError(
                f"Only 2D arrays are supported; query_data has shape"
                f"{query_data.shape}"
            )
        n_points = query_data.shape[0]
        output_distance = np.empty((n_points,), dtype=self._data.dtype)
        output_index = np.empty((n_points,), dtype=int)

        query_data_deg = np.rad2deg(query_data)

        for i in range(n_points):
            distances = compute_great_circle_distance(
                lat1=self._data[:, 0],
                lon1=self._data[:, 1],
                lat2=query_data_deg[i, 0],
                lon2=query_data_deg[i, 1],
            )
            output_index[i] = np.argmin(distances)
            output_distance[i] = distances[output_index[i]]

        return output_distance, output_index


if _SKLEARN_AVAILABLE:
    BallTree = sklearn.neighbors.BallTree
else:
    BallTree = BallTreeFallBack

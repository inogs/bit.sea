import inspect
from logging import getLogger
from warnings import warn

from bitsea.utilities.optional_dependencies import OptionalDependencyMissing

LOGGER = getLogger(__name__)

try:
    import geographiclib.geodesic

    GEOGRAPHICLIB_AVAILABLE = True
    _IMPORT_ERROR = None
except ModuleNotFoundError as e:
    LOGGER.debug('Module "geographiclib" not found.', exc_info=e)
    GEOGRAPHICLIB_AVAILABLE = False
    _IMPORT_ERROR = e


class GeodesicFallback:
    def __init__(self, major_axis: float, flattening: float):
        self._major_axis = major_axis
        self._flattening = flattening

    def Inverse(self, *args, **kwargs):
        raise _IMPORT_ERROR

    @property
    def DISTANCE(self):
        raise _IMPORT_ERROR


_GEOGRAPHICLIB_WARNINGS = set()


def warn_about_missing_geographiclib(warning_message: str) -> None:
    calling_function = inspect.stack()[1][3]
    if calling_function in _GEOGRAPHICLIB_WARNINGS:
        return
    warn(warning_message, OptionalDependencyMissing)
    _GEOGRAPHICLIB_WARNINGS.add(calling_function)


if GEOGRAPHICLIB_AVAILABLE:
    Geodesic = geographiclib.geodesic.Geodesic
else:
    Geodesic = GeodesicFallback

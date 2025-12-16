from functools import wraps
from logging import getLogger
from warnings import warn

from bitsea.utilities.optional_dependencies import OptionalDependencyMissing


LOGGER = getLogger(__name__)

try:
    import numba

    _NUMBA_AVAILABLE = True
except ModuleNotFoundError as e:
    LOGGER.debug('Module "numba" not found.', exc_info=e)
    _NUMBA_AVAILABLE = False

# When the user calls a function decorated with @jit, if numba is not
# available, we should raise a warning. On the other side, we don't want to
# raise thousands of warnings for the same reason. Therefore, we save the
# names of the functions for which we already have raised a warning to avoid
# raising it again if it has already been raised.
_NUMBA_WARNINGS = set()


def jit_fallback(*_args, **_kwargs):
    def jit_decorator(func):
        f_name = func.__name__

        @wraps(func)
        def f_wrapped(*args, **kwargs):
            if f_name not in _NUMBA_WARNINGS:
                _NUMBA_WARNINGS.add(f_name)
                warn(
                    f'Function "{f_name}" is not compiled with Numba because '
                    "it is not available on the Python environment. This may "
                    "lead to significant slowdowns. Install Numba to suppress "
                    "this warning",
                    OptionalDependencyMissing,
                )
            return func(*args, **kwargs)

        return f_wrapped

    return jit_decorator


if _NUMBA_AVAILABLE:
    njit = numba.njit
    jit = numba.jit
else:
    njit = jit = jit_fallback

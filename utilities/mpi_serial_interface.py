"""
This module facilitates the creation of simple MPI scripts that can also run in
serial.

It provides an interface for writing code that can execute both serially and in
parallel using MPI. The goal is to minimize code branches by writing code once
that functions for both parallel and serial executions.

Due to the risk of encountering a SEGFAULT on certain HPC infrastructures
simply by attempting to import mpi4py, this module refrains from importing
mpi4py and defers the responsibility to the main code to handle the importing
process.

After having imported the mpi4py module (for parallel execution), you can call
the `get_mpi_communicator()` function from this module to access the standard
MPI COMM_WORLD communicator.

If the mpi4py module is not loaded, the get_mpi_communicator() function returns
a DummyCommunicator object instead. This object offers an interface similar to
MPI communicators but doesn't perform actual communication.

Thus, the following code snippet works seamlessly for both parallel and serial
executions:

>>> import mpi_serial_interface
...
... communicator = get_mpi_communicator()
... rank = communicator.Get_rank()
...
... if rank == 0:
...     # do something here
...     pass
...
... communicator.barrier()

"""

from sys import modules
import logging

LOGGER = logging.getLogger(__name__)


class DummyCommunicator:
    """
    This class provides an interface that mimics part of the interface of MPI
    communicator modules.

    An object of this class always has a rank of 0, and its size is always 1.

    It offers certain collective operations (such as `bcast`, `scatter`,
    `gather`, and `allgather`) that behave as if they were called from a
    real MPI communicator of size 1.

    It also includes a `barrier` method, which is a dummy method included
    to maintain compatibility with the MPI communicator interface.
    """
    def Get_rank(self) -> int:
        return 0

    def bcast(self, *args, **kwargs):
        return args[0]

    def scatter(self, *args, **kwargs):
        return args[0][0]

    def gather(self, *args, **kwargs):
        return [args[0]]

    def allgather(self, arg):
        return [arg]

    @property
    def size(self) -> int:
        return 1

    def Barrier(self):
        pass

    def barrier(self):
        return self.Barrier()


def get_mpi_communicator():
    """
    Return the MPI.COMM_WORLD communicator if mpi4py has been imported.
    Otherwise, it returns a DummyCommunicator object.
    """
    if 'mpi4py' in modules:
        LOGGER.debug('mpi4py has been imported. Returning a real communicator')
        from mpi4py import MPI
        return MPI.COMM_WORLD

    LOGGER.debug(
        'mpi4py has *NOT* been imported. Returning a dummy communicator'
    )
    return DummyCommunicator()


def is_a_dummy_communicator(communicator):
    return isinstance(communicator, DummyCommunicator)

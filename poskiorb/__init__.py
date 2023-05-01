from pkg_resources import DistributionNotFound, get_distribution

from . import binary, constants, utils

# set __version__ variable
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    pass  # package is not installed

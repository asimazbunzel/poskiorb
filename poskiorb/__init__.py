from pkg_resources import DistributionNotFound, get_distribution

from . import constants
from .binary import BinarySystem
from .kicks import kick_velocity_distribution, phi_distribution, theta_distribution
from .plotter import make_grid_plot, make_scatter_plot, plot_1D_distribution
from .utils import (
    P_to_a,
    a_to_f,
    a_to_P,
    binary_orbits_after_kick,
    make_grid_of_orbital_configurations,
    set_seed,
    v_orb,
)

# set __version__ variable
try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    pass  # package is not installed

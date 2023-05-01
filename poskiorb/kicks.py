"""Module that computes natal kicks  from certain distributions

Provides values for random kicks directions and strengths according to results from Kalogera, 1996
"""

import numpy as np
from scipy.stats import cauchy as lorentz
from scipy.stats import maxwell, uniform

from .constants import pi

_all_ = ["theta_distribution", "phi_distribution", "kick_velocity_distribution"]


def theta_distribution(N: int = 10000) -> np.ndarray:
    """Provides random value for tetha angle

    Parameters
    ---------
    N : `int`
       Number of trials.

    Returns
    -------
    theta : `array`
       Polar angle of kick. 0 < tetha < pi.
    """

    return np.arccos(2 * uniform.rvs(size=N, scale=1) - 1)


def phi_distribution(N: int = 10000) -> np.ndarray:
    """Function that returns random value for phi angle

    Parameters
    ---------
    N : `int`
       Number of trials.

    Returns
    -------
    phi: `array`
       Azimutal angle of kick. 0 < phi < 2*pi.
    """

    return uniform.rvs(size=N, scale=2.0 * pi)


def kick_velocity_distribution(
    distribution: str = None, N: int = 10000, sigma: float = 265.0, kwargs: dict = {}
) -> np.ndarray:
    """Random kick velocity

    Parameters
    ---------
    distribution : `str`
       Name of the probability density function (pdf) to use for natal kicks.
       Valid options are:
        - Uniform
        - Maxwell
        - Delta
        - Lorentz
        - linearly-spaced
        - log-spaced.

    N : `int`
       Number of trials.

    sigma : `float`
       Maxwellian dispersion velocity in km/s (also max value for uniform distribution).

    Returns
    -------
    w : `array`
       Velocity of natal kick in km/s.
    """

    if distribution == "Maxwell":
        return maxwell.rvs(scale=sigma, size=N)
    elif distribution == "Uniform":
        return uniform.rvs(scale=sigma, size=N)
    elif distribution == "Lorentz":
        return lorentz.rvs(scale=sigma, size=N)
    elif distribution == "Delta":
        return np.full(N, sigma)
    elif distribution == "linearly-spaced":
        return np.linspace(start=kwargs["start"], stop=kwargs["stop"], num=N, endpoint=True)
    elif distribution == "log-spaced":
        return np.logspace(
            start=np.log10(kwargs["start"]), stop=np.log10(kwargs["stop"]), num=N, base=10.0
        )
    elif distribution == "NoKicks":
        return np.zeros(N)
    elif distribution is None or distribution == "":
        raise ValueError("no probability density function given for kick velocity!")

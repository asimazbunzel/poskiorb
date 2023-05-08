"""A collection of miscellaneous utility functions
"""

from typing import Tuple, Union

import numpy as np

from .constants import Msun, Rsun, one_third, pi, standard_cgrav

__all_ = [
    "a_to_f",
    "a_to_P",
    "binary_orbits_after_kick",
    "make_grid_of_orbital_configurations",
    "P_to_a",
    "set_seed",
    "v_orb",
]


def set_seed(seed: int) -> None:
    """Set seed to reproduce runs

    Parameters
    ----------
    seed : `int`
        Seed integer number
    """
    np.random.seed(seed=seed)


def v_orb(
    r: Union[float, np.ndarray],
    m1: Union[float, np.ndarray],
    m2: Union[float, np.ndarray],
    separation: Union[float, np.ndarray],
) -> Union[float, np.ndarray]:
    """Relative orbital velocity between two objects of masses m1 and m2, with a semi-major a, at
    position r

    Parameters
    ----------
    r : `float/array`
       Position in Rsun

    m1 : `float/array`
       Mass of primary star in Msun

    m2 : `flota/array`
       Mass of secondary star in Msun

    separation : `float/array`
       Binary separation in Rsun

    Returns
    -------
    v : `float/array`
       Relative orbital velocity in km/s
    """

    # numpy arrays !
    r = np.asarray(r)
    m1 = np.asarray(m1)
    m2 = np.asarray(m2)
    separation = np.asarray(separation)

    # use cgs. we don't use, e.g., m1 *= Msun as if m1 is an int, this will throw an error
    r = r * Rsun
    m1 = m1 * Msun
    m2 = m2 * Msun
    separation = separation * Rsun

    vpow2 = standard_cgrav * (m1 + m2) * ((2e0 / r) - (1e0 / separation))
    return np.sqrt(vpow2) / 1e5


def P_to_a(
    period: Union[float, np.ndarray], m1: Union[float, np.ndarray], m2: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """Binary separation from a given orbital period

    Parameters
    ----------
    period : `float/array`
       Binary period in days

    m1 : `float/array`
       Mass of primary star in Msun

    m2 : `flota/array`
       Mass of secondary star in Msun

    Returns
    -------
    separation : `float/array`
       Binary separation in Rsun
    """

    # numpy arrays !
    period = np.asarray(period)
    m1 = np.asarray(m1)
    m2 = np.asarray(m2)

    # use cgs
    period = period * 24e0 * 3600e0
    m1 = m1 * Msun
    m2 = m2 * Msun

    separation = np.power(standard_cgrav * (m1 + m2) * np.square(period / (2 * pi)), one_third)

    return separation / Rsun


def a_to_P(
    separation: Union[float, np.ndarray], m1: Union[float, np.ndarray], m2: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """Orbital period from a given separation

    Parameters
    ----------
    separation : `float/array`
       Binary separation in Rsun

    m1: `float/array`
       Mass of primary star in Msun

    m2: `float/array`
       Mass of secondary star in Msun

    Returns
    -------
    P : `float/array`
       Binary period in days
    """

    # numpy arrays !
    separation = np.asarray(separation)
    m1 = np.asarray(m1)
    m2 = np.asarray(m2)

    # use cgs
    separation = separation * Rsun
    m1 = m1 * Msun
    m2 = m2 * Msun

    period = np.power(separation * separation * separation / (standard_cgrav * (m1 + m2)), 0.5e0)
    period = (2 * pi) * period

    return period / (24e0 * 3600e0)


def a_to_f(
    separation: Union[float, np.ndarray], m1: Union[float, np.ndarray], m2: Union[float, np.ndarray]
) -> Union[float, np.ndarray]:
    """Converts semi-major axis to orbital frequency

    Parameters
    ----------
    separation : `float/array`
       Semi-major axis

    m1 : `float/array`
       Primary mass

    m2 : `float/array`
       Secondary mass

    Returns
    -------
    f_orb : `float/array`
       Orbital frequency
    """

    # numpy arrays !
    separation = np.asarray(separation)
    m1 = np.asarray(m1)
    m2 = np.asarray(m2)

    # use cgs
    separation = separation * Rsun
    m1 = m1 * Msun
    m2 = m2 * Msun

    f_orb = np.power(standard_cgrav * (m1 + m2) / separation**3, 0.5) / (2 * pi)

    return f_orb


def binary_orbits_after_kick(
    a: float,
    m1: float,
    m2: float,
    m1_remnant_mass: float,
    w: Union[float, np.ndarray],
    theta: Union[float, np.ndarray],
    phi: Union[float, np.ndarray],
    ids: Union[int, np.ndarray],
    verbose: bool = False,
) -> Tuple[
    Union[float, np.ndarray],
    Union[float, np.ndarray],
    Union[float, np.ndarray],
    Union[float, np.ndarray],
    Union[float, np.ndarray],
    Union[float, np.ndarray],
]:
    """Function to compute binary orbital parameters after an asymmetric core-collapse

    Assuming an initial circular orbit, this function calculates the binary configuration after a
    SN explosion with an asymmetric random component. Based on the work of Kalogera (1996)

    Parameters
    ----------
    a : `float`
       Pre-SN separation in Rsun

    m1 : `float`
       Mass of collapsing star pre-SN in Msun

    m2 : `float`
       Mass of companion in Msun

    m1_remnant_mass : `float`
       Gravitational mass of compact object in Msun

    w : `float/array`
       Natal kick velocity in km/s

    theta : `float/array`
       Polar angle of kick

    phi : `float/array`
       Azimutal angle of kick velocity

    verbose : `bool`
       Flag to control additional output to user

    Returns
    -------
    a_post : `float/array`
       Post-SN separation in Rsun

    P_post : `float/array`
       Post-SN orbital period in days

    e : `float/array`
       Eccentricity of binary post-SN

    cos_i : `float/array`
       Cosine of the inclination between pre & post SN orbits

    v_sys : `float/array`
       Systemic velocity post-SN in km/s
    """

    # numpy arrays !
    a = np.asarray(a)
    m1 = np.asarray(m1)
    m2 = np.asarray(m2)
    m1_remnant_mass = np.asarray(m1_remnant_mass)
    w_ = np.asarray(w)
    theta_ = np.asarray(theta)
    phi_ = np.asarray(phi)
    ids = np.asarray(ids)

    # use cgs
    a = a * Rsun
    m1 = m1 * Msun
    m2 = m2 * Msun
    m1_remnant_mass = m1_remnant_mass * Msun
    w_ = w_ * 1e5

    if verbose:
        if w_[None].ndim > 1:
            print(f"calculating post core-collapse orbits for {len(w_)} kick(s)")
        else:
            print("calculating post core-collapse orbits for 1 kick")

    # velocity pre-SN assuming circular orbit: np.sqrt(standard_cgrav * (m1 + m2) / a)
    v_pre = v_orb(r=a / Rsun, m1=m1 / Msun, m2=m2 / Msun, separation=a / Rsun) * 1e5

    # kick velocity (w) must be projected to (x,y,z)
    wx_ = w_ * np.cos(phi_) * np.sin(theta_)
    wy_ = w_ * np.cos(theta_)
    wz_ = w_ * np.sin(phi_) * np.sin(theta_)

    # eqs. (3), (4) & (5) of Kalogera (1996)
    a_post_ = (
        standard_cgrav
        * (m1_remnant_mass + m2)
        / (
            (2 * standard_cgrav * (m1_remnant_mass + m2) / a)
            - w_**2
            - v_pre**2
            - 2 * wy_ * v_pre
        )
    )
    # pay attention to cases where e is close to 0 !
    ecc2 = 1 - (
        (wz_**2 + wy_**2 + v_pre**2 + 2 * wy_ * v_pre)
        * a**2
        / (standard_cgrav * (m1_remnant_mass + m2) * a_post_)
    )
    ecc2 = np.where(np.abs(ecc2) < 1e-8, 0, ecc2)
    e_ = np.sqrt(ecc2)

    # set np.nan to values outside boundaries, i.e., replace unbounded binaries with NaNs
    a_post_ = np.where(a_post_ > 0, a_post_, np.nan)
    e_ = np.where(e_ >= 0, e_, np.nan)
    e_ = np.where(e_ < 1, e_, np.nan)
    # patch for cases where e is very close to 0 & 1
    e_ = np.where(abs(e_) < 1e-8, 0, e_)
    e_ = np.where(abs(e_ - 1) < 1e-8, np.nan, e_)

    # only interested in bounded binaries
    if w_[None].ndim > 1:
        bounded_mask = np.isfinite(e_) & np.isfinite(a_post_)
        a_post = a_post_[bounded_mask]
        e = e_[bounded_mask]
        wx = wx_[bounded_mask]
        wy = wy_[bounded_mask]
        wz = wz_[bounded_mask]
        ids_post = ids[bounded_mask]
        w = w_[bounded_mask]
        theta = theta_[bounded_mask]
        phi = phi_[bounded_mask]

        if verbose:
            print(f"\t{len(e)} binary(ies) remain bounded " f"({len(e)/len(w_)*100:5.2f} percent)")
            print(
                f"\t{len(w_)-len(e)} binaries become "
                f"unbounded ({(len(w_)-len(e))/len(w_)*100:5.2f} percent)"
            )

    else:
        if not np.isfinite(e_) or not np.isfinite(a_post_):
            return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan

        else:
            a_post = a_post_
            e = e_
            ids_post = ids
            wx = wx_
            wy = wy_
            wz = wz_
            w = w_
            theta = theta_
            phi = phi_
            if verbose:
                print(f"\tbinary remain bounded !")

    # Inclination between pre & post SN orbits. Eq. (11) in Kalogera, 1996
    cos_i = (wy + v_pre) / np.sqrt(wz**2 + (wy + v_pre) ** 2)

    # Systemic velocity post-SN = eq.(34) x eq.(34), of Kalogera, 1996
    v_sys_2 = (
        np.power((m1_remnant_mass * wx), 2)
        + np.power((m1_remnant_mass * wz), 2)
        + np.power((m1_remnant_mass * wy - (m1 - m1_remnant_mass) * m2 / (m1 + m2) * v_pre), 2)
    )
    v_sys_2 = v_sys_2 / np.power((m1_remnant_mass + m2), 2)
    v_sys = np.sqrt(v_sys_2) / 1.0e5

    # get orbital period of bounded binaries
    P_post = a_to_P(a_post / Rsun, m1_remnant_mass / Msun, m2 / Msun)

    return a_post / Rsun, P_post, e, cos_i, v_sys, w / 1e5, theta, phi, ids_post


def make_grid_of_orbital_configurations(
    x, y, z, xrange=[0.05, 0.95], yrange=[0.0, 1.0], xnum=None, ynum=None, norm=None, verbose=False
):
    """Compute probabilities on orbital parameters post natal kick and divide it according
    to a threshold

    Parameters
    ----------
    x : `array`
       Values on the xaxis. Either `P_post` or `a_post`

    y : `array`
       Values on the yaxis. Tipically eccentricty (`e_post`)

    z : `array`
       Values on the zaxis. Use for inclination values only (`cosi`)

    xrange : `array`
       Quantiles to use as borders xrange = [xmin, xmax]

    yrange : `array`
       Quantiles to use as borders yrange = [ymin, ymax]

    xnum : `int`
       Number of grid points in the xaxis

    ynum : `int`
       Number of grid points in the yaxis

    verbose : `boolean`
       Whether to output more information to user
    """

    x = np.asarray(x)
    y = np.asarray(y)
    z = np.asarray(z)

    if verbose:
        print("making grid of orbital configurations")

    # get borders
    xmin = np.quantile(x, xrange[0])
    xmax = np.quantile(x, xrange[1])
    ymin = np.quantile(y, yrange[0])
    ymax = np.quantile(y, yrange[1])

    if norm is None:
        norm = 1

    if verbose:
        print(
            f"\tborders for xaxis (P): {xmin:.2f}, {xmax:.2f}",
            f"for quantiles [{xrange[0]:.2f},{xrange[1]:.2f}]",
        )
        print(
            f"\tborders for yaxis (e): {ymin:.2f}, {ymax:.2f}",
            f"for quantiles [{yrange[0]:.2f},{yrange[1]:.2f}]",
        )
        if norm is not None:
            print(f"\tnorm used: {norm:.6f}")

    # get (x,y)-grid
    xborders = np.logspace(np.log10(xmin), np.log10(xmax), xnum + 1)
    yborders = np.linspace(ymin, ymax, ynum + 1)

    xgrids = np.sqrt(xborders[1:] * xborders[0:-1])
    ygrids = 0.5 * (yborders[1:] + yborders[0:-1])

    # loop over each rectangle to compute its probability
    probabilities = np.zeros((xnum, ynum))
    for k, (xk, yk, zk) in enumerate(zip(x, y, z)):
        for i, xgrid in enumerate(xgrids):
            if xborders[i] <= xk < xborders[i + 1]:
                for j, ygrid in enumerate(ygrids):
                    if yborders[j] <= yk < yborders[j + 1]:
                        probabilities[i, j] += 1

    # another method: use numpy for faster computation
    zmin, zmax = -1, 1  # limits on cosi
    points = np.column_stack((x, y, z))
    probabilities2 = np.zeros((xnum, ynum))
    for kx in range(len(xborders) - 1):
        xmin = xborders[kx]
        xmax = xborders[kx + 1]
        for ky in range(len(yborders) - 1):
            ymin = yborders[ky]
            ymax = yborders[ky + 1]
            ll = np.array([xmin, ymin, zmin])  # lower-left
            ur = np.array([xmax, ymax, zmax])  # upper-right
            inidx = np.all(np.logical_and(ll <= points, points <= ur), axis=1)
            inbox = points[inidx]
            probabilities2[kx, ky] = len(inbox)
            # outbox = pts[np.logical_not(inidx)]

    return xgrids, ygrids, xborders, yborders, norm * probabilities / len(x)

'''A collection of miscellaneous utility functions
'''

from typing import Union, Tuple

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from poskiorb.constants import *

# set the default font and fontsize
plt.rc('font', family='STIXGeneral')
plt.rcParams['text.usetex'] = False
params = {'figure.figsize': (12, 8),
          'font.style': 'normal',
          'font.serif': 'DejaVu Serif',
          'font.sans-serif': 'DejaVu Sans',
          'font.monospace': 'DejaVu Sans Mono',
          'mathtext.rm': 'sans',
          'mathtext.fontset': 'stix',
          'legend.fontsize': 11,
          'axes.labelsize': 11,
          'xtick.labelsize': 10,
          'ytick.labelsize': 10,
          'xtick.top': True,
          'xtick.bottom': True,
          'xtick.direction': 'inout',
          'xtick.minor.visible': True,
          'ytick.left': True,
          'ytick.right': True,
          'ytick.direction': 'inout',
          'ytick.minor.visible': True
          }
plt.rcParams.update(params)


__all_ = ['P_to_a', 'a_to_P', 'a_to_f', 'binary_orbits_after_kick',
          'make_grid_of_orbital_configurations' 'plot_1D_distribution', 'make_scatter_plot',
          'make_grid_plot']



def P_to_a(period: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Binary separation from a known period
    
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
    a : `float/array`
       Binary separation in Rsun
    '''

    period = period * 24e0 * 3600e0  # in sec
    m1 = m1 * Msun; m2 = m2 * Msun  # in g

    separation = np.power(standard_cgrav * (m1+m2) * np.square(period/(2*pi)), one_third)

    return separation / Rsun


def a_to_P(separation: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2:Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Orbital period from a known separation
    
    Parameters
    ----------
    a : `float/array`
       Binary separation in Rsun

    m1: `float/array`
       Mass of primary star in Msun

    m2: `float/array`
       Mass of secondary star in Msun

    Returns
    -------
    P : `float/array`
       Binary period in days
    '''

    separation = separation * Rsun  # in cm
    m1 = m1 * Msun; m2 = m2 * Msun   # in g

    period = np.power(separation*separation*separation / (standard_cgrav * (m1+m2)), 0.5e0)
    period = (2*pi) * period

    return period / (24e0 * 3600e0)


def a_to_f(separation: Union[float, np.ndarray], m1: Union[float, np.ndarray],
           m2: Union[float, np.ndarray]) -> Union[float, np.ndarray]:
    '''Converts semi-major axis to orbital frequency
    
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
    '''

    separation = separation * Rsun  # in cm
    m1 = m1 * Msun; m2 = m2 * Msun  # in g

    f_orb = np.power(standard_cgrav * (m1 + m2) / separation**3, 0.5) / (2*pi)

    return f_orb


def binary_orbits_after_kick(a: float, m1:float, m2:float, m1_remnant_mass: float,
                   w: Union[float, np.ndarray], theta: Union[float, np.ndarray],
                   phi:Union[float, np.ndarray], verbose: bool=False) -> Tuple[Union[float, np.ndarray],
                   Union[float, np.ndarray], Union[float, np.ndarray], Union[float, np.ndarray],
                   Union[float, np.ndarray]]:
    '''Function to compute binary orbital parameters after an asymmetric core-collapse

    Assuming an initial circular orbit, this function calculates the binary configuration after a
    SN explosion with an asymmetric random component. Based on the work of Kalogera (1996)

    Parameters
    ----------
    a : `float`
       Pre-SN separation in Rsun.

    m1 : `float`
       Mass of collapsing star pre-SN in Msun.

    m2 : `float`
       Mass of companion in Msun.

    m1_remnant_mass : `float`
       Gravitational mass of compact object in Msun.

    w : `float/array`
       Natal kick velocity in km/s.

    theta : `float/array`
       Polar angle of kick.

    phi : `float/array`
       Azimutal angle of kick velocity.

    verbose : `bool`
       Flag to control additional output to user.

    Returns
    -------
    a_post : `float/array`
       Post-SN separation in Rsun.

    P_post : `float/array`
       Post-SN orbital period in days.

    e : `float/array`
       Eccentricity of binary post-SN.

    cos_i : `float/array`
       Cosine of the inclination between pre & post SN orbits.

    v_sys : `float/array`
       Systemic velocity post-SN in km/s
    '''

    if verbose: print('calculating post core-collapse orbits for {} kicks'.format(len(w)))

    # Input values conversion to cgs
    a = a * Rsun
    m1 = m1 * Msun
    m2 = m2 * Msun
    m1_remnant_mass = m1_remnant_mass * Msun
    w = w * 1e5

    # Velocity pre-SN
    v_pre = np.sqrt(standard_cgrav * (m1 + m2) / a)

    # Kick velocity (w) must be projected to (x,y,z)
    wx = w * np.cos(phi) * np.sin(theta)
    wy = w * np.cos(theta)
    wz = w * np.sin(phi) * np.sin(theta)

    # Eqs. (3), (4) & (5) of Kalogera (1996)
    a_post = standard_cgrav * (m1_remnant_mass + m2) / \
             (2 * standard_cgrav * (m1_remnant_mass + m2)/a - w**2 - v_pre**2 - 2 * wy * v_pre)
    e = np.sqrt(1 - (wz**2 + wy**2 + v_pre**2 + 2 * wy * v_pre) * a**2 /
             (standard_cgrav * (m1_remnant_mass + m2) * a_post))


    # only interested in bounded binaries
    bounded_mask = (a_post > 0) & (e < 1)
    a_post = a_post[bounded_mask]
    e = e[bounded_mask]
    wx = wx[bounded_mask]
    wy = wy[bounded_mask]
    wz = wz[bounded_mask]

    if verbose:
        print('\t{} binaries remain bounded ({:5.2f} percent)'.format(len(e), len(e)/len(w)*100))
        print('\t{} binaries become unbounded ({:5.2f} percent)'.format(len(w)-len(e), (len(w)-len(e))/len(w)*100))
    
    # update natal kick distro after verbose due to use of len(w)
    w = w[bounded_mask]
    theta = theta[bounded_mask]
    phi = phi[bounded_mask]

    # Inclination between pre & post SN orbits. Eq. (11) in Kalogera, 1996
    cos_i = (wy + v_pre) / np.sqrt(wz**2 + (wy + v_pre)**2)

    # Systemic velocity post-SN = eq.(34) x eq.(34), of Kalogera, 1996
    v_sys_2 = np.power((m1_remnant_mass*wx), 2) + np.power((m1_remnant_mass*wz), 2) + \
             np.power((m1_remnant_mass*wy - (m1-m1_remnant_mass)*m2/(m1+m2)*v_pre), 2)
    v_sys_2 = v_sys_2 / np.power((m1_remnant_mass+m2), 2)
    v_sys = np.sqrt(v_sys_2) / 1.e5

    # get orbital period of bounded binaries
    P_post = a_to_P(a_post/Rsun, m1_remnant_mass/Msun, m2/Msun)

    return a_post/Rsun, P_post, e, cos_i, v_sys, w/1e5, theta, phi


def make_grid_of_orbital_configurations(x, y, z, xrange=[0.05, 0.95], yrange=[0.0, 1.0],
        xnum=None, ynum=None, norm=None, verbose=False):
    '''Compute probabilities on orbital parameters post natal kick and divide it according
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
    '''

    if verbose: print('making grid of orbital configurations')

    # get borders
    xmin = np.quantile(x, xrange[0])
    xmax = np.quantile(x, xrange[1])
    ymin = np.quantile(y, yrange[0])
    ymax = np.quantile(y, yrange[1])

    if norm is None: norm = 1

    if verbose:
        print('\tborders for xaxis (P): {:.2f}, {:.2f}'.format(xmin, xmax),
              'for quantiles [{:.2f},{:.2f}]'.format(xrange[0], xrange[1]))
        print('\tborders for yaxis (e): {:.2f}, {:.2f}'.format(ymin, ymax),
              'for quantiles [{:.2f},{:.2f}]'.format(yrange[0], yrange[1]))
        if norm is not None:
            print('\tnorm used: {:.6f}'.format(norm))

    # get (x,y)-grid
    xborders = np.logspace(np.log10(xmin), np.log10(xmax), xnum+1)
    yborders = np.linspace(ymin, ymax, ynum+1)

    xgrids = np.sqrt(xborders[1:] * xborders[0:-1])
    ygrids = 0.5*(yborders[1:] + yborders[0:-1])

    # loop over each rectangle to compute its probability
    probabilities = np.zeros((xnum, ynum))
    for k, (xk, yk, zk) in enumerate(zip(x, y, z)):
        for i, xgrid in enumerate(xgrids):
            if xborders[i] <= xk < xborders[i+1]:
                for j, ygrid in enumerate(ygrids):
                    if yborders[j] <= yk < yborders[j+1]:
                        probabilities[i,j] += 1

    # another method: use numpy for faster computation
    zmin, zmax = -1, 1  # limits on cosi
    points = np.column_stack((x, y, z))
    probabilities2 = np.zeros((xnum, ynum))
    for kx in range(len(xborders)-1):
        xmin = xborders[kx]
        xmax = xborders[kx+1]
        for ky in range(len(yborders)-1):
            ymin = yborders[ky]
            ymax = yborders[ky+1]
            ll = np.array([xmin, ymin, zmin])  # lower-left
            ur = np.array([xmax, ymax, zmax])  # upper-right
            inidx = np.all(np.logical_and(ll <= points, points <= ur), axis=1)
            inbox = points[inidx]
            probabilities2[kx,ky] = len(inbox)
            # outbox = pts[np.logical_not(inidx)]

    return xgrids, ygrids, xborders, yborders, norm*probabilities/len(x)
    


def plot_1D_distribution(x, weights=None, disttype='hist', fig=None, ax=None, xlabel=None,
            ylabel=None, xlim=None, ylim=None, color=None, show=True, **kwargs):
    '''plot a 1D distribution of ``x``

    This function is a wrapper for :func:`matplotlib.pyplot.hist`,
    :func:`seaborn.kdeplot` and :func:`seaborn.ecdfplot`.

    Copied from the excellent repo `The LISA Evolution and Gravitational Wave ORbit Kit`
    All credits goes to authors: https://github.com/katiebreivik/LEGWORK (visualization.py module)
    
    Parameters
    ----------
    x : `float/int array`
        Variable to plot, should be a 1D array

    weights : `float/int array`
        Weights for each variable in ``x``, must have the same shape

    disttype : `{{ 'hist', 'kde', 'ecdf' }}`
        Which type of distribution plot to use

    fig: `matplotlib Figure`
        A figure on which to plot the distribution. Both `ax` and `fig` must be
        supplied for either to be used

    ax: `matplotlib Axis`
        An axis on which to plot the distribution. Both `ax` and `fig` must be
        supplied for either to be used

    xlabel : `string`
        Label for the x axis, passed to Axes.set_xlabel()

    ylabel : `string`
        Label for the y axis, passed to Axes.set_ylabel()

    xlim : `tuple`
        Lower and upper limits for the x axis, passed to Axes.set_xlim()

    ylim : `tuple`
        Lower and upper limits for the y axis, passed to Axes.set_ylim()

    color : `string or tuple`
        Colour to use for the plot, see
        https://matplotlib.org/tutorials/colors/colors.html for details on how
        to specify a colour

    show : `boolean`
        Whether to immediately show the plot or only return the Figure and Axis

    **kwargs : `(if disttype=='hist')`
        Include values for any of `bins, range, density, cumulative, bottom,
        histtype, align, orientation, rwidth, log, label`. See
        :func:`matplotlib.pyplot.hist` for more details.

    **kwargs : `(if disttype=='kde')`
        Include values for any of `gridsize, cut, clip, legend, cumulative,
        bw_method, bw_adjust, log_scale, fill, label, linewidth, linestyle`.
        See :func:`seaborn.kdeplot` for more details.

    **kwargs : `(if disttype=='ecdf')`
        Include values for any of `stat, complementary, log_scale, legend,
        label, linewidth, linestyle`. See :func:`seaborn.edcfplot`
        for more details.

    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted

    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    '''

    # create new figure and axes is either weren't provided
    if fig is None or ax is None: fig, ax = plt.subplots()

    # possible kwargs for matplotlib.hist
    hist_args = {'bins': 'auto', 'range': None, 'density': True,
                 'cumulative': False, 'bottom': None, 'histtype': 'bar', 'align': 'mid',
                 'orientation': 'vertical', 'rwidth': None, 'log': False, 'label': None}

    # possible kwargs for seaborn.kdeplot
    kde_args = {'gridsize': 200, 'cut': 3, 'clip': None, 'legend': True,
                'cumulative': False, 'bw_method': 'scott', 'bw_adjust': 1, 'log_scale': None,
                'fill': None, 'label': None, 'linewidth': None, 'linestyle': None}

    # possible kwargs for seaborn.ecdfplot
    ecdf_args = {'stat': 'proportion', 'complementary': False, 'log_scale': None,
                 'legend': True, 'label': None, 'linewidth': None, 'linestyle': None}

    # set which ones we are using for this plot
    plot_args = hist_args if disttype == 'hist' else kde_args if disttype == 'kde' else ecdf_args

    # update the values with those supplied
    for key, value in kwargs.items():
        if key in plot_args:
            plot_args[key] = value
        else:
            # warn user if they give an invalid kwarg
            print('Warning: keyword argument `{}`'.format(key),
                  'not recognised for disttype `{}`'.format(disttype), 'and will be ignored')

    # create whichever plot was requested
    if disttype == 'hist':
        ax.hist(x, weights=weights, color=color, **plot_args)
    elif disttype == 'kde':
        sns.kdeplot(x=x, weights=weights, ax=ax, color=color, **plot_args)
    elif disttype == 'ecdf':
        sns.ecdfplot(x=x, weights=weights, ax=ax, color=color, **plot_args)

    # update axis labels
    if xlabel is not None: ax.set_xlabel(xlabel)
    if ylabel is not None: ax.set_ylabel(ylabel)

    # update axis limits
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show: plt.show()

    # return the figure and axis for further plotting
    return fig, ax


def make_scatter_plot(x, y, fig=None, ax=None, xlabel=None, ylabel=None, xlim=None, ylim=None,
        s=None, color=None, marker=None, show=True, xlogscale=True, ylogscale=False, **kwargs):
    '''make scatter plot
    
    This function is a wrapper for :func:`matplotlib.pyplot.scatter`,
    
    Parameters
    ---------
    x : `float/int array`
        Variable to plot on xaxis, should be a 1D array
    
    y : `float/int array`
        Variable to plot on yaxis, should be a 1D array

    fig: `matplotlib Figure`
        A figure on which to plot the distribution. Both `ax` and `fig` must be
        supplied for either to be used

    ax: `matplotlib Axis`
        An axis on which to plot the distribution. Both `ax` and `fig` must be
        supplied for either to be used

    xlabel : `string`
        Label for the x axis, passed to Axes.set_xlabel()

    ylabel : `string`
        Label for the y axis, passed to Axes.set_ylabel()

    xlim : `tuple`
        Lower and upper limits for the x axis, passed to Axes.set_xlim()

    ylim : `tuple`
        Lower and upper limits for the y axis, passed to Axes.set_ylim()

    s : `integer`
        Size of the scatter points

    color : `string or tuple`
        Colour to use for the plot, see
        https://matplotlib.org/tutorials/colors/colors.html for details on how
        to specify a colour

    marker : `string`
        The marker style. See matplotlib.markers for more information about marker styles.

    show : `boolean`
        Whether to immediately show the plot or only return the Figure and Axis

    xlogscale : `boolean`
        Whether to use log scale on the xaxis
    
    ylogscale : `boolean`
        Whether to use log scale on the yaxis

    **kwargs : `(if disttype=='hist')`
        Include values for any of the rest of arguments to pass to matplotlib scatter.
        See matplotlib scatter doc for more info.ç

    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted
    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    '''

    # create new figure and axes is either weren't provided
    if fig is None or ax is None: fig, ax = plt.subplots()

    if xlogscale: ax.set_xscale('log')
    if ylogscale: ax.set_yscale('log')
    
    # create scatter plot
    ax.scatter(x, y, s=s, c=color, marker=marker, **kwargs)

    # update axis labels
    if xlabel is not None: ax.set_xlabel(xlabel)
    if ylabel is not None: ax.set_ylabel(ylabel)

    # update axis limits
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show: plt.show()

    # return the figure and axis for further plotting
    return fig, ax


def make_grid_plot(xgrid, ygrid, xborders, yborders, annotations, fig, ax, xlim=None, ylim=None,
            show=True):
    '''Fill 2D binary parameter space after kick with a grid based on some probability

    Parameters
    ----------
    xgrid : `array`
        X coordinate of center of each rectangle in the grid
    
    ygrid : `array`
        Y coordinate of center of each rectangle in the grid
        
    xborders : `array`
        X coordinate of border of each rectangle in the grid. Sort of bins in xaxis
    
    yborders : `array`
        Y coordinate of border of each rectangle in the grid. Sort of bins in yaxis

    annotations : `array`
        What to annotate in (xgrid, ygrid). Tipically np.log10(probability)

    fig: `matplotlib Figure`
        A figure on which to plot the distribution. Both `ax` and `fig` must be
        supplied for either to be used

    ax: `matplotlib Axis`
        An axis on which to plot the distribution. Both `ax` and `fig` must be
        supplied for either to be used
    
    xlim : `tuple`
        Lower and upper limits for the x axis, passed to Axes.set_xlim()

    ylim : `tuple`
        Lower and upper limits for the y axis, passed to Axes.set_ylim()

    show : `boolean`
        Whether to immediately show the plot or only return the Figure and Axis

    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted
    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    '''
    
    for xb in xborders: ax.axvline(xb, ls='--', color='gray', zorder=99, alpha=0.45)
    for yb in yborders: ax.axhline(yb, ls='--', color='gray', zorder=99, alpha=0.45)

    for i, xg in enumerate(xgrid):
        for j, yg in enumerate(ygrid):
            prob = annotations[i,j]
            if prob > 0.01:
                ax.text(xg, yg, '{:.1f}'.format(np.log10(prob)), c='black', ha='center',
                        va='center', fontsize=7,
                        bbox=dict(facecolor='white', edgecolor='C0', alpha=0.75, pad=2))

    # update axis limits
    if xlim is not None: ax.set_xlim(xlim)
    if ylim is not None: ax.set_ylim(ylim)
    
    # immediately show the plot if requested
    if show: plt.show()

    return fig, ax

"""A collection of plotting functions
"""

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

__all__ = [
    "make_grid_plot",
    "make_scatter_plot",
    "plot_1D_distribution",
]

# set the default font and fontsize
plt.rc("font", family="STIXGeneral")
plt.rcParams["text.usetex"] = False
params = {
    "figure.figsize": (12, 8),
    "font.style": "normal",
    "font.serif": "DejaVu Serif",
    "font.sans-serif": "DejaVu Sans",
    "font.monospace": "DejaVu Sans Mono",
    "mathtext.rm": "sans",
    "mathtext.fontset": "stix",
    "legend.fontsize": 11,
    "axes.labelsize": 11,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "xtick.top": True,
    "xtick.bottom": True,
    "xtick.direction": "inout",
    "xtick.minor.visible": True,
    "ytick.left": True,
    "ytick.right": True,
    "ytick.direction": "inout",
    "ytick.minor.visible": True,
}
plt.rcParams.update(params)


def plot_1D_distribution(
    x,
    weights=None,
    disttype="hist",
    fig=None,
    ax=None,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    color=None,
    show=True,
    **kwargs,
):
    """plot a 1D distribution of ``x``

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

    disttype : `{{ "hist", "kde", "ecdf" }}`
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

    **kwargs : `(if disttype=="hist")`
        Include values for any of `bins, range, density, cumulative, bottom,
        histtype, align, orientation, rwidth, log, label`. See
        :func:`matplotlib.pyplot.hist` for more details.

    **kwargs : `(if disttype=="kde")`
        Include values for any of `gridsize, cut, clip, legend, cumulative,
        bw_method, bw_adjust, log_scale, fill, label, linewidth, linestyle`.
        See :func:`seaborn.kdeplot` for more details.

    **kwargs : `(if disttype=="ecdf")`
        Include values for any of `stat, complementary, log_scale, legend,
        label, linewidth, linestyle`. See :func:`seaborn.edcfplot`
        for more details.

    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted

    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    """

    # create new figure and axes is either weren"t provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    # possible kwargs for matplotlib.hist
    hist_args = {
        "bins": "auto",
        "range": None,
        "density": True,
        "cumulative": False,
        "bottom": None,
        "histtype": "bar",
        "align": "mid",
        "orientation": "vertical",
        "rwidth": None,
        "log": False,
        "label": None,
    }

    # possible kwargs for seaborn.kdeplot
    kde_args = {
        "gridsize": 200,
        "cut": 3,
        "clip": None,
        "legend": True,
        "cumulative": False,
        "bw_method": "scott",
        "bw_adjust": 1,
        "log_scale": None,
        "fill": None,
        "label": None,
        "linewidth": None,
        "linestyle": None,
    }

    # possible kwargs for seaborn.ecdfplot
    ecdf_args = {
        "stat": "proportion",
        "complementary": False,
        "log_scale": None,
        "legend": True,
        "label": None,
        "linewidth": None,
        "linestyle": None,
    }

    # set which ones we are using for this plot
    plot_args = hist_args if disttype == "hist" else kde_args if disttype == "kde" else ecdf_args

    # update the values with those supplied
    for key, value in kwargs.items():
        if key in plot_args:
            plot_args[key] = value
        else:
            # warn user if they give an invalid kwarg
            print(
                f"Warning: keyword argument `{key}`",
                f"not recognised for disttype `{disttype}`",
                "and will be ignored",
            )

    # create whichever plot was requested
    if disttype == "hist":
        ax.hist(x, weights=weights, color=color, **plot_args)
    elif disttype == "kde":
        sns.kdeplot(x=x, weights=weights, ax=ax, color=color, **plot_args)
    elif disttype == "ecdf":
        sns.ecdfplot(x=x, weights=weights, ax=ax, color=color, **plot_args)

    # update axis labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    # update axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show:
        plt.show()

    # return the figure and axis for further plotting
    return fig, ax


def make_scatter_plot(
    x,
    y,
    fig=None,
    ax=None,
    xlabel=None,
    ylabel=None,
    xlim=None,
    ylim=None,
    s=None,
    color=None,
    marker=None,
    show=True,
    xlogscale=True,
    ylogscale=False,
    **kwargs,
):
    """make scatter plot

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

    **kwargs : `(if disttype=="hist")`
        Include values for any of the rest of arguments to pass to matplotlib scatter.
        See matplotlib scatter doc for more info.ç

    Returns
    -------
    fig : `matplotlib Figure`
        The figure on which the distribution is plotted
    ax : `matplotlib Axis`
        The axis on which the distribution is plotted
    """

    # create new figure and axes is either weren"t provided
    if fig is None or ax is None:
        fig, ax = plt.subplots()

    if xlogscale:
        ax.set_xscale("log")
    if ylogscale:
        ax.set_yscale("log")

    # create scatter plot
    ax.scatter(x, y, s=s, c=color, marker=marker, **kwargs)

    # update axis labels
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

    # update axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show:
        plt.show()

    # return the figure and axis for further plotting
    return fig, ax


def make_grid_plot(
    xgrid, ygrid, xborders, yborders, min_z, annotations, fig, ax, xlim=None, ylim=None, show=True
):
    """Fill 2D binary parameter space after kick with a grid based on some probability

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

    min_z : `float`
        Minimum Z value above which `annotations` will be added to the plot

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
    """

    for xb in xborders:
        ax.axvline(xb, ls="--", color="gray", zorder=99, alpha=0.45)
    for yb in yborders:
        ax.axhline(yb, ls="--", color="gray", zorder=99, alpha=0.45)

    for i, xg in enumerate(xgrid):
        for j, yg in enumerate(ygrid):
            prob = annotations[i, j]
            if prob > min_z:
                ax.text(
                    xg,
                    yg,
                    f"{np.log10(prob):.1f}",
                    c="black",
                    ha="center",
                    va="center",
                    fontsize=7,
                    bbox=dict(facecolor="white", edgecolor="C0", alpha=0.75, pad=2),
                )

    # update axis limits
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)

    # immediately show the plot if requested
    if show:
        plt.show()

    return fig, ax

'''A collection of classes for binary systems at core-collapse stage
'''

from typing import Any, Callable, Tuple, Union

import numpy as np

from biaswise import utils
from biaswise import kicks

__all__ = ['BinarySytem']



class BinarySystem(object):
    '''Class containing the status of a binary system at core-collapse

    It contains information on the binary such as orbital period, separation, eccentricty and
    masses, but also holds information of the collapsing star (core mass, fallback fraction
    and remnant mass). All this is then used to populate the post-collapse binary parameters
    2D grid of period/separation and eccentricity assuming an instantaneous asymmetric kick
    distribution

    Parameters
    ----------
    m1 : `float`
       Mass of the collapsing star in Msun units.
    
    m1_core_mass : `float`
       Mass of the carbon-oxygen core of the collapsing star in Msun units.

    m1_fallback_fraction : `float`
       Fraction of mass that fallsback during core-collapse to the compact object.

    m1_remnant_mass : `float`
       Remnant mass of the compact object left after core-collapse in Msun units.

    m2 : `float`
       Mass of the companion star in Msun units.

    P : `float`
       Binary orbital period just before collapse (either `P` or `a` must be supplied) in days
       units.

    a : `float`
       Binary separation just before collapse (either `a` or `P` must be supplied) in Rsun units.

    Attributes
    ----------
    f_orb : `float`
       Orbital frequency just before collapse in units of frequency
    '''

    def __init__(self, m1: float, m1_core_mass: float, m1_fallback_fraction: float,
            m1_remnant_mass: float, m2: float, P: float=None, a: float=None) -> None:

        if P is None and a is None: raise ValueError('either `P` or `a` must be specified')

        # a is not given, compute using Kepler
        if a is None: a = utils.P_to_a(P, m1, m2)
        # P is not given, compute using Kepler
        if P is None: P = utils.a_to_P(a, m1, m2)

        # compute orbital frequency with separation
        f_orb = utils.a_to_f(a, m1, m2)

        self.m1 = float(m1)
        self.m1_core_mass = float(m1_core_mass)
        self.m1_fallback_fraction = float(m1_fallback_fraction)
        self.m1_remnant_mass = float(m1_remnant_mass)
        self.m2 = float(m2)
        self.P = float(P)
        self.a = float(a)
        self.f_orb = float(f_orb)


    def set_natal_kick_distribution(self, n_trials: int=1, distribution_id: str=None,
            seed: int=22, kick_sigma: float=265e0, min_v_kick: float=0e0, max_v_kick: float=1e99,
            kick_scaling: Callable[[Any],Any]=lambda x: x) -> None:
        '''Set distribution in the magnitude of the natal kick
        
        The angle distribution of the asymmetric natal kick is assumed to be isotropic

        Parameters
        ----------
        n_trials : `integer`
           Number of trials to perform draws from a distribution.

        distribution_id : `string`
           Name of the distribution of natal kicks.

        seed : `integer`
           Seed to start random number generator.

        kick_sigma : `float`
           Dispersion velocity for a Maxwellian distribution in km/s.

        min_v_kick : `float`
           Min kick velocity in km/s (only used for linearly-spaced and log-spaced cases).
        
        max_v_kick : `float`
           Max kick velocity in km/s (only used for linearly-spaced and log-spaced cases).

        kick_scaling : `function`
           Some function to apply to natal kick strength as a scaling factor. Default is
           to leave kick strength as it is.
        '''
        
        # n_trials must be integer > 0
        if n_trials <= 0: raise ValueError('`n_trials` needs to be a positive integer')

        if distribution_id is None: raise ValueError('`distribution_id` needs to be supplied')

        possible_distribution_ids = ('Maxwell', 'Uniform', 'Lorentz', 'Delta', 'linearly-spaced',
                                     'log-spaced', 'NoKicks')
        if distribution_id not in possible_distribution_ids:
            msg = 'cannot find `distribution_id` in list: '
            for id in possible_distribution_ids: msg += '{}, '.format(id)
            raise ValueError(msg[:-2])

        if kick_sigma < 0: raise ValueError('`kick_sigma must be positive`')

        self.natal_kick_info = {
                'n_trials': n_trials,
                'seed': seed,
                'kick_distribution': distribution_id,
                'kick_sigma': kick_sigma,
                'min_v_kick': min_v_kick,
                'max_v_kick': max_v_kick,
                'kick_scaling': kick_scaling
                }


    def get_natal_kick_distribution(self) -> None:
        '''Compute random asymmetric natal kicks from a given distribution

        The orientation of the kick, controlled by the (theta, phi) variables, is assumed
        to be isotropically distributed
        '''

        # call to random kick distributions
        theta = kicks.theta_distribution(N=self.natal_kick_info['n_trials'])
        
        phi = kicks.phi_distribution(N=self.natal_kick_info['n_trials'])

        w = kicks.kick_velocity_distribution(distribution=self.natal_kick_info['kick_distribution'],
                N=self.natal_kick_info['n_trials'], sigma=self.natal_kick_info['kick_sigma'],
                kwargs={self.natal_kick_info['min_v_kick'], self.natal_kick_info['max_v_kick']})

        # apply kick scaling function to randomly drawn `w`
        w = self.natal_kick_info['kick_scaling'](w)

        # if fallback_fraction = 1, no asymmetric kick, only Blauuw kick
        if np.all(w == 0):
            w = np.array([0])
            theta = np.array([0])
            phi = np.array([0])

        self.theta = theta
        self.phi = phi
        self.w = w


    def plot_kick_distribution(self, xattr: str, **kwargs) -> None:
        '''Plot utility for natal kick distributions (w, theta, phi)
        
        Parameters
        ----------
        xattr : `string`
           Name of variable to plot on xaxis

        **kwargs : `various`
           Contains additional stuff to make histogram/KDE of kick distributions
        '''

        axis_map = {'w': self.w, 'theta': self.theta, 'phi': self.phi, 'w_post': self.w_post,
                    'theta_post': self.theta_post, 'phi_post': self.phi_post, 'cos_i': self.cosi}

        labels = {'w': '$v_{\\rm kick}$ [km s$^{-1}$]', 'theta': '$\\theta$', 'phi': '$\\phi$',
                  'w_post': '$v_{\\rm kick}$ [km s$^{-1}$]', 'theta_post': '$\\theta$',
                  'phi_post': '$\\phi$', 'cos_i': '$\\cos\\,(i)$'}

        if xattr not in axis_map.keys():
            msg = '`xattr` must be one of: ' + ', '.join(['`{}`'.format(k) for k in
                list(axis_map.keys())])
            raise ValueError(msg)

        x = axis_map[xattr]
        if x is None: raise ValueError('`{}` cannot be None'.format(xattr))
        if len(x) == 1: raise ValueError('cannot plot single point into a distribution')

        # pass xlabel in the kwargs
        if 'xlabel' not in kwargs.keys(): kwargs['xlabel'] = labels[xattr]

        utils.plot_1D_distribution(x=x, **kwargs)


    def get_orbital_distribution(self, verbose: bool=False) -> None:
        '''Evaluate orbital parameter distribution based on a distribution of natal kicks

        Parameters
        ----------
        verbose : `boolean`
           Flag to control the output of additional info to user
        '''

        # compute new binary orbital parameters
        a_post, P_post, e_post, cosi, v_sys, w, theta, phi = utils.binary_orbits_after_kick(
                a=self.a, m1=self.m1, m2=self.m2, m1_remnant_mass=self.m1_remnant_mass, w=self.w,
                theta=self.theta, phi=self.phi, verbose=verbose)

        # update object with binaries bounded after asymmetric kick
        self.a_post = a_post
        self.P_post = P_post
        self.e_post = e_post
        self.cosi = cosi
        self.v_sys = v_sys
        self.w_post = w
        self.theta_post = theta
        self.phi_post = phi


    def plot_post_kick_orbital_configurations(self, xattr: str, yattr: str,
            **kwargs) -> Tuple[Any, Any]:
        '''TBD
        '''

        axis_map = {'P': self.P_post, 'e': self.e_post, 'a': self.a_post}

        labels = {'P': '$P_{\\rm orb}$ [d]', 'e': 'eccentricity', 'a': '$a$ [R$_\\odot$]'}

        for attr in (xattr, yattr):
            if attr not in axis_map.keys():
                msg = '`xattr` and `yattr` must be one of: ' + ', '.join(['`{}`'.format(k) for k in
                    list(axis_map.keys())])
                raise ValueError(msg)

        x = axis_map[xattr]
        if x is None: raise ValueError('`{}` cannot be None'.format(xattr))
        if len(x) == 1: raise ValueError('cannot plot single point')

        # pass xlabel in the kwargs
        if 'xlabel' not in kwargs.keys(): kwargs['xlabel'] = labels[xattr]

        y = axis_map[yattr]
        if y is None: raise ValueError('`{}` cannot be None'.format(yattr))
        if len(y) == 1: raise ValueError('cannot plot single point')

        # pass xlabel in the kwargs
        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = labels[yattr]

        return utils.make_scatter_plot(x, y, **kwargs)


    def get_post_kick_grid(self, xnum: int=10, ynum: int=10,
            xquantiles: Tuple[float, float]=[0.05, 0.95],
            yquantiles: Tuple[float, float]=[0.00, 1.00],
            min_prob: float=0.01, use_unbounded_for_norm: bool=False,
            verbose: bool=False) -> None:
        '''Compute a 2D grid of orbital configurations of binaries surviving asymmetric kicks

        Based on two arrays of the same length, it divides the 2D plane into rectangular grids,
        computes the probability of each rectangle (MonteCarlo approach, a simple summation)
        and returns the grid and the probabilities.
        
        Parameters
        ----------
        xnum : `integer`
           Number of bins for the xattr

        ynum : `integer`
           Number of bins for the yattr

        xquantiles : `array`
           Compute bins on xaxis between (xquantiles[0], xquantiles[1])
           
        yquantiles : `array`
           Compute bins on yaxis between (yquantiles[0], yquantiles[1])

        min_prob : `float`
           Minimum probability value to consider it a possible candidate rectangle to further
           explore with a detailed evolutionary simulation. Default: 1% (min_prob=0.01)

        use_unbounded_for_norm : `boolean`
           Whether to consider the complete sample of binaries when evaluating probabilities

        verbose : `boolean`
           Output more information for user
        '''

        axis_map = {'P': self.P_post, 'e': self.e_post, 'a': self.a_post}
        
        x = self.P_post
        if x is None: raise ValueError('`P_post` cannot be None')
        y = self.e_post
        if y is None: raise ValueError('`e_post` cannot be None')

        # check that quantiles are OK
        if xquantiles[0] > xquantiles[1]:
            raise ValueError('xquantiles must be in increasing order')
        if yquantiles[0] > yquantiles[1]:
            raise ValueError('yquantiles must be in increasing order')

        # also, make sure xnum and ynum are positive integer
        if xnum < 0: raise ValueError('xnum must be a positive integer')
        if ynum < 0: raise ValueError('ynum must be a positive integer')

        norm = None
        if use_unbounded_for_norm:
            norm = len(self.w_post) / len(self.w)

        xgrid, ygrid, xborders, yborders, probs = utils.make_grid_of_orbital_configurations(x, y,
                xquantiles, yquantiles, xnum, ynum, norm, verbose)

        self._P_post_grid = xgrid
        self._a_post_grid = utils.P_to_a(xgrid, self.m1_remnant_mass, self.m2)
        self._P_post_borders = xborders
        self._a_post_borders = utils.P_to_a(xborders, self.m1_remnant_mass, self.m2)
        self._e_post_grid = ygrid
        self._e_post_borders = yborders
        self._post_probabilities = probs

        # now get probability limit for post kick binary parameters, using meshgrid to vectorize
        # computation
        self._min_prob = min_prob
        X, Y = np.meshgrid(self._P_post_grid, self._e_post_grid)
        Z = self._post_probabilities

        # put np.nan where probability is below threshold
        ZZ = np.where(Z.T > min_prob, Z.T, np.nan)
        XX = np.where(Z.T > min_prob, X, np.nan)
        YY = np.where(Z.T > min_prob, Y, np.nan)

        self.P_post_grid = XX[~np.isnan(XX)]
        self.a_post_grid = utils.P_to_a(self.P_post_grid, self.m1_remnant_mass, self.m2)
        self.e_post_grid = YY[~np.isnan(YY)]
        self.prob_grid = ZZ[~np.isnan(ZZ)].ravel()


    def show_post_kick_with_grid(self, xattr: str, yattr: str, min_prob: float=0.01, **kwargs):
        '''Plot grid of post kick orbital configurations
        '''

        scatter_map = {'P': self.P_post, 'e': self.e_post, 'a': self.a_post}

        labels = {'P': '$P_{\\rm orb}$ [d]', 'e': 'eccentricity', 'a': '$a$ [R$_\\odot$]'}

        for attr in (xattr, yattr):
            if attr not in scatter_map.keys():
                msg = '`xattr` and `yattr` must be one of: ' + ', '.join(['`{}`'.format(k) for k in
                    list(axis_map.keys())])
                raise ValueError(msg)

        x = scatter_map[xattr]
        if x is None: raise ValueError('`{}` cannot be None'.format(xattr))
        if len(x) == 1: raise ValueError('cannot plot single point')

        # pass xlabel in the kwargs
        if 'xlabel' not in kwargs.keys(): kwargs['xlabel'] = labels[xattr]

        y = scatter_map[yattr]
        if y is None: raise ValueError('`{}` cannot be None'.format(yattr))
        if len(y) == 1: raise ValueError('cannot plot single point')

        # pass xlabel in the kwargs
        if 'ylabel' not in kwargs.keys(): kwargs['ylabel'] = labels[yattr]

        fig, ax = utils.make_scatter_plot(x, y, show=False, **kwargs)

        grid_map = {'P': self._P_post_grid, 'a': self._a_post_grid, 'e': self._e_post_grid}
        borders_map = {'P': self._P_post_borders, 'a': self._a_post_borders,
                       'e': self._e_post_borders}

        xgrid = grid_map[xattr]
        ygrid = grid_map[yattr]

        xborders = borders_map[xattr]
        yborders = borders_map[yattr]

        xlims = [0.9*xborders[0], 1.1*xborders[-1]]
        ylims = [yborders[0], yborders[-1]]

        return utils.make_grid_plot(xgrid, ygrid, xborders, yborders, self._post_probabilities,
                fig, ax, xlims, ylims)


    def save_target_grid(self, fname: str='grid.data') -> None:
        '''Save orbital parameters conforming 2D grid of Porb (or separation) and e

        It contains a header with the following columns:

        # asymmetric natal-kick id    orbital period [days]   separation [Rsun]   eccentricity

        It saves some info on the distribution used for the natal kicks with important values.
        In addition, the probability set to define grid is also present
        
        Parameters
        ----------
        fname : `str`
            Name of file where data will be saved.
        '''

        def format_string(value):
            if isinstance(value, str) or isinstance(value, int):
                return '{:>19}'.format(value)
            else:
                return '{:>19E}'.format(value)

        msg = '# Target grid of orbital parameters\n'
        msg += '# Binary at core-collapse\n'
        msg += '#'
        header_names = ['m1 [Msun]', 'm2 [Msun]', 'P [days]', 'a [Rsun]', 'm1_core [Msun]',
                        'm1_fb', 'm1_remnant [Msun]']
        for name in header_names: msg += '{}'.format(format_string(name))
        msg += '\n'
        msg += '#'
        header_values = [self.m1, self.m2, self.P, self.a, self.m1_core_mass,
                         self.m1_fallback_fraction, self.m1_remnant_mass]
        for value in header_values: msg += '{}'.format(format_string(value))
        msg += '\n'

        msg += '# Asymmetric natal kick parameters\n'
        msg += '#'
        header_names = ['distribution', 'sigma', 'min_w', 'max_w', 'N_trials', 'min_prob']
        for name in header_names: msg += '{}'.format(format_string(name))
        msg += '\n'
        msg += '#'
        header_values = [self.natal_kick_info['kick_distribution'],
                         self.natal_kick_info['kick_sigma'],
                         self.natal_kick_info['min_v_kick'],
                         self.natal_kick_info['max_v_kick'],
                         self.natal_kick_info['n_trials'],
                         self._min_prob]
        for value in header_values: msg += '{}'.format(format_string(value))
        msg += '\n\n'

        column_names = ['natal kick id', 'period [days]', 'separation [Rsun]', 'eccentricity',
                        'probability']
        for name in column_names: msg += '{}'.format(format_string(name))

        # natal kick id should have a length according to the number of kicks
        n = len(str(len(self.P_post_grid)))
        string = '{:0' + str(n) + 'd}'
        id_names = [string.format(k+1) for k in range(len(self.P_post_grid))]

        msg += '\n'
        for k in range(len(self.P_post_grid)):
            msg += '{}{}{}{}{}\n'.format(format_string(id_names[k]),
                    format_string(self.P_post_grid[k]), format_string(self.a_post_grid[k]),
                    format_string(self.e_post_grid[k]), format_string(self.prob_grid[k]))

        with open(fname, 'w') as f:
            f.write(msg)


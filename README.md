# BIASWISE

Compute binary parameters after a core-collapse with an asymmetric kick

## Requirements

This package has three dependencies: `numpy`, `matplotlib` and `seaborn`. They are installed (or updated)
when this package is installed (see below)

## Installation

Clone (using `git clone <repo-name>`) or download this repository. Then `cd` into it and run `pip install .`
This will install dependencies as well as the module `biaswise`

## Usage

In order to use it, you just simply need to import the module (`import biaswise`) when using the interactive
python shell, via a script or a jupyter notebook

## Example

In an interactive python shell: 

```
$ python
```

Load module:

```
import biaswise
```

Create a binary system class with a star of 12.81 Msun collapsing into a compact object of 8 Msun, a
companion of 7 Msun with an orbital period pre-collapse of 5 days:

```
binary = biaswise.binary.BinarySystem(m1=12.81, m1_core_mass=9, m1_remnant_mass=8, m1_fallback_fraction=0.50, m2=7, P=5)
```

Configurate the distribution of natal kicks to be randomly computed. In this case, 50$\,$000 kicks from a
Maxwellian distribution in its strength while scaling this magnitude by the amount of fallback
($v_{\rm kick} = v_{\rm kick} * (1 - f_{\rm fb})$):

```
binary.set_natal_kick_distribution(n_trials=50000, distribution_id='Maxwell',
                                   kick_scaling=lambda x: (1-binary.m1_fallback_fraction)*x)
```


Then, compute this asymmetric kicks:

```
binary.get_natal_kick_distribution()
```

After having all the kicks, compute the post-kick orbital configurations (only those which remain bound):

```
binary.get_orbital_distribution(verbose=True)
```

Once the orbital parameters are known, divide the parameter space of separation (or orbital period) and
eccentricity according to the probability of finding binaries in different 2D bins: 

```
binary.get_post_kick_grid(use_unbounded_for_norm=True, verbose=True)
```

If needed, show this in a matplotlib figure:

```
binary.show_post_kick_with_grid(xattr='P', yattr='e', s=1)
```

And save the grid of post-orbital parameters together with information on the probabilities for each bin:

```
binary.save_target_grid(fname='grid.data')
```

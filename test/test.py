
import poskiorb


binary = poskiorb.binary.BinarySystem(m1=12.81, m1_core_mass=9, m1_remnant_mass=8, m1_fallback_fraction=0.50, m2=7, P=5)

binary.set_natal_kick_distribution(n_trials=50000, distribution_id='Maxwell',
                                   kick_scaling=lambda x: (1-binary.m1_fallback_fraction)*x)

binary.get_natal_kick_distribution()

binary.get_orbital_distribution(verbose=True)

binary.get_post_kick_grid(use_unbounded_for_norm=True, verbose=True)

binary.save_complete_grid(kick_fname='kick-grid.data', orbit_fname='orbit.data')

import unittest

import numpy as np

import poskiorb


class TestBinary(unittest.TestCase):
    """Test binary module"""

    def test_blaauw_kick_no_mass_loss(self):
        """Blaauw kick"""

        # load a known case
        bS = poskiorb.BinarySystem(
            m1=50, m1_core_mass=50, m1_remnant_mass=50, m1_fallback_fraction=0, m2=50, P=2
        )
        bS.set_single_natal_kick(w=0, theta=0, phi=0)

        bS.get_orbital_distribution()

        self.assertTrue(np.isclose(bS.P_post, bS.P), "unknown changes detected for Porb")
        self.assertTrue(np.isclose(bS.e_post, 0), "unknown changes detected for e")

    def test_blaauw_kick_with_mass_loss(self):
        """Blaauw kick with mass loss"""

        # repeat, but now assuming a 10 Msun lost during collapse
        bS = poskiorb.BinarySystem(
            m1=50, m1_core_mass=50, m1_remnant_mass=40, m1_fallback_fraction=0, m2=50, P=2
        )
        bS.set_single_natal_kick(w=0, theta=0, phi=0)

        bS.get_orbital_distribution()

        self.assertTrue(np.isclose(bS.P_post, 2.52, atol=0.01), f"changes detected in Porb: {bS.P_post}")
        self.assertTrue(np.isclose(bS.e_post, 0.11, atol=0.01), f"changes detected in e: {bS.e_post}")

    def test_kick_distributions(self):
        """Kick distribution"""

        # load a known case
        bS = poskiorb.BinarySystem(
            m1=12.81, m1_core_mass=9, m1_remnant_mass=8, m1_fallback_fraction=0.50, m2=7, P=5
        )
        bS.set_natal_kick_distribution(
            n_trials=10,
            distribution_id="Maxwell",
            seed=1000,
            kick_scaling=lambda x: (1 - bS.m1_fallback_fraction) * x,
        )

        bS.get_natal_kick_distribution()

        # this are the true values we look
        w_should_be = [
            332.57518829,
            117.39083776,
            241.72212351,
            294.38411899,
            153.10462065,
            167.582432,
            195.94083464,
            210.24296413,
            174.55265852,
            407.11161159,
        ]
        theta_should_be = [
            1.25856886,
            2.44961571,
            0.4497272,
            1.60642106,
            0.73033783,
            2.18380955,
            2.73527078,
            1.77788441,
            2.13380819,
            0.81827496,
        ]
        phi_should_be = [
            1.30113674,
            4.66507366,
            2.46397705,
            1.14515149,
            4.67179593,
            0.43719711,
            5.56273771,
            5.98564129,
            5.85054675,
            2.61022966,
        ]

        self.assertTrue(np.allclose(bS.w, w_should_be), "kick strength not equal")
        self.assertTrue(np.allclose(bS.theta, theta_should_be), "theta not equal")
        self.assertTrue(np.allclose(bS.phi, phi_should_be), "phi not equal")

    def test_grid(self):
        """Grid of binaries after kick"""

        # load a known case
        bS = poskiorb.BinarySystem(
            m1=12.81, m1_core_mass=9, m1_remnant_mass=8, m1_fallback_fraction=0.50, m2=7, P=5
        )
        bS.set_natal_kick_distribution(
            n_trials=10,
            distribution_id="Maxwell",
            seed=1000,
            kick_scaling=lambda x: (1 - bS.m1_fallback_fraction) * x,
        )

        bS.get_natal_kick_distribution()
        bS.get_orbital_distribution(verbose=True)
        bS.get_post_kick_grid(use_unbounded_for_norm=True, verbose=True)

        grid_should_be = [
            0.1,
            0.1,
            0.1,
        ]

        self.assertTrue(np.allclose(bS.prob_grid, grid_should_be), f"grid not equal: {bS.prob_grid}")

import unittest

import numpy as np

import poskiorb


class TestKicks(unittest.TestCase):
    """Test kicks module"""

    def test_kicks_distributions(self):
        """Test on \vec{w} = (w, theta, phi)"""

        poskiorb.set_seed(seed=100)

        # number of kicks
        N = 50000

        w = poskiorb.kick_velocity_distribution(distribution="Maxwell", N=N, sigma=265e0)
        theta = poskiorb.theta_distribution(N=N)
        phi = poskiorb.phi_distribution(N=N)

        self.assertTrue(np.mean(w) - 439.23, f"mean w distribution not matching: {np.mean(w)}")
        self.assertTrue(
            np.mean(theta) - 1.58, f"mean theta distribution not matching: {np.mean(theta)}"
        )
        self.assertTrue(np.mean(phi) - 3.03, f"mean phi distribution not matching: {np.mean(phi)}")

import unittest

import numpy as np

import poskiorb


class TestUtils(unittest.TestCase):
    """Test utils module"""

    Msun = poskiorb.constants.Msun
    Mearth = poskiorb.constants.m_earth
    Rsun = poskiorb.constants.Rsun
    au = poskiorb.constants.au
    e = 0.0167
    mean_v_earth = 29.78
    max_v_earth = 30.29
    min_v_earth = 29.29

    def test_v_orb(self):
        """Orbital velocity"""

        mean_v_orb = poskiorb.v_orb(
            r=self.au / self.Rsun, m1=self.Mearth / self.Msun, m2=1, separation=self.au / self.Rsun
        )
        max_v_orb = poskiorb.v_orb(
            r=(self.au * (1 - self.e)) / self.Rsun,
            m1=self.Mearth / self.Msun,
            m2=1,
            separation=self.au / self.Rsun,
        )
        min_v_orb = poskiorb.v_orb(
            r=(self.au * (1 + self.e)) / self.Rsun,
            m1=self.Mearth / self.Msun,
            m2=1,
            separation=self.au / self.Rsun,
        )

        mean_rel_diff = abs(mean_v_orb - self.mean_v_earth) / self.mean_v_earth
        max_rel_diff = abs(max_v_orb - self.max_v_earth) / self.max_v_earth
        min_rel_diff = abs(min_v_orb - self.min_v_earth) / self.min_v_earth

        self.assertTrue(
            mean_rel_diff < 1e-3, f"earth orbital velocity not matching: rel = {mean_rel_diff}"
        )
        self.assertTrue(
            max_rel_diff < 1e-3, f"earth orbital velocity not matching: rel = {max_rel_diff}"
        )
        self.assertTrue(
            min_rel_diff < 1e-3, f"earth orbital velocity not matching: rel = {min_rel_diff}"
        )

    def test_kepler(self):
        sep = poskiorb.P_to_a(period=365.25, m1=self.Mearth / self.Msun, m2=1)
        porb = poskiorb.a_to_P(separation=self.au / self.Rsun, m1=self.Mearth / self.Msun, m2=1)

        sep_mean_rel_diff = abs(sep - (self.au / self.Rsun)) / (self.au / self.Rsun)
        porb_mean_rel_diff = abs(porb - 365.25) / 365.25

        self.assertTrue(
            sep_mean_rel_diff < 1e-3,
            f"earth-sun separation not matching: rel = {sep_mean_rel_diff}",
        )
        self.assertTrue(
            porb_mean_rel_diff < 1e-3,
            f"earth-sun orbital period not matching: rel = {porb_mean_rel_diff}",
        )

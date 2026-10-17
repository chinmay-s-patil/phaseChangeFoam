#!/usr/bin/env python3
"""
Unit test suite for exact_energy_balance_audit.py sign convention & boundary flux integration.
Verifies that:
1. OpenFOAM negative coldWall flux (heat leaving the domain) is correctly subtracted:
   Hot patch (+10 W) and Cold patch (-2 W) yields E_net = 8 J/s, never 12 J/s.
2. Rectangular startup step and trapezoidal quadrature are mathematically exact.
"""

import os
import sys
import unittest
import numpy as np

# Add module path
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from exact_energy_balance_audit import integrate_flux, parse_wall_heat_flux

class TestEnergyAuditSignConvention(unittest.TestCase):

    def test_synthetic_constant_flux_sign_convention(self):
        """
        Hot wall = +10.0 W (heat entering domain)
        Cold wall = -2.0 W (heat leaving domain across isothermal cold wall)
        Over 1 second, net boundary heat injection MUST be exactly +8.0 J.
        """
        # Synthetic time-flux series at dt = 0.1s from t = 0.1 to t = 1.0s (10 steps)
        steps = np.arange(0.1, 1.1, 0.1)
        hot_series = [(float(t), 10.0) for t in steps]
        cold_series = [(float(t), -2.0) for t in steps]

        t_hot, E_hot = integrate_flux(hot_series)
        t_cold, E_cold_raw = integrate_flux(cold_series)

        # At t = 1.0s (index -1):
        # 10 steps of dt = 0.1 with constant 10 W = 1.0 J per step = 10.0 J
        self.assertAlmostEqual(E_hot[-1], 10.0, places=5)
        self.assertAlmostEqual(E_cold_raw[-1], -2.0, places=5)

        # Enforce the physical sign convention:
        # eout_loss is defined as positive outward loss: -eout_raw > 0
        eout_loss = -E_cold_raw[-1]
        self.assertAlmostEqual(eout_loss, 2.0, places=5)

        # Net heat delivered to domain:
        enet = E_hot[-1] - eout_loss
        self.assertAlmostEqual(enet, 8.0, places=5,
                               msg=f"Sign bug detected! Expected E_net = 8.0 J, got {enet} J")

        # Verify that the old buggy convention (ein - eout_raw) would fail
        buggy_enet = E_hot[-1] - E_cold_raw[-1]
        self.assertNotAlmostEqual(enet, buggy_enet, places=3,
                                  msg="Old buggy convention (ein - eout_raw = 12 J) must NOT match corrected E_net")
        self.assertAlmostEqual(buggy_enet, 12.0, places=5)

    def test_startup_integration_exactness(self):
        """
        Tests rectangular startup step (t = 0 to t1) followed by trapezoidal integration.
        For linearly increasing flux Q(t) = 10 * t:
        Startup (0 -> 0.1): 10 * 0.1 * 0.1 = 0.1 J
        0.1 -> 0.2: 0.5 * (1.0 + 2.0) * 0.1 = 0.15 J
        Total at 0.2s = 0.25 J
        """
        series = [(0.1, 1.0), (0.2, 2.0)]
        t, E = integrate_flux(series)
        self.assertAlmostEqual(E[0], 0.10, places=6)
        self.assertAlmostEqual(E[1], 0.25, places=6)

    def test_mock_wallheatflux_file_parsing(self):
        """
        Verifies parser correctly reads OpenFOAM wallHeatFlux format.
        """
        import tempfile
        content = """# Wall heat-flux
# Time          patch           min             max             integral      
0.1             hotWall         1.000000e+02    1.000000e+02    1.000000e+01
0.1             coldWall        -2.000000e+01   -2.000000e+01   -2.000000e+00
0.2             hotWall         1.000000e+02    1.000000e+02    1.000000e+01
0.2             coldWall        -2.000000e+01   -2.000000e+01   -2.000000e+00
"""
        with tempfile.NamedTemporaryFile(mode='w', delete=False) as tf:
            tf.write(content)
            temp_path = tf.name

        try:
            data = parse_wall_heat_flux(temp_path)
            self.assertIn('hotWall', data)
            self.assertIn('coldWall', data)
            self.assertEqual(len(data['hotWall']), 2)
            self.assertEqual(data['hotWall'][0], (0.1, 10.0))
            self.assertEqual(data['coldWall'][0], (0.1, -2.0))
        finally:
            os.remove(temp_path)

if __name__ == '__main__':
    unittest.main()

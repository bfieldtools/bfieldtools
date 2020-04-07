"""
Disc validation
=========================

This example validates the thermal noise computations against the analytical solution for a thin disc.
"""

import numpy as np


kB = 1.38064852e-23  # Boltzman constant
mu0 = 4 * np.pi * 1e-7
d = 100e-9  # Film thickness in meters
T = 273 + 160  # Temperature in C


r = 100e-6

z = np.linspace(0.2e-3, 1.7e-3)


def plat_sigma(T):
    ref_T = 293  # K
    ref_rho = 1.06e-7  # ohm*meter
    alpha = 0.00392  # 1/K

    rho = alpha * (T - ref_T) * ref_rho + ref_rho

    return 1 / rho


mu0 * np.sqrt((3 * plat_sigma(T) * kB * T * d) / (2048 * np.pi)) * (2 * r) / (z ** 2)

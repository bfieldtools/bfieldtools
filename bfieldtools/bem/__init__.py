"""Boundary elements for piecewise homogeneous magnetostatic media.

Potentials follow ``H = -grad(U)``. Surface densities are constant per face;
they are distinct from the vertex stream functions used elsewhere in bfieldtools.
"""

from .boundary import MagneticBoundary
from .operators import (
    adjoint_double_layer,
    single_layer_gradient,
    single_layer_potential,
)

__all__ = [
    "MagneticBoundary",
    "adjoint_double_layer",
    "single_layer_gradient",
    "single_layer_potential",
]

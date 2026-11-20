import numpy as np

from bfieldtools.line_magnetics import magnetic_field
from bfieldtools.line_magnetics import scalar_potential
from bfieldtools.line_magnetics import vector_potential


def create_3d_grid(xx, yy, zz):
    """Creates a direct product grid from three 1D arrays (xx, yy and zz)
    that is appropriately formated for `scalar_potential` and `magnetic_field`.
    """
    X, Y, Z = np.meshgrid(xx, yy, zz, indexing="ij")

    x = X.ravel()
    y = Y.ravel()
    z = Z.ravel()

    return np.array([x, y, z]).T


def test_field_calculations_in_grid():
    """
    The test checks the relation B = -mu0*grad(U) between the scalar potential and
    the magnetic field, and B = rot(A) between the vector potential and the magnetic
    field produced by a rectangular current loop. It compares magnetic fields
    computed by three different methods and prints out a relative error and a
    warning if the error is too large.

    magnetic_field, scalar_potential and vector_potential use independent analytic
    formulas and it is useful to check their consistency.

    The test also compares the magnetic field of a rectangular loop to values
    computed externally and prints a notification if the discrepancy is found
    to be too large.
    """
    # Defines a rectangular current loop offset by ly from the xz plane
    length = 0.3  # (m)
    width = 6e-2  # (m)
    ly = 0  # (m)
    loop_points = np.array(
        [
            [width / 2, ly, length / 2],
            [width / 2, ly, -length / 2],
            [-width / 2, ly, -length / 2],
            [-width / 2, ly, length / 2],
            [width / 2, ly, length / 2],
        ]
    )

    # Defines a grid inside a 3D box
    nx = ny = nz = 10

    # Box sizes
    lx = ly = lz = 0.4  # (m)

    xx = np.linspace(-lx / 2, lx / 2, nx)
    yy = np.linspace(-ly / 2, ly / 2, ny)
    zz = np.linspace(-lz / 2, lz / 2, nz)

    grid = create_3d_grid(xx, yy, zz)

    d = 1e-6  # Numerical differentiation step (m)

    # Defines offseted grids for numerical differentiation
    grid_xp = create_3d_grid(xx + d, yy, zz)
    grid_yp = create_3d_grid(xx, yy + d, zz)
    grid_zp = create_3d_grid(xx, yy, zz + d)

    grid_xm = create_3d_grid(xx - d, yy, zz)
    grid_ym = create_3d_grid(xx, yy - d, zz)
    grid_zm = create_3d_grid(xx, yy, zz - d)

    # The scalar potential at r+-dx, r+-dy and r+-dz
    U_xp = scalar_potential(loop_points, grid_xp)
    U_xm = scalar_potential(loop_points, grid_xm)
    U_yp = scalar_potential(loop_points, grid_yp)
    U_ym = scalar_potential(loop_points, grid_ym)
    U_zp = scalar_potential(loop_points, grid_zp)
    U_zm = scalar_potential(loop_points, grid_zm)

    # The vector potential at r+-dx, r+-dy and r+-dz
    A_xp = vector_potential(loop_points, grid_xp)
    A_xm = vector_potential(loop_points, grid_xm)
    A_yp = vector_potential(loop_points, grid_yp)
    A_ym = vector_potential(loop_points, grid_ym)
    A_zp = vector_potential(loop_points, grid_zp)
    A_zm = vector_potential(loop_points, grid_zm)

    # The calculation of magnetic field directly and via scalar potential and
    # compare the two methods

    mu0 = 4 * np.pi * 1e-7

    # Compute the magnetic field from scalar potential, B = -mu0*U
    B = -(mu0 / (2 * d)) * np.vstack([(U_xp - U_xm), (U_yp - U_ym), (U_zp - U_zm)]).T

    # Compute the madnetic field directly
    B_2 = magnetic_field(loop_points, grid)

    # Compute the magnetic field from vector potential, B = rot(A)
    B_3 = (
        1
        / (2 * d)
        * np.vstack(
            [
                (A_yp[:, 2] - A_ym[:, 2]) - (A_zp[:, 1] - A_zm[:, 1]),
                -(A_xp[:, 2] - A_xm[:, 2]) + (A_zp[:, 0] - A_zm[:, 0]),
                (A_xp[:, 1] - A_xm[:, 1]) - (A_yp[:, 0] - A_ym[:, 0]),
            ]
        ).T
    )

    error_bounds = [5e-8, 2e-7, 1e-9]

    # Print the maximum relative deviation between the fields calculated by
    # first two methods.
    err0 = np.max(np.abs((B_2 - B) / (np.abs(B) + np.abs(B_2))))

    if err0 < error_bounds[0]:
        print("The scalar potential test is passed.")
    else:
        print(
            "The scalar potential test is failed. The relative error is expected "
            "to be around 3.6e-8."
        )
    print(f"Relative error: {err0}\n")

    err1 = np.max(np.abs((B_2 - B_3) / (np.abs(B_3) + np.abs(B_2))))

    if err1 < error_bounds[1]:
        print("The vector potential test is passed.")
    else:
        print(
            "The vector potential test is failed. The relative error is expected "
            "to be around 1e-7."
        )
    print(f"Relative error: {err1}\n")

    # Calsulates the magnetic field on the diagonal of the grid.
    B_lin = magnetic_field(loop_points, np.vstack([xx, yy, zz]).T)

    # B_ref is magnetic field of the same configuration of currents computed using
    # a different external program.
    B_ref = [
        [5.0405815897640994e-8, 5.692103512792891e-9, 3.786395675311007e-8],
        [1.130657032079328e-7, 1.776261166601295e-8, 7.131528630113071e-8],
        [3.2247298177571e-7, 6.873374013470318e-8, 1.3724305960053128e-7],
        [1.2321797579182563e-6, 3.0821577841363554e-7, 1.7259933917883818e-7],
        [6.628725320025185e-6, 6.310902639454654e-6, 3.273181280102712e-8],
        [6.628725320025185e-6, 6.310902639454654e-6, 3.273181280102712e-8],
        [1.2321797579182563e-6, 3.0821577841363554e-7, 1.7259933917883818e-7],
        [3.2247298177571e-7, 6.873374013470318e-8, 1.3724305960053128e-7],
        [1.130657032079328e-7, 1.776261166601295e-8, 7.131528630113071e-8],
        [5.0405815897640994e-8, 5.692103512792891e-9, 3.786395675311007e-8],
    ]

    err2 = np.max(np.abs((B_lin - B_ref) / (np.abs(B_lin) + np.abs(B_ref))))

    if err2 < error_bounds[2]:
        print("The magnetic field test is passed.")
    else:
        print(
            "The magnetic field test is failed. The relative error is expected "
            "to be below 1e-9."
        )
    print(f"Relative error: {err2}")

    assert err0 < error_bounds[0] and err1 < error_bounds[1] and err2 < error_bounds[2]

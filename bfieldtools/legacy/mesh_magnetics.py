import time
import numpy as np

from ..mesh_calculus import gradient_matrix
from .integrals import omega, gamma0


def magnetic_field_coupling_analytic(mesh, r, Nchunks=None):
    """
    **DEPRECATED**
    Given a mesh, computes the "C matrix" which gives the magnetic field at
    some target points due to currents (stream function) on a surface mesh.

    Parameters
    ----------

    mesh: Trimesh mesh object
            mesh describing the geometry of the field source
    r: ndarray (Np, 3)
        evaluation points

    Returns
    -------
    C: ndarray (Np, 3, Nvertices)
        Coupling matrix corresponding to a mapping from a stream function
        on the mesh to B-field at the evaluation points

    **DEPRECATED**
    """

    coef = 1e-7

    print(
        "Computing magnetic field coupling matrix analytically, %d vertices by %d target points... "
        % (len(mesh.vertices), len(r)),
        end="",
    )
    start = time.time()

    if Nchunks is None:
        if r.shape[0] > 1000:
            Nchunks = r.shape[0] // 100
        else:
            Nchunks = 1

    #    ta = mesh.area_faces
    tn = mesh.face_normals

    # Nfaces, 3, 3
    rfaces = mesh.vertices[mesh.faces]
    # Neval, Nfaces, xyz
    coeffs = np.zeros((r.shape[0:1] + mesh.faces.shape[:-1] + (3,)))
    # Nfaces, Nedges, 3
    edges = np.roll(rfaces, 1, -2) - np.roll(rfaces, 2, -2)
    #    grad = np.cross(tn[:, None, :], edges, axis=-1)/(2*ta[:, None, None])
    Gx, Gy, Gz = gradient_matrix(mesh, rotated=False)
    Rx, Ry, Rz = gradient_matrix(mesh, rotated=True)
    solid_angle = np.zeros((r.shape[0], rfaces.shape[0]))
    # Calculate potentials and related coefficients
    for n in range(Nchunks):
        RRchunk = r[n::Nchunks, None, None, :] - rfaces[None, :, :, :]
        # Neval, Nfaces, Nedges
        #        result = -np.sum(np.sum(gamma0(RRchunk)[..., None]*edges,
        #                               axis=-2)[...,None, :]*edges, axis=-1)
        #        result *= 1/(2*ta[..., :, None])
        solid_angle[n::Nchunks] = omega(RRchunk)
        coeffs[n::Nchunks] = -np.einsum("...i,...ik->...k", gamma0(RRchunk), edges)

    #    # Accumulate the elements
    C = np.zeros((3, r.shape[0], mesh.vertices.shape[0]))
    #    for ind_f, f in enumerate(mesh.faces):
    #        C[:, f, :] += coeffs[:, ind_f, :, None]*tn[ind_f]
    #        C[:, f, :] += -solid_angle[:, ind_f:ind_f+1, None]*grad[ind_f]

    G = (Gx, Gy, Gz)
    for i in range(3):
        cc = coeffs * tn[:, i, None]
        C[i] += cc[:, :, 0] @ Rx + cc[:, :, 1] @ Ry + cc[:, :, 2] @ Rz
        C[i] += -solid_angle @ G[i]

    duration = time.time() - start
    print("took %.2f seconds." % duration)

    C *= coef
    #    return np.moveaxis(C, 0, 1)
    return np.swapaxes(C, 0, 1)

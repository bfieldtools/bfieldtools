"""
Legacy code for sphtools
"""


def compute_sphcoeffs_mesh(mesh, lmax):
    """
    Computes multipole moment (spherical harmonics coefficient) transformation
    from the mesh.

    Parameters
    ----------
    mesh: mesh object - the surface mesh
    lmax: int
        maximum degree l of the fit

    Returns
    -------
    alm: (lmax*(lmax+2)xNvertices array
          transformation from the mesh to alm coefficients (r**l-terms)
    blm: (lmax*(lmax+2)xNvertices array
          transformation from the mesh to blm coefficients (r**(-l)-terms)

    """

    Gx, Gy, Gz = gradient_matrix(mesh, rotated=True)
    tri_normals, tri_areas = tri_normals_and_areas(mesh.vertices, mesh.faces)

    centers = np.mean(mesh.vertices[mesh.faces], axis=1)
    centers_sp = cartesian2spherical(centers)

    alm = np.zeros((lmax * (lmax + 2), mesh.vertices.shape[0]))
    blm = np.zeros((lmax * (lmax + 2), mesh.vertices.shape[0]))

    idx = 0
    for l in range(1, lmax + 1):
        for m in range(-1 * l, l + 1):
            derylm = np.sqrt(l * (l + 1)) * Blm(
                l, m, centers_sp[:, 1], centers_sp[:, 2]
            )
            derylm = sphvec2cart(centers_sp, derylm)
            # FIX: gradient of Ylm has also 1/r multiplier
            derylm /= centers_sp[:, 0:1]

            crossp = np.cross(centers, derylm)
            alm_terms = crossp.T * (centers_sp[:, 0] ** l) * tri_areas
            blm_terms = crossp.T * (centers_sp[:, 0] ** (-l - 1)) * tri_areas
            #            integral_alm = np.sum(G*alm_terms[:,:,None], axis=(0,1))
            #            integral_blm = np.sum(G*blm_terms[:,:,None], axis=(0,1))

            integral_alm = alm_terms[0] @ Gx + alm_terms[1] @ Gy + alm_terms[2] @ Gz
            integral_blm = blm_terms[0] @ Gx + blm_terms[1] @ Gy + blm_terms[2] @ Gz

            # FIX: l-coefficients here were swapped
            blm[idx] = 1 / ((2 * l + 1) * l) * integral_blm
            # FIX: REMOVED MINUS SIGN
            alm[idx] = (
                -1 / ((2 * l + 1) * (l + 1)) * integral_alm
            )  # THIS SIGN SHOULD BE CHECKED TOO
            #            for i in range(mesh.vertices.shape[0]):
            #                G = np.zeros(crossp.shape)
            #                G[:,0] = Gx[:,i]
            #                G[:,1] = Gy[:,i]
            #                G[:,2] = Gz[:,i]
            #                dotp = np.sum(G*crossp,axis=1)
            #                integral_blm = np.sum(dotp*centers_sp[:,0]**l*tri_areas)
            #                blm[idx,i] = -1/((2*l+1)*(l+1))*integral_blm
            #
            #                integral_alm = np.sum(dotp*centers_sp[:,0]**(-l-1)*tri_areas)
            #                alm[idx,i] = 1/((2*l+1)*l)*integral_alm

            idx += 1
        print("l = %d computed" % (l))

    return alm, blm

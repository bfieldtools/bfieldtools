import numpy as np
from numba import jit, prange
from utils import tri_normals_and_areas

def get_neighbour_vertices(vertices, edges):
    '''
    Computes and returns the neighbor vertex indices for all vertices
    '''
#    return [edge[edge[0] == index] for edge in edges if index in edge]
    vi = []
    for vert_idx in range(len(vertices)):
        vi.append([])
        for edge in edges:
            if edge[0] == vert_idx:
                vi[vert_idx].append(edge[1])
            elif edge[1] == vert_idx:
                vi[vert_idx].append(edge[0])

        vi[vert_idx] = np.unique(vi[vert_idx])

    return vi


def get_vert_links(verts, tris):
    '''
    Computes and returns the triangles that each vertex corresponds to.
    '''

    vert_links = []
    for m in range(len(verts)):
        vert_links.append(np.where(tris == m)[0])

    return vert_links


def create_basis(verts, tris, centre=np.array([0, 0, 0]), tri_areas=None, tri_normals=None):
    '''
    Calculate basis function for each triangle
    '''
    n_verts = len(verts)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_normals, type(None)) or isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)

#    vert_links = get_neighbour_vertices(verts, edges)

    vmi = []
    for m in range(n_verts):

        vmi.append([])

        vert_links, vert_idx = np.where(tris == m)

        for tri_idx in vert_links:

            #indices within triangle
            m_vert_idx = np.where(tris[tri_idx] == m)[0][0]
            a_vert_idx = (m_vert_idx + 1) % 3
            b_vert_idx = (m_vert_idx + 2) % 3


            #A - coordinates of the second node of the i_th face of the m_th node
            Ami = verts[a_vert_idx]
            #B - coordinates of the third node of the i_th face of the m_th node
            Bmi = verts[b_vert_idx]
            #C - coordinates of the first node of the i_th face of the m_th node
            Cmi = verts[m]


            vectorOA = Ami - centre

            dotprod = np.dot(tri_normals[tri_idx], vectorOA)


            if dotprod > 0:
                vmi[m].append((Ami - Bmi)/tri_areas[tri_idx])
            else:
                vmi[m].append((Bmi - Ami)/tri_areas[tri_idx])

    return dict(v=vmi, A=Ami, B=Bmi, C=Cmi)






@jit
def compute_L(verts1, tris1, basis1=None, verts2=None, tris2=None, basis2=None,
              vert_links1=None, vert_links2=None,
              tri_areas1=None, tri_areas2=None,
              mu0=4*np.pi*1e-7):
    '''
    Compute and return mutual inductance matrix.
    '''

    #If second mesh not given, compute L within surface
    if verts2 is None:
        verts2 = verts1
        tris2 = tris1
        basis2 = basis1

        same_surface = 1

    else:
        same_surface = 0

    #Compute vertex links if not provided
    if vert_links1 is None:
        vert_links1 = get_vert_links(verts1, tris1)

        if same_surface:
            vert_links2 = vert_links1
        else:
            vert_links2 = get_vert_links(verts2, tris2)

    #Compute triangle areas if not provided
    if tri_areas1 is None:
        tri_normals1, tri_areas1 = tri_normals_and_areas(verts1, tris1)

        if same_surface:
            tri_normals2, tri_areas2 = tri_normals1, tri_areas1
        else:
            tri_normals2, tri_areas2 = tri_normals_and_areas(verts2, tris2)

    #Compute basis if not provided
    if basis1 is None:
        basis1 = create_basis(verts1, tris1, tri_areas=tri_areas1)

        if same_surface:
            basis2 = basis1
        else:
            basis2 = create_basis(verts2, tris2, tri_areas=tri_areas2)


    # Gauss Legendre integration
    cl = ck = triGaussPoints(2) #FIX QUADRATURE POINTS

    n_integration_points = len(ck)

    coef = mu0 / (4*np.pi)

    n_verts1 = len(verts1)
    n_verts2 = len(verts2)

    L = np.full((n_verts1, n_verts2), fill_value=np.nan)

    # Iterate though matrix, utilizing symmetry if possible
    for m in prange(n_verts1):
        if same_surface:
            start_second_loop = m
        else:
            start_second_loop = 0

        for n in range(start_second_loop, n_verts2):

            superBigSum = 0

            for i in range(vert_links1[m]):

                currentTriangle_i = vert_links1[m][i]

                vmi = basis1['v'][m][i]

                r_o = basis1[m][i]['r_o'] #FIX

                for j in range(vert_links2[n]):

                    currentTriangle_j = vert_links2[n][j]

                    vnj = basis2['v'][n][j]

                    integral = 0

                    #If the current triangles are not the same,
                    #then we will not share any basis function. We can thus solve that easily
                    if currentTriangle_i != currentTriangle_j or same_surface == 0:
                        r_p = bases2[n][j]['r_o']

                        for k in range(n_integration_points):
                            for l in range(n_integration_points):
                                integral += cl[l]/np.linalg.norm(r_o[k] - r_p[l])

                            integral *= ck[k]

                        bigSum = integral * (2 * tri_areas2[currentTriangle_j] *
                                             2 * tri_areas1[currentTriangle_i])

                    # If the 2 triangle are similar, we use an approximate
                    # calculation, according to page 72 equ 5.40 of Poole's
                    # thesis

                    else:
                        Ami = basis1['A'][m][i]
                        Bmi = basis1['B'][m][i]
                        Cmi = basis1['C'][m][i]

                        a = np.dot(Cmi - Ami, Cmi - Ami)
                        b = np.dot(Cmi - Ami, Cmi - Bmi)
                        c = np.dot(Cmi - Bmi, Cmi - Bmi)

                        sa = np.sqrt(a)
                        sc = np.sqrt(c)
                        ss = np.sqrt(a - 2*b + c)
                        sac = np.sqrt(a * c)

                        bigSum = (1 * (4 * tri_areas1[currentTriangle_i]**2))                   \
                        * (1 / (6 * sa) * np.log(((a - b + sa * ss) * (b + sac)) /              \
                                               ((-b + sac) * (-a + b +sa*ss)))                  \
                        + 1 / (6 * sc) * np.log(((b + sac) * (-b + c + sc * ss)) /              \
                                                ((b - c + sc * ss) * (-b + sac)))               \
                        + 1 / (6 * ss) * np.log(((a - b + sa * ss) * (-b + c + sc * ss)) /      \
                                                ((b - c + sc * ss) * (-a + b + sa * ss) )))

                    vmi_vnj = np.dot(vmi, vnj)

                    superBigSum += vmi_vnj*bigSum

            superBigSum *= coef

        L[m] = superBigSum


#@jit
def compute_R(verts, tris, basis=None, vert_links=None, tri_areas=None, rho=1.68*1e-8, t=1e-4):
    '''
    Computes resitivity matrix for surface mesh made of material with resistivity rho and thickness t

    Equation come from the thesis of Michael Poole "Improved Equipment and
    %Techniques for Dynamic Shimming in High Field MRI
    page 66
    and G. N. Peeren, Stream function approach for determining optimal surface
    currents


    '''

    n_verts = len(verts)
    #Initialize sparse matrix
    R = np.full((n_verts, n_verts), fill_value=np.nan)


    surface_resistance = rho/t

    #Compute vertex links if needed
    if vert_links is None:
        vert_links = get_vert_links(verts, tris)

    #Compute basis functions if needded
    if basis is None:
        basis = create_basis(verts, tris)

    # Triangle normals and areas, compute if not provided
    if isinstance(tri_areas, type(None)):
        tri_normals, tri_areas = tri_normals_and_areas(verts, tris)


    # Iterate though upper triangular matrix including diagonal
    for m in prange(n_verts):
        for n in range(m, n_verts):
            bigSum = 0
            for i in range(len(vert_links[m])):

                currentTriangle_m = vert_links[m][i]


                for j in range(len(vert_links[n])):

                    currentTriangle_n = vert_links[n][j]


                    #if share a common node
                    if currentTriangle_m == currentTriangle_n:

                        vmi = basis['v'][m][i]

                        vnj = basis['v'][n][j]

                        bigSum += np.dot(vmi, vnj)*tri_areas[m]

            R[m, n] = surface_resistance * bigSum


    #Fill in lower triangular matrix, R is symmetric
    i_lower = np.tril_indices(R.shape[0], -1)
    R[i_lower] = R.T[i_lower]

    return R


def compute_C(mesh, basis, r):
    return C




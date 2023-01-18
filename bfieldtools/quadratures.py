"""

Quadrature rules for numerical integration

Refs:
D.A. Dunavant https://doi.org/10.1002/nme.1620210612
https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature
R. Womersley https://web.maths.unsw.edu.au/~rsw/Sphere/EffSphDes/sf.html

"""

import numpy as np
from dataclasses import dataclass


@dataclass
class QuadratureRule:
    points: np.ndarray
    weights: np.ndarray

    def __init__(self, points, weights):
        self.points = np.asarray(points)
        self.weights = np.asarray(weights)


def symmetric_triangle_points(w):
    return np.ones((3, 3)) * w + np.eye(3) * (1 - 3 * w)


dunavant_rules = (
    QuadratureRule([[1 / 3, 1 / 3, 1 / 3]], [1]),
    QuadratureRule(
        [[2 / 3, 1 / 6, 1 / 6], [1 / 6, 2 / 3, 1 / 6], [1 / 6, 1 / 6, 2 / 3]],
        [1 / 3] * 3,
    ),
    QuadratureRule(
        [
            [1 / 3, 1 / 3, 1 / 3],
            [3 / 5, 1 / 5, 1 / 5],
            [1 / 5, 3 / 5, 1 / 5],
            [1 / 5, 1 / 5, 3 / 5],
        ],
        [-0.5625] + [0.520833333333333] * 3,
    ),
    QuadratureRule(
        [
            *symmetric_triangle_points(0.445948490915965),
            *symmetric_triangle_points(0.091576213509771),
        ],
        [0.223381589678011] * 3 + [0.109951743655322] * 3,
    ),
)


def gauss_legendre_points(root):
    p = 0.5 * (1 - root)  # Map from [-1,1] to point weights
    return [p, 1 - p]


gauss_legendre_rules = (
    QuadratureRule([gauss_legendre_points(0)], [1]),
    QuadratureRule(
        [gauss_legendre_points(1 / np.sqrt(3)), gauss_legendre_points(-1 / np.sqrt(3))],
        [0.5, 0.5],
    ),
    QuadratureRule(
        [
            gauss_legendre_points(0),
            gauss_legendre_points(3 / np.sqrt(5)),
            gauss_legendre_points(-3 / np.sqrt(5)),
        ],
        [8 / 18, 5 / 18, 5 / 18],
    ),
)


line_rules = {
    "midpoint": QuadratureRule([[0.5, 0.5]], [1]),
    "trapezoid": QuadratureRule([[1, 0], [0, 1]], [0.5, 0.5]),
    "simpson": QuadratureRule([[1, 0], [0.5, 0.5], [0, 1]], [1 / 6, 4 / 6, 1 / 6]),
    "gauss_legendre": gauss_legendre_rules,
}


def get_triangle_rule(method):
    name, degree = method
    rules = {"dunavant": dunavant_rules}
    return rules[name][degree]


def get_line_rule(method):
    if isinstance(method, str):
        return line_rules[method]
    name, degree = method
    return line_rules[name][degree]


def get_spherical_quad_points():
    """
    https://web.maths.unsw.edu.au/~rsw/Sphere/Points/SF/SF29-Nov-2012/sf011.00072
    """

    points = np.array(
        [
            [0.0, 0.0, 1.0],
            [0.378287016561, 0.0, 0.92568835636],
            [0.28407948391, 0.419547319235, 0.862136238505],
            [-0.099456444681, 0.466447295733, 0.878939893232],
            [-0.37739769508, 0.140853065159, 0.915276676084],
            [-0.273773355171, -0.341438850481, 0.899148297769],
            [0.133827977641, -0.434471906342, 0.890687506929],
            [0.49969539963, -0.429434648675, 0.752256864446],
            [0.731487037715, -0.109100761224, 0.673070380833],
            [0.656856370862, 0.306973246447, 0.688699596359],
            [0.389726059908, 0.721724436417, 0.572037967365],
            [0.018273797033, 0.81276113642, 0.582310401304],
            [-0.396760437216, 0.628063259952, 0.669415937186],
            [-0.676876198655, 0.226871042497, 0.700262909035],
            [-0.62077531104, -0.241528059766, 0.745856695048],
            [-0.398201099405, -0.664199428633, 0.632672903984],
            [-0.003640875597, -0.752688389465, 0.658366867627],
            [0.359014250301, -0.788403950909, 0.499527755256],
            [0.709162025463, -0.584840414266, 0.393765046039],
            [0.915878280266, -0.260023450188, 0.305867260571],
            [0.92471224043, 0.163019134789, 0.34399423555],
            [0.7334498165, 0.568768547289, 0.372228029963],
            [0.443104014007, 0.882357646087, 0.15844184158],
            [0.055285932806, 0.980960651193, 0.186171067688],
            [-0.351395954503, 0.878824716701, 0.322781660686],
            [-0.714034318096, 0.559274656785, 0.421149440056],
            [-0.91000883805, 0.100250718404, 0.402285605175],
            [-0.834360560251, -0.31822795447, 0.450081575374],
            [-0.646520204723, -0.685495983611, 0.334823657078],
            [-0.227405618009, -0.929945482636, 0.288943046676],
            [0.140561845255, -0.972464383399, 0.185890803106],
            [0.553653283943, -0.830542709124, 0.060554516769],
            [0.834840868173, -0.545742953019, -0.072148139665],
            [0.985344838353, -0.140559763018, -0.096635927845],
            [0.957169570037, 0.288327528355, -0.026337247151],
            [0.752131277497, 0.658031296075, -0.035963242314],
            [0.363490194323, 0.895363912555, -0.257290385994],
            [-0.049181686978, 0.975957950784, -0.212337561368],
            [-0.443415079946, 0.895144639151, -0.045816393086],
            [-0.77239733488, 0.63268954546, 0.055734155921],
            [-0.965048679778, 0.257714332397, 0.047585381526],
            [-0.979932414691, -0.196278110146, 0.0347471742],
            [-0.835279746206, -0.549741306073, -0.009604268536],
            [-0.476541494622, -0.878259023643, -0.039614281438],
            [-0.093714804004, -0.974027403192, -0.206126546912],
            [0.32742996315, -0.911435565395, -0.24914820763],
            [0.68344656136, -0.601865155405, -0.413109104807],
            [0.87497731068, -0.149788219502, -0.460410898104],
            [0.857144230238, 0.289510274013, -0.426013579372],
            [0.640750246807, 0.626804787704, -0.443345101842],
            [0.215763026167, 0.752918708661, -0.621739281924],
            [-0.219047994337, 0.819550210165, -0.529486004721],
            [-0.601092091964, 0.698332831106, -0.388612344089],
            [-0.847915646298, 0.391639391132, -0.357292099096],
            [-0.937539741887, -0.01713695109, -0.347455835021],
            [-0.810716425464, -0.381225542397, -0.444303908722],
            [-0.59208984103, -0.714685113192, -0.372363812865],
            [-0.145967579893, -0.826087402211, -0.544309718386],
            [0.314460467875, -0.731499963639, -0.604997865566],
            [0.597428805095, -0.348651928186, -0.722163870472],
            [0.644583268257, 0.11020334144, -0.75654982243],
            [0.437627186678, 0.430285095091, -0.789517056448],
            [-0.056264905419, 0.529292394121, -0.846571805545],
            [-0.463756814509, 0.481889253054, -0.743446275656],
            [-0.696504925305, 0.155396382186, -0.700523271155],
            [-0.602768636458, -0.218482694245, -0.767421190232],
            [-0.381534762176, -0.571014446067, -0.726893202357],
            [0.096892046972, -0.531430102956, -0.841542617403],
            [0.311557196419, -0.171117473011, -0.93469295696],
            [0.142557019271, 0.174502828245, -0.974282432969],
            [-0.288549869318, 0.124950315214, -0.949276772941],
            [-0.174781523944, -0.243366400441, -0.95405671426],
        ]
    )

    # Don't know if this is needed
    # Some cases there might be singularity on (0,0,1)
    # so this avoids evaluating directly at that point
    points = np.roll(points, 1, axis=1)

    return QuadratureRule(points, np.ones(len(points)) / len(points))


def get_quad_points(verts, tris, method, return_ref_coords=False):
    """Get quad points and weights from quadrature rules for triangles
       in a mesh (verts, tris)

    Parameters
    ----------
    verts: array-like [Nverts x 3]
    tris: array-like [Ntris x 3]
    method: tuple (method_name: str, degree: int)

    Returns
    -------
    w: array-like  (Nquad, )
        quadrature weights
    qp: array-like (Ntris, Nquad, xyz)
        quadrature points in each triangle

    """

    rule = get_triangle_rule(method)
    qp = np.einsum("ij,kjl->kil", rule.points, verts[tris])

    if return_ref_coords:
        return rule.weights, qp, rule.points

    return rule.weights, qp


def get_line_quad_points(line_vertices, method="midpoint", index=None):
    """Get quad points and weights from quadrature rules for edges in
       a closed loop


    Parameters
    ----------
    line_vertices: array-like [Nverts x 3]
        Assumes vertices are connected according to index, last index connects to first
    method: str | tuple(str, int)
        One of 'midpoint', 'trapezoid',  'simpson'
        or ('gauss_legendre', 0 | 1 | 2)

    Returns
    -------
    w: array-like  (Nquad, )
        quadrature weights
    qp: array-like (Nverts, Nquad)
        quadrature points for each edge connecting the vertices

    """
    rule = get_line_rule(method)

    qp = np.zeros((len(line_vertices), len(rule.weights), 3))

    p0 = line_vertices
    p1 = np.roll(line_vertices, -1)

    B = np.array([p0, p1])
    qp = np.einsum("ij,jkl->kil", rule.points, B)

    return rule.weights, qp

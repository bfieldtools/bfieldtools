from . import mesh_conductor
from . import mesh_impedance
from . import mesh_magnetics
from . import mesh_calculus
from . import line_magnetics
from . import line_conductor
from . import coil_optimize
from . import sphtools
from . import suhtools
from . import thermal_noise
from . import utils
from . import integrals
from . import legacy
from . import flatten_mesh

try:
    from .version import __version__
except ModuleNotFoundError:
    print("version.py not present, did you install the package?")
    print("Attempting to get version through importlib")
    from importlib.metadata import version

    __version__ = version("bfieldtools")


__all__ = [
    "line_magnetics",
    "line_conductor",
    "coil_optimize",
    "mesh_calculus",
    "mesh_magnetics",
    "mesh_conductor",
    "mesh_impedance",
    "sphtools",
    "suhtools",
    "thermal_noise",
    "utils",
    "integrals",
    "legacy",
    "flatten_mesh",
]

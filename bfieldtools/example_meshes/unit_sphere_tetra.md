# Unit sphere tetrahedral mesh

`unit_sphere_tetra.npz` stores dimensionless `vertices` and zero-based
`tetrahedra` arrays, loaded with `volume.load_example_tetra_mesh()`.
Scale the vertices explicitly to assign a radius in metres.

Generated with `trimesh.creation.icosphere(subdivisions=2)` and
`volume.tetrahedralize(surface, max_volume=0.01)` using the optional meshing
dependencies. The archive contains 226 vertices and 808 positive tetrahedra;
the tree–cotree current space has 649 independent coefficients. The minimum
mean-ratio cell quality is approximately 0.261. The boundary is a faceted
approximation to a unit sphere. This is a demonstration mesh, not a spatial
convergence reference. Storing the arrays makes example runs independent of
TetGen versions and meshing nondeterminism.

Regenerate from the repository root with:

```python
import numpy as np
import trimesh
from bfieldtools.volume import tetrahedralize

mesh = tetrahedralize(trimesh.creation.icosphere(subdivisions=2), max_volume=0.01)
np.savez_compressed(
    "bfieldtools/example_meshes/unit_sphere_tetra.npz",
    vertices=mesh.vertices,
    tetrahedra=mesh.tetrahedra,
)
```

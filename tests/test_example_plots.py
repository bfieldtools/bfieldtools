"""Geometry and data ordering used by the PyVista examples."""

import numpy as np
import pyvista as pv
from numpy.testing import assert_allclose, assert_array_equal

from bfieldtools._example_plots import (
    mesh_actor,
    polyline,
    structured_field,
    triangle_mesh,
    vector_glyphs,
)


def test_surface_interpolation():
    surface = triangle_mesh([0, 1, 0], [0, 0, 1], [0, 0, 0], [[0, 1, 2]])
    surface["scalars"] = [0.0, 1.0, 2.0]
    sampled = pv.PolyData([[0.25, 0.25, 0.0]]).sample(surface)
    assert_allclose(sampled["scalars"], [0.75])


def test_structured_field_point_order():
    x, y, z = np.meshgrid([0.0, 1.0], [0.0, 2.0], [0.0, 3.0], indexing="ij")
    grid = structured_field(x, y, z, x + y, y + z, z + x)
    px, py, pz = grid.points.T
    assert_allclose(grid["vectors"], np.column_stack((px + py, py + pz, pz + px)))
    assert_allclose(grid["magnitude"], np.linalg.norm(grid["vectors"], axis=1))


def test_vector_glyph_scale_and_colours():
    glyphs = vector_glyphs(
        [0, 0], [0, 1], [0, 0], [1, 2], [0, 0], [0, 0], scale_factor=0.5, scalars=[3, 7]
    )
    assert_allclose(glyphs.bounds[1], 1.0)
    assert_array_equal(np.unique(glyphs["scalars"]), [3, 7])
    zeros = vector_glyphs([0], [0], [0], [0], [0], [0])
    assert np.isfinite(zeros.points).all()


def test_polyline_connections():
    line = polyline([0, 1, 2], [0, 1, 0], [0, 0, 0], [1, 2, 3])
    assert line.n_lines == 1
    assert_array_equal(line.lines, [3, 0, 1, 2])
    assert_array_equal(line["scalars"], [1, 2, 3])


def test_mesh_actor_accepts_plotter_and_actor():
    plotter = pv.Plotter(off_screen=True)
    try:
        actor = plotter.add_mesh(pv.Sphere())
        assert mesh_actor(plotter) is actor
        assert mesh_actor(actor) is actor
    finally:
        plotter.close()


def test_standalone_html_export(tmp_path):
    import pytest
    from html.parser import HTMLParser

    pytest.importorskip("trame")

    class Scripts(HTMLParser):
        def __init__(self):
            super().__init__()
            self.external = []

        def handle_starttag(self, tag, attrs):
            if tag == "script":
                self.external.extend(value for key, value in attrs if key == "src")

    plotter = pv.Plotter(off_screen=True)
    try:
        plotter.add_mesh(pv.Sphere())
        filename = tmp_path / "scene.html"
        plotter.export_html(filename)
        html = filename.read_text()
        assert "</html>" in html.lower()
        scripts = Scripts()
        scripts.feed(html)
        assert not scripts.external
    finally:
        plotter.close()

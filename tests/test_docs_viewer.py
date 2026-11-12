"""Local-file compatibility of the generated documentation viewers."""

import base64
from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest


def test_embed_gallery_scenes(tmp_path):
    pytest.importorskip("trame_vtk")
    extension = runpy.run_path(
        str(Path(__file__).parents[1] / "docs/_ext/standalone_viewer.py")
    )
    page = tmp_path / "auto_examples/section/example.html"
    page.parent.mkdir(parents=True)
    scene = tmp_path / "_images/section/scene.vtksz"
    scene.parent.mkdir(parents=True)
    scene.write_bytes(b"example scene data")
    page.write_text(
        '<iframe src="../../_static/static_viewer.html?fileURL='
        '../_images/section/scene.vtksz"></iframe>\n'
        '<iframe src="https://example.com/viewer"></iframe>'
    )
    app = SimpleNamespace(outdir=tmp_path, builder=SimpleNamespace(format="html"))
    extension["embed_gallery_scenes"](app, None)

    html = page.read_text()
    assert 'src="../../_images/section/scene.html"' in html
    assert 'src="https://example.com/viewer"' in html
    viewer = scene.with_suffix(".html").read_text()
    assert base64.b64encode(scene.read_bytes()).decode() in viewer
    assert "OfflineLocalView.load(container, { base64Str })" in viewer
    assert "fileURL=" not in html

    extension["embed_gallery_scenes"](app, None)
    assert page.read_text() == html

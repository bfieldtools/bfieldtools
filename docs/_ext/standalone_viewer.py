"""Embed PyVista scene data so gallery viewers also work with file:// URLs."""

import os
from pathlib import Path
import re
from urllib.parse import parse_qs, unquote, urlsplit


IFRAME_SOURCE = re.compile(
    r"""(?P<prefix><iframe\b[^>]*\bsrc=)(?P<quote>['"])(?P<url>[^'"]+)(?P=quote)"""
)


def embed_gallery_scenes(app, exception):
    if exception is not None or app.builder.format != "html":
        return

    from trame_vtk.tools.vtksz2html import embed_data_to_viewer_file

    output = Path(app.outdir)
    exported = set()
    for page in output.rglob("*.html"):
        if page.relative_to(output).parts[0] in ("_images", "_static"):
            continue

        def replace_viewer(match):
            url = urlsplit(match["url"])
            assets = parse_qs(url.query).get("fileURL")
            if (
                url.scheme
                or url.netloc
                or not url.path.endswith("/static_viewer.html")
                or not assets
            ):
                return match[0]
            viewer = page.parent / unquote(url.path)
            scene = (viewer.parent / assets[0]).resolve()
            embedded = scene.with_suffix(".html")
            if scene not in exported:
                embed_data_to_viewer_file(scene.read_bytes(), embedded)
                exported.add(scene)
            relative = Path(os.path.relpath(embedded, page.parent)).as_posix()
            return f'{match["prefix"]}{match["quote"]}{relative}{match["quote"]}'

        original = page.read_text(encoding="utf-8")
        updated = IFRAME_SOURCE.sub(replace_viewer, original)
        if updated != original:
            page.write_text(updated, encoding="utf-8")


def setup(app):
    app.connect("build-finished", embed_gallery_scenes)
    return {"parallel_read_safe": True, "parallel_write_safe": True}

# Using uv with bfieldtools

Install [uv](https://docs.astral.sh/uv/getting-started/installation/), then run
these commands from the repository root.

## Development

```bash
uv sync --locked --extra dev
uv run --locked --extra dev pytest tests/
uv run --locked --extra dev black --check .
```

`uv sync` creates `.venv`, installs Python 3.11 as specified in
`.python-version` if needed, and installs bfieldtools in editable mode.
`--locked` uses the committed `uv.lock` and reports an error if it is out of
date. No manual environment creation or activation is needed.

For the package without development tools, use `uv sync --locked`.

## Documentation

```bash
uv sync --locked --extra docs
uv run --locked --extra docs make -C docs html
```

To keep both development and documentation tools installed, pass
`--extra dev --extra docs` to `uv sync` and `uv run`.

## Updating dependencies

Edit `dependencies` or `project.optional-dependencies` in `pyproject.toml`, then:

```bash
uv lock
uv sync --locked --extra dev
```

Commit both `pyproject.toml` and `uv.lock`. To upgrade all dependencies within
the declared constraints, use `uv lock --upgrade`.

## SciPy installation

The project configures uv to install SciPy from a prebuilt wheel, avoiding
source builds that require BLAS/LAPACK libraries. Use the default Python 3.11
environment; overriding it with a newer Python can leave no compatible wheel
for the project's `scipy<1.14.0` constraint.

`uv pip install -e ".[dev]"` also respects this wheel setting, but does not use
`uv.lock`. Prefer `uv sync --locked --extra dev` for reproducible installation.

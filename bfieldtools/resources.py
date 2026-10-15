from importlib import resources


def files(module: str, file: str | None = None):
    """Small wrapper for pkg_resources compatibility"""
    if file is None:
        return resources.files(module)
    return resources.files(module) / file

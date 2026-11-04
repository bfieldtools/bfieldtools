"""
Pytest configuration for bfieldtools tests.

This configuration sets up PyVista for headless/off-screen rendering
to prevent interactive windows from opening during test execution.

To enable interactive plotting during tests, set the environment variable:
    BFIELDTOOLS_TEST_PLOT=1 pytest tests/test_viz.py
"""
import pytest
import os


@pytest.fixture(scope="session", autouse=True)
def configure_pyvista():
    """Configure PyVista for headless rendering during tests."""
    import pyvista as pv

    # Check if interactive plotting is requested
    enable_plotting = os.environ.get("BFIELDTOOLS_TEST_PLOT", "0") == "1"

    if not enable_plotting:
        # Set PyVista to use off-screen rendering (no interactive windows)
        pv.OFF_SCREEN = True

        # Start Xvfb if needed (for CI/CD environments)
        try:
            pv.start_xvfb()
        except Exception:
            # If Xvfb is not available or not needed, continue silently
            pass

    yield

    # Cleanup after all tests
    if not enable_plotting:
        pv.close_all()

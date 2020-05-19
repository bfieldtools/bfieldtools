"""
Module with functions for exporting polylines to the PCB design software KiCAD
"""

__all__ = ["python_to_kicad"]


def python_to_kicad(
    loops, filename, plane_axes, origin, layer, net, scaling=1, trace_width=0.2
):
    """
    Appends polylines to KiCAD PCB files as traces.

    Paramaters
    ----------
        loops: N_loops long list of (Nv_loop, 3) arrays
            Each corresponding to the vertex locations for widning loops
        filename: str
            filename/path where the ifno is written (append mode is used)
        plane_axes: tuple with length 2
            specifies the X- and Y-dimensions used for the PCB
        origin: (2, ) array-like
            Origin shift applid in original units
        layer: str
            Which layer to write to
        scaling: float
            Scaling factor applied to the loops
        trace_width: float
            The trace width in mm
        net: int
            Specifies which net number the loops are assigned to

    Returns
    -------
        None

    """

    with open(filename, "a") as file:
        for loop in loops:
            for seg_idx in range(1, len(loop)):
                x_start = loop[seg_idx - 1, plane_axes[0]] + origin[0]
                y_start = loop[seg_idx - 1, plane_axes[1]] + origin[1]

                x_end = loop[seg_idx, plane_axes[0]] + origin[0]
                y_end = loop[seg_idx, plane_axes[1]] + origin[1]

                file.write(
                    "    (segment (start %.2f %.4f) (end %.2f %.2f) (width %.2f) (layer %s) (net %d))\n"
                    % (
                        x_start * scaling,
                        y_start * scaling,
                        x_end * scaling,
                        y_end * scaling,
                        trace_width,
                        layer,
                        net,
                    )
                )

    return

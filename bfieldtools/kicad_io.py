
def python_to_kicad(loops, filename, plane_axes, origin, layer, net, scaling=1, trace_width=0.2):
    '''

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

    '''

    with open(filename, 'a') as file:
        for loop in loops:
            for seg_idx in range(1, len(loop)):
                    x_start = loop[seg_idx - 1, plane_axes[0]] + origin[0]
                    y_start = loop[seg_idx - 1, plane_axes[1]] + origin[1]

                    x_end = loop[seg_idx, plane_axes[0]] + origin[0]
                    y_end = loop[seg_idx, plane_axes[1]] + origin[1]

                    file.write('    (segment (start %.2f %.4f) (end %.2f %.2f) (width %.2f) (layer %s) (net %d))\n'%
                               (x_start*scaling, y_start*scaling, x_end*scaling, y_end*scaling, trace_width, layer, net))

    return


def kicad_to_python(file_directory, layers, stand_off, scale_factor = 10e2):

    '''
    Author: Elias Lius

    Brings modified coils back to python from kicad-files.
    Inputs:
        file_directory: string
            Directory to the kicad-files
        layers: (1,N-layers) string array
            On which layers this coil is in kicad-file, e.g., ['F.Cu', 'In1.Cu']
        stand_off: float
            Amount of plane's stand off from the origin
        scale_factor: float
            As KiCad handles everything in mm, the wanted scale must be adjusted. Default value is 10e2.

    Outputs:
        vert_start: starts of the line segments in new plane imported from kicad
        vert_end: end of the line segments in new plane imported from kicad

        (Starts and ends help calculation of bfield, as linesegments drawn in kicad go to the end of the segment lists.)
    '''

    import numpy as np
    seg_start = np.zeros((1,3), dtype=np.float64)
    seg_end = np.zeros((1,3), dtype=np.float64)
    filedir = file_directory
    with open(filedir, 'r') as infile:
        k = 0
        for line in infile:
            line = line.replace( ')', '')
            #line = line.replace('(', '')
            txt = line.split()

            if len(txt) > 0 and txt[0] == '(segment' and (txt[-3] == layers[0]  or txt[-3] == layers[1]):
                seg_start[:,1] = stand_off
                seg_end[:,1] = stand_off

                seg_start[k,0] = float(txt[2])/scale_factor
                seg_end[k,0] = float(txt[5])/scale_factor
                seg_start[k,2] = float(txt[3])/scale_factor
                seg_end[k,2] = float(txt[6])/scale_factor


                seg_start = np.append(seg_start, [[0,0,0]], axis = 0)
                seg_end = np.append(seg_end, [[0,0,0]], axis = 0)

                k += 1

        seg_start = seg_start[:-1]
        seg_end = seg_end[:-1]


        return seg_start, seg_end



def shifted_kicad_to_python(file_directory, layers, orig_plane, scale_factor = 10e2):

    '''
    Author: Elias Lius

    Brings modified coils back to python from kicad-files.
    Inputs:
        file_directory: string
            Directory to the kicad-files
        layers: (1,N-layers) string array
            On which layers this coil is in kicad-file, e.g., ['F.Cu', 'In1.Cu']
        orig_plane: (N_loops, 3) array
            Points of the original plane, from which kicad-file was originally made from
            (This is needed to shift the plane to right place after importing it back from kicad)
        scale_factor: float
            As KiCad handles everything in mm, the wanted scale must be adjusted. Default value is 10e2.

    Outputs:
        new_points_starts: starts of the line segments in new plane imported from kicad including regularization of the place sifting
        new_points_ends: end of the line segments in new plane imported from kicad including regularization of the place sifting

        (Starts and ends help calculation of bfield, as linesegments drawn in kicad go to the end of the segment lists.)
    '''

    import numpy as np
    points = orig_plane
    stand_off = points[0][1]
    seg_start = np.zeros((1,3), dtype=np.float64)
    seg_end = np.zeros((1,3), dtype=np.float64)
    filedir = file_directory
    with open(filedir, 'r') as infile:
        k = 0
        for line in infile:
            line = line.replace( ')', '')
            #line = line.replace('(', '')
            txt = line.split()

            if len(txt) > 0 and txt[0] == '(segment' and (txt[-3] == layers[0]  or txt[-3] == layers[1]):
                seg_start[:,1] = stand_off
                seg_end[:,1] = stand_off

                seg_start[k,0] = float(txt[2])/scale_factor
                seg_end[k,0] = float(txt[5])/scale_factor
                seg_start[k,2] = float(txt[3])/scale_factor
                seg_end[k,2] = float(txt[6])/scale_factor


                seg_start = np.append(seg_start, [[0,0,0]], axis = 0)
                seg_end = np.append(seg_end, [[0,0,0]], axis = 0)

                k += 1

        seg_start = seg_start[:-1]
        seg_end = seg_end[:-1]

        #Find the amount of shift of the points by matching the top-left corners of the initial planes and modified planes
        sum_orig = np.sum(points, axis = 1)
        corner_ind_orig = np.where(sum_orig == sum_orig.max())
        sum_seg = np.sum(seg_start, axis = 1)
        corner_ind_new = np.where(sum_seg == sum_seg.max())

        er = seg_start[corner_ind_new] - points[corner_ind_orig]

        new_seg_starts = seg_start - er
        new_seg_ends = seg_end - er

        return new_seg_starts, new_seg_ends



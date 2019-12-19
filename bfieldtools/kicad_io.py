
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


def kicad_to_python(file_directory, layers, stand_off, stack_spacing, stack, plane_axes, scale_factor = 10e2, origin = (0,0), net = None):

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
        stack_spcing: float
            spacing between the leyers
        stack: int
            the index of the stack
        plane_axes: tuple with length 2
            specifies the X- and Y-dimensions used for the PCB
        stack_spacing: float
            specifies the distance between the stacks
        scale_factor: float
            As KiCad handles everything in mm, the wanted scale must be adjusted. Default value is 10e2.
        origin: (2, ) array-like
            Cancelling the origin shift done when coils were converted to kicad-file. Default (0,0,0)
        Net: Int
            Defines the net-class from which the segments are searched. Default None

    Outputs:
        vert_start: starts of the line segments in new plane imported from kicad
        vert_end: end of the line segments in new plane imported from kicad

        (Starts and ends help calculation of bfield, as linesegments drawn in kicad go to the end of the segment lists.)
    '''

    import numpy as np
    seg_start = np.zeros((1,3), dtype=np.float64)
    seg_end = np.zeros((1,3), dtype=np.float64)
    filedir = file_directory
    norm_dir = np.array([1,2,3])[np.logical_not(np.isin([1,2,3], plane_axes))][0]
    with open(filedir, 'r') as infile:
        k = 0
        for line in infile:
            line = line.replace( ')', '')
            txt = line.split()
            if len(txt) > 0 and txt[0] == '(segment' and txt[-2] == '(tstamp': #Checks if there is something weird written on the segment lines
                txt = txt[0:-2]
            if not net:
                if len(txt) > 0 and txt[0] == '(segment' and (txt[-3] == layers[0]  or txt[-3] == layers[1]):
                    if txt[-3] == layers[1]:
                        seg_start[k,norm_dir] =  stand_off + stack*stack_spacing + 0.5*stack_spacing #for back layer of the stack the spacing is 1.5 timest the defined stack_spacing + standoff
                        seg_end[k,norm_dir] =   stand_off + stack*stack_spacing + 0.5*stack_spacing
                    else:
                        seg_start[k,norm_dir] = stand_off + stack*stack_spacing   #for back layer of the stack the spacing is stand_off + stack_spacing
                        seg_end[k,norm_dir] = stand_off + stack*stack_spacing
    
                    seg_start[k,plane_axes[0]] = float(txt[2])/scale_factor - origin[0]
                    seg_end[k,plane_axes[0]] = float(txt[5])/scale_factor - origin[0]
                    seg_start[k,plane_axes[1]] = float(txt[3])/scale_factor - origin[1]
                    seg_end[k,plane_axes[1]] = float(txt[6])/scale_factor - origin[1]
    
    
                    seg_start = np.append(seg_start, [[0,0,0]], axis = 0)
                    seg_end = np.append(seg_end, [[0,0,0]], axis = 0)
    
                    k += 1
            else:
                if len(txt) > 0 and txt[0] == '(segment' and (txt[-3] == layers[0]  or txt[-3] == layers[1] or txt[-1] == str(net)):
                    if txt[-3] == layers[1]:
                        seg_start[k,norm_dir] =  stand_off + stack*stack_spacing + 0.5*stack_spacing #for back layer of the stack the spacing is 1.5 timest the defined stack_spacing + standoff
                        seg_end[k,norm_dir] =   stand_off + stack*stack_spacing + 0.5*stack_spacing
                    else:
                        seg_start[k,norm_dir] = stand_off + stack*stack_spacing   #for back layer of the stack the spacing is stand_off + stack_spacing
                        seg_end[k,norm_dir] = stand_off + stack*stack_spacing
    
                    seg_start[k,plane_axes[0]] = float(txt[2])/scale_factor - origin[0]
                    seg_end[k,plane_axes[0]] = float(txt[5])/scale_factor - origin[0]
                    seg_start[k,plane_axes[1]] = float(txt[3])/scale_factor - origin[1]
                    seg_end[k,plane_axes[1]] = float(txt[6])/scale_factor - origin[1]
    
    
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




def segment_looper(seg_starts, seg_ends, start_ind, up_low):
    """
    Author: Elias Lius
    Arranges the current segments according to the current direction.
    Inputs:
        seg_starts: (N, 3) array
            Start points of the segments
        seg_ends: (N, 3) array
            End points of the segments
        start_ind: int
            Index of the starting segment of the coil
        up_low: int
            Determines wether the segments are from upper or lower coil
            (1 == upper, 0 == lower)
    Outputs:
        match_order: list
            List of indices defining the correct order of the segments
        no_match: list
            List of segments with no matches with current specifications
    """
    
    
    match_order = [] #initialize match_order-list
    no_match = []
    indices = np.arange(len(seg_starts)) #initialize the list of indices of the segments
    match_order.append(start_ind) #add the starting index to mathc_order
    indices = np.delete(indices, start_ind) #delete the used index from indices to prevent using the index again
    
    matcher = seg_ends[start_ind] #The point to which we will next find a matching point
    matcher_layer = 'front' #indicates is the point on outer or inner layer of the stack (e.g. F.Cu or In1.Cu)
    
    l_spacings = np.sort(np.abs(np.unique(seg_starts[:,2]))) #Coordinate of the two layers towards the plane norm
    
    #Next we loop over all the remaining segments and find their correct order
    #We try to find the correct order by matching starting or ending points of the segment to previous segments ending point (ending point meaning the correct ending point according to current direction).
    #Note that some segments may be drawn in wrong direction in kiCAD and hence we need to try to find the match from segments' ending and starting points.
    #As we also have two layers for each stack, we need to take that into consideration aswell. 
    while len(indices):
        #First we divide the remaining segments to inner and outer layers, i.e., which segments are on front layer and which on back
        match_list = seg_starts[indices]
        if up_low == 1:
            front_inds = indices[np.where(match_list[:, 2] == l_spacings[0])[0]] #indices for segments on front layer
            back_inds = indices[np.where(match_list[:, 2] == l_spacings[1])[0]] #indices for segments on back layer
        if up_low == 0:
            front_inds = indices[np.where(match_list[:, 2] == -l_spacings[0])[0]] 
            back_inds = indices[np.where(match_list[:, 2] == -l_spacings[1])[0]] 
    
        
        #Next we will wind the closest matches from front and back layers and see which one is closer, and add that segment's index to match_order next.
        #We use try-except condition to avoid errors in case of empty front_inds or back_inds.
        try:
            min_s_f = np.min(np.linalg.norm((seg_starts[front_inds,0:2] - matcher[0:2]),axis=1)) #Min to matcher distance from starting points on front layer
            min_e_f = np.min(np.linalg.norm((seg_ends[front_inds,0:2] - matcher[0:2]),axis=1)) #Min to matcher disntance from ending points on front layer
            
            min_front = np.min(np.array([min_s_f, min_e_f])) #See which one has smaller distance to current matcher-point
    
            
            #here we select the correct index of the point with minimum distance to the mathcer as match_ind_f (match_index_front)
            if min_s_f == min_front:
                match_ind_f = front_inds[(np.linalg.norm((seg_starts[front_inds,0:2] - matcher[0:2]),axis=1)).argmin()]
                matcher_front = seg_ends[match_ind_f]
                
            if min_e_f == min_front:
                match_ind_f = front_inds[(np.linalg.norm((seg_ends[front_inds,0:2] - matcher[0:2]),axis=1)).argmin()]
                matcher_front = seg_starts[match_ind_f]
    
        #if error occurs we select min distance on front layer to be infinity
        except:
            min_front = np.inf
            
        
        #Here we do the same thing for the back layer
        try:
            min_s_b = np.min(np.linalg.norm(seg_starts[back_inds,0:2] - matcher[0:2],axis=1))
            min_e_b = np.min(np.linalg.norm(seg_ends[back_inds,0:2] - matcher[0:2],axis=1))
            
            min_back = np.min(np.array([min_s_b, min_e_b]))
    
            #Selection of the correct index of the point with smallest distance to matcher on back-layer as match_ind_b (match_index_back)
            if min_s_b == min_back:
                match_ind_b = back_inds[(np.linalg.norm(seg_starts[back_inds,0:2] - matcher[0:2],axis=1)).argmin()]
                matcher_back = seg_ends[match_ind_b]
                
            if min_e_b == min_back:
                match_ind_b = back_inds[(np.linalg.norm(seg_ends[back_inds,0:2] - matcher[0:2],axis=1)).argmin()]
                matcher_back = seg_starts[match_ind_b]
            
        #if error occurs we select min distance on back layer to be infinity     
        except:
            min_back = np.inf
    
    
    
    
        #in case of equal distances of the closest points on front and back layer, we select the point which is on the same layer as the matcher
        if min_front<0.000127 and min_back<0.000127:
            if matcher_layer == 'front':
                matcher = matcher_front
                match_order.append(match_ind_f)
                indices = np.delete(indices, np.where(indices == match_ind_f))
                matcher_layer = 'front'
            if matcher_layer == 'back':
                matcher = matcher_back
                match_order.append(match_ind_b)
                indices = np.delete(indices, np.where(indices == match_ind_b))
                matcher_layer = 'back'
                
        #If both distances are over the limit, it might be due to layer switching by via. Thus, we check if the closest point to matcher lies on the other layer than matcher itself, and if so, accept
        #it as a match. Otherwise, we move both closest points (front and back) to the no_match, list.
        elif min_front > 0.000127 and min_back > 0.000127:
            if min_front < min_back and matcher_layer == 'back':
                matcher = matcher_front
                match_order.append(match_ind_f)
                indices = np.delete(indices, np.where(indices == match_ind_f))
                matcher_layer = 'front'
            elif min_back < min_front and matcher_layer == 'front':
                matcher = matcher_back
                match_order.append(match_ind_b)
                indices = np.delete(indices, np.where(indices == match_ind_b))
                matcher_layer = 'back'                
            else:
                no_match.append(match_ind_b)
                indices = np.delete(indices, np.where(indices == match_ind_b))
                no_match.append(match_ind_f)
                indices = np.delete(indices, np.where(indices == match_ind_f))
            
        
        #Now we can compare, which layer had the closest match to the matcher and select that point's index as the next value of match_order
        elif min_front < min_back:
            if min_front > 0.000127:
                no_match.append(match_ind_f)
                indices = np.delete(indices, np.where(indices == match_ind_f))
            else:
                matcher = matcher_front
                match_order.append(match_ind_f)
                indices = np.delete(indices, np.where(indices == match_ind_f))
                matcher_layer = 'front'
    
            
        elif min_front > min_back:
            if min_back > 0.000127:
                no_match.append(match_ind_b)
                indices = np.delete(indices, np.where(indices == match_ind_b))
            else:
                matcher = matcher_back
                match_order.append(match_ind_b)
                indices = np.delete(indices, np.where(indices == match_ind_b))
                matcher_layer = 'back'

    return match_order, no_match






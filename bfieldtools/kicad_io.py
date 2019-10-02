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
        
        

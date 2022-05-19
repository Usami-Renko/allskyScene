# from __future__ import print_function

cimport cython
# from libc.math cimport sin, cos, asin, acos, tan, atan
cimport numpy as cnp
import numpy as np
from libc.math cimport abs

@cython.boundscheck(False)
@cython.wraparound(False)
cpdef monotonicalize(float[:,:] lon, float[:,:] lat, float boundmin, float boundmax, float fillvalue):
    '''
    lon (nlines, nscans): longitude
    boundmin: smaller side of cyclic topology, in Geography is Longitude -180.
    boundmax: larger side of cyclic topology, in Geography is Longitude 180.
    '''
    cdef int nlines = lon.shape[0]
    cdef int nscans = lon.shape[1]
    cdef float offset
    cdef float thresh = 5.0
    cdef float delta = 1e-6
    cdef float thislon
    cdef float nextlon
    cdef int iline
    cdef int iscan 
    cdef cnp.ndarray Templine = np.zeros([nscans], dtype=np.float32)
    cdef cnp.ndarray Tempoffset = np.zeros([nlines], dtype=np.float32)

    # 1. deal with fillvalues
    # Known BUG: case when an area of GeoLocation is missing
    for iline in range(nlines):
        for iscan in range(nscans):
            if (abs(lon[iline, iscan] - fillvalue) < delta) | (abs(lat[iline, iscan] - fillvalue) < delta):
                if iline == 0:
                    lon[iline, iscan] = lon[iline+1,iscan] + delta
                    lat[iline, iscan] = lat[iline+1,iscan] + delta
                elif iline == nlines-1:
                    lon[iline, iscan] = lon[iline-1,iscan] + delta
                    lat[iline, iscan] = lat[iline-1,iscan] + delta
                else:
                    lon[iline, iscan] = (lon[iline-1,iscan] + lon[iline+1,iscan])/2.0 + delta
                    lat[iline, iscan] = (lat[iline-1,iscan] + lat[iline+1,iscan])/2.0 + delta

    # 2. monotonicalize every lines
    for iline in range(nlines):
        offset = 0.0
        Templine[0] = lon[iline, 0]
        for iscan in range(nscans-1):
            thislon = lon[iline, iscan]
            nextlon = lon[iline, iscan+1]
            if (boundmax - thislon < thresh) and (nextlon - boundmin < thresh):
                # stride over the right border 180. -> -180.
                #                           thislon -> nextlon
                offset += boundmax - boundmin
            elif (thislon - boundmin < thresh) and (boundmax - nextlon < thresh):
                # stride over the left border -180. -> 180.
                #                           thislon -> nextlon
                offset -= boundmax - boundmin
            Templine[iscan+1] = nextlon + offset 
        for iscan in range(nscans):
            lon[iline, iscan] = Templine[iscan]

    # 3. monotonicalize the lines in a row
    # lon[iline,0] have not changed after step 2, so still in range [-180.,180.] 
    Tempoffset[:] = 0.0
    offset = 0.0
    for iline in range(nlines - 1):
        thislon = lon[iline,   0]
        nextlon = lon[iline+1, 0]
        if (boundmax - thislon < thresh) and (nextlon - boundmin < thresh):
            # stride over the right border 180. -> -180.
            #                           thislon -> nextlon
            offset += boundmax - boundmin
        elif (thislon - boundmin < thresh) and (boundmax - nextlon < thresh):
            # stride over the left border -180. -> 180.
            #                           thislon -> nextlon
            offset -= boundmax - boundmin
        Tempoffset[iline+1] = offset
    
    for iline in range(nlines):
        for iscan in range(nscans):
            lon[iline, iscan] += Tempoffset[iline]


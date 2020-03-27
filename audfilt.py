import numpy as np

def audfilt(rl, ru, asize, sampfreq):
    """Calculate an auditory filter array. Direct translation of internal function audfilt from MSBG hearing loss model."""

    #  rl = broadening factor on the lower side
    #  ru = broadening factor on the upper side
    # asize = 256; # default, changed from size
    asize = int(asize)
    filter = np.zeros((asize,asize));
    filter[0,0] = 1.;
    
    # Dividing by the erb to remove spectral tilt from the excitation pattern
    filter[0,0] = filter[0,0]/((rl+ru)/2);
    g = np.zeros(asize);
    
    for i in np.linspace(1,asize-1,asize-1,dtype=int):
         fhz = i * np.divide(sampfreq,(2 * asize))
         erbhz = 24.7 * ((fhz * .00437) + 1.0)
         pl = 4.0 * fhz/(erbhz*rl)
         pu = 4.0 * fhz/(erbhz*ru)
         jj = np.arange(0,i,1,dtype=int) # generate indices changed j to jj
          # for lower side of the filter
         g[jj] = abs((i-jj)/i)
         filter[i,jj] = (1+(pl * g[jj])) * np.exp(-1.0 * pl * g[jj])
         jj = np.arange(i,asize);                            # for upper side of the filter and center
         g[jj] = abs((i-jj)/i)
         filter[i,jj] = (1+(pu * g[jj])) * np.exp(-1.0 * pu * g[jj])
         filter[i,] = np.divide(filter[i,], (erbhz * (rl+ru)/(2 * 24.7)))
    
    return filter

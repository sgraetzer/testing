# to track a certain percentage of frames in order to get measure of rms
# derivative of generate_key
# wavlevel.m to measure the rms of a 16 kHz monaural wav file and return mean and stddev
# adaptively sets threshold after looking at histogram of whole recording
# MAS 25-9-00
# NB threshold is in dB
import numpy as np
import math
    
def generate_key_percent(sig, thr_dB, winlen):
    """Generate key percent. Direct translation of internal function generate_key_percent from MSBG hearing loss model."""
    
    winlen = int(winlen)
    sig = sig.flatten()
    if winlen != math.floor(winlen): # whoops on fractional indexing: 7-March 2002
        winlen = math.floor(winlen)
        print('\nGenerate_key_percent:\tWindow length must be integer: now', winlen)    
    
    siglen = len(sig)
    
    # new Dec 2003. Possibly track percentage of frames rather than fixed threshold
    if np.size(thr_dB) > 1:
        track_percent = 1
        percent2track = thr_dB[2]
        expected = thr_dB[0]
        print('\nGenerate_key_percent:\ttracking %.1f percentage of frames ' % percent2track)
    else:
        track_percent = 0
        expected = thr_dB
        print('\nGenerate_key_percent:\ttracking fixed threshold')
    
    non_zero = np.power(10,(expected-30)/10)  # put floor into histogram distribution

    nframes = -1;
    totframes = math.floor(siglen/winlen)
    every_dB = np.zeros((1,totframes)).flatten()

    for ix in np.arange(0,winlen*totframes-1,winlen):
        nframes = nframes + 1
        this_sum = np.sum(np.power(sig[ix:(ix+winlen-1)],2)) # sum of squares
        every_dB[nframes] = 10*np.log10(non_zero + this_sum/winlen)

    nframes = nframes + 1
    every_dB = every_dB[:nframes-1] # from now on save only those analysed
    # error here
    nbins, lvls = np.histogram(every_dB[0:nframes-1],139) # Bec 2003, was 100 to give about a 0.5 dB quantising of levels
    
    if track_percent: # new 1-Dec-2003
        inactive_bins = (100-percent2track)*nframes/100 # min number of bins to use
        nlvls = len(lvls)
        inactive_ix = 0; ixcnt = 0
        for ix in np.arange(0,nlvls,1):
            inactive_ix = inactive_ix + nbins[ix]
            if inactive_ix > inactive_bins: 
                break
            else: 
                ixcnt = ixcnt + 1
        if ix == 1:
            print('\nGenerate_key_percent:\tCounted every bin.........')
        elif ix == nlvls:
            print('\nGenerate_key_percent:\tErrrrr, no levels to count')
        # Improve error message here
        expected = lvls[max(1,ixcnt)]
        
    # set new threshold conservatively to include more bins than desired (rather than fewer)
    used_thr_dB = expected

    # histogram should produce a two-peaked curve: thresh should  be set in valley
    # between the two peaks, and set threshold a bit above that, as it heads for main peak
    frame_index = np.nonzero(every_dB >= expected)
    if np.size(frame_index,0) > 1:
        frame_index = frame_index[1]
    else:
        frame_index = frame_index[0]
    valid_frames = len(frame_index)
    key = np.zeros((1,valid_frames*winlen))[0]
    
    # Up to here
    # convert frame numbers into indices for sig   
    for ix in np.arange(0,valid_frames,1):
        meas_span = np.arange(((frame_index[ix])*winlen),(frame_index[ix]+1)*winlen,1)
        key_span  = np.arange(((ix)*winlen),(ix+1)*winlen,1)
#        if min(key_span) < 1:
#            print('\n\t\tkey_span: Trapped erroneous indexing %d:%d: PAUSED' % 1+((ix-1)*winlen) % ix*winlen)
#            pause
        key[key_span] = meas_span
        key = key.flatten()

    return key, used_thr_dB



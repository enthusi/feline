#
# FILE: masklib.py
# AUTHOR: C. Herenz
# DATE: 2011-07
# INFO: convinience functions for masking (currently used by mask3D & signal_to_noise_mask.py

import numpy as np

def grow_cycle(mask):
    """
    one grow cycle - every filled pixel with three non-empty neighbours gets set to zero
    """
    shape1=mask.shape[0]
    shape2=mask.shape[1]
    dmask = np.ones_like(mask)
    # brute force neighbour testing:
    for i in xrange(shape1):
        for j in xrange(shape2):
            k=0
            # same topology as in conways game of life
            if mask[i,j] == 1:
                if i>0 and j>0:
                    if mask[i-1,j-1] == 0:
                        k+=1
                if i>0:
                    if mask[i-1,j] == 0:
                        k+=1
                if i>0 and j+1<shape2:
                    if mask[i-1,j+1] == 0:
                        k+=1
                if j+1<shape2:
                    if mask[i,j+1] == 0:
                        k+=1
                if i+1<shape1 and j+1<shape2:
                    if mask[i+1,j+1] == 0:
                        k+=1
                if i+1<shape1:
                    if mask[i+1,j] == 0:
                        k+=1
                if i+1<shape1 and j>0:
                    if mask[i+1,j-1] == 0:
                        k+=1
                if j>0:
                    if mask[i,j-1] == 0:
                        k+=1
            if k >= 3: # has more than three neighbors -> grow
                dmask[i,j]=0
    return dmask 

def shrink_cycle(mask):
    """
    one shrink cycle - every full pixel with three empty neighbours gets  set to zero
    """
    shape1=mask.shape[0]
    shape2=mask.shape[1]
    dmask = np.zeros_like(mask)
    # brute force neighbour testing:
    for i in xrange(shape1):
        for j in xrange(shape2):
            k=0
            # same topology as in conways game of life
            if mask[i,j] == 0:
                if i>0 and j>0:
                    if mask[i-1,j-1] == 1:
                        k+=1
                if i>0:
                    if mask[i-1,j] == 1:
                        k+=1
                if i>0 and j+1<shape2:
                    if mask[i-1,j+1] == 1:
                        k+=1
                if j+1<shape2:
                    if mask[i,j+1] == 1:
                        k+=1
                if i+1<shape1 and j+1<shape2:
                    if mask[i+1,j+1] == 1:
                        k+=1
                if i+1<shape1:
                    if mask[i+1,j] == 1:
                        k+=1
                if i+1<shape1 and j>0:
                    if mask[i+1,j-1] == 1:
                        k+=1
                if j>0:
                    if mask[i,j-1] == 1:
                        k+=1
            if k >= 3: # has more than three empty neighbors -> shrink
                dmask[i,j] = 1                
#    mask = dmask+mask # new mask ... shrinked one cycle
    return dmask


def evolve(Ncycles,mask):
    """
    evolve the <mask> in <Ncycles>
    """
    # After some experimentation I found that the
    # following method works best :
    # - define number of total cycles
    # - 2 cycles will be shrink cycles - the second one and the last one
    if Ncycles > 2:
        Nshrink = 2
        # - grow cylces are the rest
        Ngrow = Ncycles - Nshrink
    else:
        Nshrink = 0
        Ngrow = Ncycles
    print('Working on the Mask with '+str(Nshrink)+' shrink cycles and '+str(Ngrow)+' grow cycles')
    shrinks=0 ; grows=0
    for i in xrange(Ncycles):
        if grows < Ngrow:
            print('Grow cylcle: '+str(grows+1))
            mask *= grow_cycle(mask)  # just multiply for growing some weed (1*0 = 0)
            grows += 1 
        if shrinks < Nshrink-1:
            print('Shrink cylcle: '+str(shrinks+1))
            mask += shrink_cycle(mask) # just add for shrink (0+1 = 1)
            shrinks += 1
    print('Shrink cylcle: '+str(shrinks+1)+' (last)')
    mask += shrink_cycle(mask) # last shrink cycle

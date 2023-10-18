import matplotlib as mpl
mpl.use('TkAgg')
import numpy as np
import scipy
import scipy.ndimage as ndimage
import scipy.ndimage.filters as filters
import matplotlib.pyplot as plt
import array
import sys
import pyfits

import struct
from astropy import wcs
from astropy.io import fits
from math import *
import os
import pylab as p
from skimage.feature import peak_local_max
#from skyimage.morphology import local_maxima
from os.path import isfile

from scipy.misc import imsave
from scipy.misc import imread
fig=plt.figure()
manuals=[]

#qsotag=sys.argv[4].ljust(14)
#print "#looking for tag: '%s'" % qsotag

def onclick(event):
    print('#button=%d, x=%d, y=%d, xdata=%f, ydata=%f' %
          (event.button, event.x, event.y, event.xdata, event.ydata))
    
    plt.plot(event.xdata,event.ydata, 'yo')
    manuals.append(( event.ydata, event.xdata))
    plt.show()
    

def world_to_pix(coord, rad):

        #print pix

        radarray = np.array([[rad[0], rad[1], 0]], np.float_)
        #print pixarray
        world = coord.wcs_world2pix(radarray, 0)
        x = world[0][0]
        y = world[0][1]
        return x, y

    
def pix_to_world(coord, pix):

        #print pix

        pixarray = np.array([[pix[0], pix[1], 0]], np.float_)
        #print pixarray
        world = coord.wcs_pix2world(pixarray, 0)
        ra = world[0][0]
        dec = world[0][1]
        return ra, dec


          
#yd=461#322#329
#xd=454#323#331  
#if len(sys.argv)<4:
  #print "Syntax: %s s2ncube threshold neighborhood_size" % sys.argv[0]
  #sys.exit(1)
  
#neighborhood_size = 7
#threshold = 80


hdu = fits.open(sys.argv[1])
coord = wcs.WCS(hdu[0].header)

#neighborhood_size = int(sys.argv[2])
#threshold = int(sys.argv[3])

#print "#using threshold: %d and neighborhood_size: %d" % (threshold,neighborhood_size)

with open("raw_reordered_s2ncube.dat", 'rb') as f:
    header=f.read()[:16]
    
dz= struct.unpack('f', header[0:4])[0]
xd= struct.unpack('f', header[4:8])[0]
yd= struct.unpack('f', header[8:12])[0]

xd=int(xd)
yd=int(yd)
dz=int(dz)
print "#Cube dimensions (z,y,x): %d, %d, %d" % (dz,xd,yd)

#for i in header:
  #print ord(i)
#sys.exit(0)

atoms=[[6564.61],[4862.72],[4341.68],[4102.89],[3727.09,3729.88],[4960.30,5008.24],\
[6549.86,6585.27],[6718.29,6732.67],[3869.81,3968.53],[1908.10,1906.05],[1215.67]]


atom_id={}
atom_id[6564.61]="Ha"
atom_id[4862.72]="Hb"
atom_id[4341.68]="Hg"
atom_id[4102.89]="Hd"

atom_id[3727.09]="OIIa"
atom_id[3729.88]="OIIb"

atom_id[4960.30]="[OIII]a"
atom_id[5008.24]="[OIII]b"

atom_id[6549.86]="[NII]a"
atom_id[6585.27]="[NII]b"

atom_id[6718.29]="[SII]a"
atom_id[6732.67]="[SII]b"

atom_id[3869.81]="NeIIIa"
atom_id[3968.53]="NeIIIb"

atom_id[1908.10]="CIIIa"
atom_id[1906.05]="CIIIb"

atom_id[1215.45]="Lya"
mpl.rcParams["savefig.directory"] = "."

def gauss2d(xy, amp, x0, y0, a, b, c):
    x, y = xy
    inner = a * (x - x0)**2
    inner += 2 * b * (x - x0)**2 * (y - y0)**2
    inner += c * (y - y0)**2
    return amp * np.exp(-inner)


def twoD_Gaussian((x, y), amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
    return g.ravel()
  
def print_lines(toggle,z):
  for k in range(len(atoms)):
        #is k in the template?
        if toggle&0x1==0:
          toggle=toggle/2
          continue
        
        #ok, we consider this atom/transition
        toggle=toggle/2
        atom=atoms[k]
        
        for emission in atom:
          pos=emission*(z+1)
          name=atom_id[emission]
          print "%s (%.1f)," % (name,pos),
  print

          
isize=xd*yd
size=isize
        
data=np.fromfile('float32_array_omp4.raw',dtype='float32')
plane,redshift,template,used=np.split(data,4)



plane.resize((xd,yd))
redshift.resize((xd,yd))
template.resize((xd,yd))
used.resize((xd,yd))

imsave("image.png", plane)
if isfile("imagemask.png"):
    mymask=imread("imagemask.png")

    ny,nx= mymask.shape
    for iy in range(ny):
        for ix in range(nx):
            if mymask[iy,ix]==0xff:
                plane[iy,ix]=0

else:
    pass
    #print 'CAUTION! no valid imagemask.png defined yet!'
    
            

data=ndimage.gaussian_filter(plane, sigma=1)




#data=plane
coordinates = peak_local_max(data, min_distance=1,threshold_abs=50,exclude_border=10,num_peaks=300)
#print len(coordinates)
#sys.exit(1)
#remove border
#for yi in range(yd-1):
#    for xi in range(xd):
#        print xd,yd,xi,yi,plane[xi,yi]
#        #data[xi,yi]=0.0
#        data[1:-1,1:-1] = 0

#print len(xy[:,0])


#remove borders from results
width=0
bleft=width
btop=width
bright=xd-width
bbottom=yd-width

hdu_muse = fits.open("image00.fits",
                     memmap=False)  
data_muse = hdu_muse[1].data
nan_sel = np.isnan(data_muse)

xy=[]
ay=[]
ax=[]
for val in coordinates:
    nx,ny =val
    #print nx,ny
    #plt.plot(ny,nx,'go')
    ax.append(nx)
    ay.append(ny)
    xy.append((nx,ny))
    
print "#",len(xy) 
    

          

  
test=plt.imshow(data,vmin=30,vmax=500, interpolation='nearest',cmap='jet')

plt.colorbar()
#plt.savefig('data.png', bbox_inches = 'tight')

plt.autoscale(False)
#plt.plot(xy[:, 1], xy[:, 0], 'ro')

#plot other algo

#plt.plot(exclude_xy[:, 1], exclude_xy[:, 0], 'w.')
#read in lorrie
#fitscat=[]
#catfilename="../lowz_sources_v1.2.fits"
#hdulist= pyfits.open(catfilename,ext=1)
#print hdulist.info()
#prihdr = hdulist[1].data

#ra=prihdr[(prihdr['QSO']==qsotag)]['RA']
#dec=prihdr[(prihdr['QSO']==qsotag)]['DEC']
#print ra,prihdr['QSO']
#try:
#    x=prihdr[(prihdr['QSO']==qsotag)]['X']
#    y=prihdr[(prihdr['QSO']==qsotag)]['Y']
#except:
#    print "#computing x,y from ra,dec",  len(ra)
#    x=[]
#    y=[]
#    for j in range(len(ra)):
#        xj,yj=world_to_pix(coord, (ra[j],dec[j]))
#        tr,td=pix_to_world(coord,(xj,yj))
#        #print ra[j],dec[j],xj,yj,tr,td
#        x.append(xj)
#        y.append(yj)
#        
#qop=prihdr[(prihdr['QSO']==qsotag)]['QOP_MARZ']
#z=prihdr[(prihdr['QSO']==qsotag)]['Z_MUSE']
#zplate=prihdr[(prihdr['QSO']==qsotag)]['Z_PLATEFIT']
#b=prihdr[(prihdr['QSO']==qsotag)]['b_kpc']
#infile=open('sextractor_sources_v0.8_t1_crop.txt','r')

#for i in range(len(x)):
# #   print x[i],y[i],ra[i],dec[i],qop[i],zplate[i],z[i],b[i]
#    
##for row in infile:
#    lx=x[i]
#    ly=y[i]
#    lz=zplate[i]
#    qual=qop[i]
#    if lz>1.5: continue
#    if lz<0.01: continue
#    #print qual
#    if qual<0.9: continue
#    fitscat.append((ly,lx))
#    plt.plot(lx,ly,'y*')
#%https://matplotlib.org/users/event_handling.html
#cid = fig.canvas.mpl_connect('button_press_event', onclick)

run_id=0
for hit in xy: #y,x
  #continue
  y,x=hit
  if x<1: continue
  y=int(y)
  x=int(x)
  u_i=int(used[y,x])
  z_i=redshift[y,x]
  q_i=plane[y,x]
  t_i=int(template[y,x])
  
  print "%d %d %d %1.6f %d %d %d" % (run_id,int(y),int(x),z_i,q_i,u_i,t_i),
  ra, dec = pix_to_world(coord, (x,y))
  print "\t%.6f %.6f" % (ra,dec),
  print_lines(t_i,z_i)
  run_id+=1
  
#now add the fake sources
"""
catalog="sources_fake1_vconsol.csv"
refcat=open(catalog)
refy=[]
refx=[]
for entry in refcat:
            continue
            if '#' in entry: continue
            if 'xpos' in entry: continue
        
            x=float(entry.split(',')[0])
            y=float(entry.split(',')[1])
            x=int(x)
            y=int(y)                
            #rest=3728.85169431
            #wl=float(entry.split(',')[2])
            #rz=wl/rest-1.0
            u_i=int(used[y,x])
            z_i=redshift[y,x]
            q_i=plane[y,x]
            t_i=int(template[y,x])
            print "%d %d %d %1.6f %d %d %d" % (run_id,int(y),int(x),z_i,q_i,u_i,t_i),
            ra, dec = pix_to_world(coord, (x,y))
            print "\t%.6f %.6f" % (ra,dec),
            print_lines(t_i,z_i)
            run_id+=1

            refx.append(y)
            refy.append(x)
            
plt.plot(refy,refx,'gx',markersize=4)
"""    


plt.plot(ay,ax,'rx',markersize=2)

plt.title("%d sources " % (len(xy))   )
plt.show()
plt.savefig('result.png', bbox_inches = 'tight')

  

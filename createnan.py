from mpdaf.obj import Cube
from mpdaf.obj import Image

import matplotlib.pyplot as plt
from math import *
import sys
import numpy as np

first=True

plot=True

if first:
  name=sys.argv[1]
  #infile=CubeDisk(name)
  #infile.info()
  #image=Cube(name,ext=4)
  #image.write('image00.fits')
  #sys.exit(0)
  data=Cube(name,ext=1).sum(axis=0)
  data.info()
  print "data loaded"
  stat=Cube(name,ext=2).sum(axis=0)
  stat.info()
  print "stat loaded"
  #data.plot()
  #raw_input('Press enter')
  #stat.plot()
  #raw_input('Press enter')
  snimage=data.copy()
  snimage.data=data.data/stat.data
  snimage.info()
  snimage.write('image00.fits')
  
else: 
  snimage=Image('image00.fits')
  

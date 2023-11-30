from mpdaf.obj import Cube
import sys
import struct
import numpy.ma
import matplotlib.pyplot as plt

zqso=1.23
pre_select_sn = 2
ignore_below=1.5 #in the s2n cube
significance=2.5 #lines*scale*significance
scale=10
#samples=400
zstep=0.0004#1/3000.0
lmin=4750
lmax=9250
step=1.25
z=0.0
skipped=0

cubefile=sys.argv[1]
c0=Cube(cubefile)
dz,dy,dx=c0.shape
start=c0.wave.get_crval()
print dz,dy,dx,start

fout = open('raw_reordered_s2ncube.dat', 'wb')
fout.write(struct.pack('f', dz))
fout.write(struct.pack('f', dy))
fout.write(struct.pack('f', dx))
fout.write(struct.pack('f', start))


c1=c0 > pre_select_sn
i1=c1.sum(axis=0)


if sys.argv[2]=="plot":
  i1.plot()
  plt.show()
  raw_input("press any key")


for y in range(dy):
  if not (y%10): print "%.1f done.\r" % (float(y)/float(dy)*100.0)
  for x in range(dx):
    spec=c0.data[:,y,x]
    myfmt='f'*len(spec.data)
    fout.write(struct.pack(myfmt,*spec.data))

print
fout.close()
    
                

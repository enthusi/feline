import os
import sys
import struct
import project_path_config
import mpdaf.obj
import matplotlib.pyplot as plt


filename = sys.argv[1]
file = os.path.join(project_path_config.DATA_PATH_PROCESSED, filename)

c0 = mpdaf.obj.Cube(file)
dz, dy, dx = c0.shape
start = c0.wave.get_crval()
print(f"{dz}, {dy}, {dx}, {start}")

with open(os.path.join(project_path_config.DATA_PATH_PROCESSED, "raw_reordered_s2ncube.dat"), "wb") as fout:
	fout.write(struct.pack('f', dz))
	fout.write(struct.pack('f', dy))
	fout.write(struct.pack('f', dx))
	fout.write(struct.pack('f', start))

	for y in range(dy):
		if not (y % 10):
			print(f"{float(y) / float(dy) * 100:.1f} done.\r")
		for x in range(dx):
			spec = c0.data[:, y, x]
			myfmt = "f" * len(spec.data)
			fout.write(struct.pack(myfmt, *spec.data))

# brauchen wir das noch?
pre_select_sn = 2

c1 = c0 > pre_select_sn
i1 = c1.sum(axis=0)

if sys.argv[2] == "plot":
	i1.plot()
	plt.show()
	input("press any key")

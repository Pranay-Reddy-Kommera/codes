
import pylab as plt
import numpy as np

from math import pi


# read data from files
raw_x = np.loadtxt("./data/grid_xs.dat",)
raw_z = np.loadtxt("./data/grid_zs.dat")
raw_ght = np.loadtxt("./data/mtn_grid_xz.dat")
raw_data = np.loadtxt("./data/data_Ycut.dat")

# reshape data
meshx,meshz = np.meshgrid(raw_x,raw_z)
x = meshx.reshape(-1)/1000.0
z = meshz.reshape(-1)/1000.0
d = raw_data
ght = raw_ght/1000.0

# plot
f,ax = plt.subplots(2,1, sharex=True, sharey=True)
levels = np.linspace(0.1, 1.0, 10)

ax[0].tricontourf(x, z, d, levels, alpha=0.75, cmap=plt.cm.Reds)
ax[0].tricontour(x, z, d, levels, colors='k', linewidth=0.33)
ax[0].set_title('computational (x,zeta) grid')
#ax[0].plot(x,z,'k.')

ax[1].tricontourf(x, ght, d, levels, alpha=0.75, cmap=plt.cm.Reds)
ax[1].tricontour(x, ght, d, levels, colors='k', linewidth=0.33)
ax[1].set_title('physical (x,z) grid')
#ax[1].plot(x,ght,'k.')

# agnesi
#mountain = 3.e3 / (((raw_x - max(raw_x)/2)/12.e3)**2 + 1)**1.5
# schar
mountain = 5.e3 * np.exp(-((raw_x - max(raw_x)/2)/10.e3)**2) * (np.cos(pi*raw_x/6.e3))**2
ax[1].fill_between(raw_x/1000, 0, mountain/1000, color='k')

plt.savefig("plot.pdf")


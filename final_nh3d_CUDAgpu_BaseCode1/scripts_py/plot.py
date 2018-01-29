
import pylab as pl
import numpy as np

from math import pi
from scipy.interpolate import RectBivariateSpline


def prepare_orig_data(x,y,data):
    mx,my = np.meshgrid(x,y)
    return mx,my,data

def prepare_spline_data(x,y,data):
    interp = RectBivariateSpline(x,y,data.T)
    newxs = np.linspace(x[0],x[-1],4*x.size)
    newys = np.linspace(y[0],y[-1],4*y.size)
    mx,my = np.meshgrid(newxs,newys)
    newdata = interp(newxs,newys,grid=True).T
    return mx,my,newdata

# streamlines appear to work best when using a UNIFORM underlying grid..
# to avoid any complications, do not try to reuse the data grid
def plot_streamlines(xin,yin):
    x = np.linspace(xin[0], xin[-1], 64)
    y = np.linspace(yin[0], yin[-1], 64)
    X,Y = np.meshgrid(x,y)
    relX = (X - xin[0]) / (xin[-1] - xin[0])
    relY = (Y - xin[0]) / (yin[-1] - yin[0])
    A = 2*pi
    B = pi
    U = 1 - 4 * np.cos(B*relY) * np.sin(A*relX)**2 * np.sin(B*relY)
    V = 4 * np.cos(A*relX) * np.sin(B*relY)**2 * np.sin(A*relX)
    S = np.sqrt(U*U + V*V)
    pl.streamplot(X, Y, U, V,
                  linewidth=np.sqrt(S),
                  color="lightsteelblue",
                  density=1.1,
                  minlength=0.9,
                  arrowstyle='-')


# read data from files
x = np.loadtxt("./data/grid_xs.dat",)
#y = np.loadtxt("./data/grid_ys.dat")
z = np.loadtxt("./data/grid_zs.dat")
data = np.loadtxt("./data/data_Ycut.dat")

# cut raw volume data and prepare for plotting
data = data.reshape(z.size, x.size)
plotX,plotZ,plotD = prepare_spline_data(x, z, data)

# contours
levels = np.linspace(0.1, 1.0, 10)
c1 = pl.contourf(plotX, plotZ, plotD, levels, alpha=.75, cmap=pl.cm.Reds)
c2 = pl.contour(plotX, plotZ, plotD, levels, colors='black', linewidth=.33)
cbar = pl.colorbar(c1, orientation='horizontal')
cbar.add_lines(c2)

# streamlines
#plot_streamlines(x,z)

pl.axes().set_aspect('equal')
pl.savefig("plot.pdf")


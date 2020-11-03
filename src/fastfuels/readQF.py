
"""
Read and plot QUICFIRE input files

lucaswells 9-17-2020
"""

import numpy as np
import pyvista as pv # pip3 install pyvista
import matplotlib.pyplot as plt

def read_input(file_name, dim_xyz):

    nx, ny, nz = dim_xyz
    size = nx*ny*nz
    data = np.fromfile(file_name, dtype=np.float32, count=size)

    # return 3D view of flattened input
    return data.reshape((nz, ny, nx))

def show_3d(data, theta=None):

    data = data.T

    data[data == 0] = -1
    grid = pv.UniformGrid()
    grid.dimensions = np.array(data.shape) + 1
    grid.spacing = (2,2,1)
    grid.cell_arrays['values'] = data.flatten(order='F')

    grid = grid.threshold(0)

    # show it
    pv.set_plot_theme('document')

    # change cmap by importing matplotlib.pyplot and using ppyplot.cm
    grid.plot(cmap=plt.cm.summer)


# example; change fname and dims
data = read_input('../../output/bulk_density.dat', [225, 225, 100])
show_3d(data, 1e-6)

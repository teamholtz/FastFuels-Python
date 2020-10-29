
"""
Model agnostic 3D storage structure for CFD fire model inputs
"""

# built-in


# external
import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt


__author__     = "Lucas Wells"
__copyright__  = "Copyright 2020, Holtz Forestry LLC"
__version__    = "0.0.1"
__maintainer__ = "Lucas Wells"
__email__      = "lucas@holtzforestry.com"
__status__     = "Prototype"


class ParameterArray:

	def __init__(self, dim, res):

		# force integer dim and res elements
		self.dim = dim #[int(i) for i in dim]
		self.res = res #[int(i) for i in res]

		self.dim_x, self.dim_y, self.dim_z = self.dim
		self.res_x, self.res_y, self.res_z = self.res

		self.n_cells_x = int(self.dim_x / self.res_x)
		self.n_cells_y = int(self.dim_y / self.res_y)
		self.n_cells_z = int(self.dim_z / self.res_z)

		self.domain = np.zeros((self.n_cells_y, self.n_cells_x,
			self.n_cells_z), np.float32)

	def _cast(self, pos):

		x, y, z = [float(i) for i in pos]

		x = int(x/self.res_x)
		y = int(y/self.res_y)
		z = int(z/self.res_z)

		return [x, y, z]

	def is_empty(self, pos, nodata=0):

		x, y, z = self._cast(pos)

		val = self.domain[y, x, z]

		if val == nodata:
			return 1
		else:
			return 0

	def in_bounds(self, x, y, z):

		rtn = True

		if not ((x > 0) and (x < self.n_cells_x)):
			rtn = False
		if not ((y > 0) and (y < self.n_cells_y)):
			rtn = False
		if not ((z > 0) and (z < self.n_cells_z)):
			rtn = False

		return rtn

	def insert(self, pos, value):

		x, y, z = self._cast(pos)

		if self.in_bounds(x, y, z):
			self.domain[y, x, z] = value
		else:
			print(f"x: {x}, y: {y}, z: {z} out of bounds")

class View:

	def __init__(self, param_array):

		self.param_array = param_array
		self.data = param_array.domain

	def slice(self, idx, dim):

		if dim == 'x':
			plt.imshow(self.data[:,idx,:])
			plt.show()
		if dim == 'y':
			plt.imshow(self.data[idx,:,:])
			plt.show()
		if dim == 'z':
			plt.imshow(self.data[:,:,idx])
			plt.show()

	def show3d(self):

		data = self.data
		data[data == 0] = -1
		grid = pv.UniformGrid()
		grid.dimensions = np.array(data.shape) + 1
		grid.spacing = self.param_array.res
		grid.cell_arrays['values'] = data.flatten(order='F')

		grid = grid.threshold(0)
		grid.plot()

if __name__ == '__main__':

	dim = [1000,1000,100]
	res = [2,2,1]

	test = ParameterArray(dim, res)
	for i in range(10000):
		x = np.random.random()*dim[0]
		y = np.random.random()*dim[1]
		z = np.random.random()*dim[2]
		val = np.random.normal(100, 10)
		test.insert([x,y,z], val)
	viewer = View(test)
	viewer.slice(10, 'y')
	viewer.show3d()

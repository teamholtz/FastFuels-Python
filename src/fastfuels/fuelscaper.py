
"""
Generator for canopy fuels
"""

# built-in
import time
import json
import sys
import math

# external
#import cupy as cp
import cv2
import scipy.ndimage as ndimage
import numpy as np
import pyvista as pv
import gdal
sys.path.append('../../../qFIA/src/')
import qfia

# package-level
import domain
import projections

__author__     = "Lucas Wells"
__copyright__  = "Copyright 2020, Holtz Forestry LLC"
__version__    = "0.0.1"
__maintainer__ = "Lucas Wells"
__email__      = "lucas@holtzforestry.com"
__status__     = "Prototype"


""" Temporary: generates species reference json file
ref_sp = {}
data = np.genfromtxt('../../data/ref_sp.csv', delimiter=',')
for i in range(data.shape[0]):
    code = int(data[i,0])
    group = -1
    if not math.isnan(data[i,1]):
        group = int(data[i,1])
    sapadj = -1
    if not math.isnan(data[i,2]):
        sapadj = data[i,2]
    sorh = -1
    if not math.isnan(data[i,3]):
        sorh = int(data[i,3])
    ref_sp[code] = [group, sapadj, sorh]

with open('../../data/ref_sp.txt', 'w') as f:
    json.dump(ref_sp, f)
"""

# Globals

# load species group and sapling adjustment factor database
# idx 0 = species group; idx 1 = sapling adjustment factor
# sapling adjustment applies to diameter range 1.0-4.9 inches
# No estimates for trees less than 1.0 inches
with open('../../data/ref_sp.txt', 'r') as f:
    SPECIES_REF = json.load(f)

# parameter for above ground biomass estimators by sp. group
JENKINS_AGB_PARAMS = {
    '1':  [-2.0339, 2.2592], # cedar/larch
    '2':  [-2.2304, 2.4435], # Douglas fir
    '3':  [-2.5384, 2.4813], # true fir/hemlock
    '4':  [-2.5356, 2.4349], # pine
    '5':  [-2.0773, 2.3323], # spruce
    '6':  [-2.2094, 2.3867], # aspen/alder/cottonwood/willow
    '7':  [-1.9123, 2.3651], # soft maple/birch
    '8':  [-2.4800, 2.4835], # mixed hardwood
    '9':  [-2.0127, 2.4342], # hard maple/oak/hickory/beech
    '10': [-0.7152, 1.7029], # juniper/oak/mesquite
}

# parameters for component ratio estimators
# 0 = hardwood; 1 = softwood
JENKINS_COMPONENT_PARAMS = {
    '0': {
        'foliage' : [-4.0813,  5.8816],
        'roots'    : [-1.6911,  0.8160],
        'bark'    : [-2.0129, -1.6805],
        'wood'    : [-0.3065, -5.4240]
    },
    '1': {
        'foliage' : [-2.9584,  4.4766],
        'roots'    : [-1.5619,  0.6614],
        'bark'    : [-2.0980, -1.1432],
        'wood'    : [-0.3737, -1.8055]
    }
}

# canopy profile beta parameters by Jenkins sp-grp [a,b,c]
BETA_CANOPY_PARAMS = {
    '1':  [1.2405, 1.5580, 0.1286], # cedar/larch
    '2':  [1.2405, 1.5580, 0.1286], # Douglas fir
    '3':  [1.1250, 1.6973, 0.0718], # true fir/hemlock
    '4':  [1.1821, 1.4627, 0.1528], # pine
    '5':  [1.1250, 1.6973, 0.0718], # spruce
    '6':  [2.3, 1.6, 0.2], # aspen/alder/cottonwood/willow
    '7':  [2.3, 1.6, 0.2], # soft maple/birch
    '8':  [2.3, 1.6, 0.2], # mixed hardwood
    '9':  [2.3, 1.6, 0.2], # hard maple/oak/hickory/beech
    '10': [1.75, 1.6, 0.1] # juniper/oak/mesquite
}

# fuel parameters for the Scott and Burgan 40
# [name, loading (tons/ac), SAV (1/ft), ext. MC (percent), bed depth (ft)]
SB40_PARAMS = {
	91:  ['NB1', 0.00, 0.0000, 0.00, 0.00],
	92:  ['NB2', 0.00, 0.0000, 0.00, 0.00],
	93:  ['NB3', 0.00, 0.0000, 0.00, 0.00],
	98:  ['NB8', 0.00, 0.0000, 0.00, 0.00],
	99:  ['NB9', 0.00, 0.0000, 0.00, 0.00],
	101: ['GR1', 0.40, 2054.0, 15.0, 0.4],
	102: ['GR2', 1.10, 1820.0, 15.0, 1.0],
	103: ['GR3', 1.60, 1290.0, 30.0, 2.0],
	104: ['GR4', 2.15, 1826.0, 15.0, 2.0],
	105: ['GR5', 2.90, 1631.0, 40.0, 1.5],
	106: ['GR6', 3.50, 2006.0, 40.0, 1.5],
	107: ['GR7', 6.40, 1834.0, 15.0, 3.0],
	108: ['GR8', 7.80, 1302.0, 30.0, 4.0],
	109: ['GR9', 10.0, 1612.0, 40.0, 5.0],
	121: ['GS1', 1.35, 1832.0, 15.0, 0.9],
	122: ['GS2', 2.10, 1827.0, 15.0, 1.5],
	123: ['GS3', 3.00, 1614.0, 40.0, 1.8],
	124: ['GS4', 12.4, 1674.0, 40.0, 2.1],
	141: ['SH1', 1.70, 1674.0, 15.0, 1.0],
	142: ['SH2', 5.20, 1672.0, 15.0, 1.0],
	143: ['SH3', 6.65, 1371.0, 40.0, 2.4],
	144: ['SH4', 3.40, 1682.0, 30.0, 3.0],
	145: ['SH5', 6.50, 1252.0, 15.0, 6.0],
	146: ['SH6', 4.30, 1144.0, 30.0, 2.0],
	147: ['SH7', 6.90, 1233.0, 15.0, 6.0],
	148: ['SH8', 6.40, 1386.0, 40.0, 3.0],
	149: ['SH9', 13.1, 1378.0, 40.0, 4.4],
	161: ['TU1', 1.30, 1606.0, 20.0, 0.6],
	162: ['TU2', 1.15, 1767.0, 30.0, 1.0],
	163: ['TU3', 2.85, 1611.0, 30.0, 1.3],
	164: ['TU4', 6.50, 2216.0, 12.0, 0.5],
	165: ['TU5', 7.00, 1224.0, 25.0, 1.0],
	181: ['TL1', 1.00, 1716.0, 30.0, 0.2],
	182: ['TL2', 1.40, 1806.0, 25.0, 0.2],
	183: ['TL3', 0.50, 1532.0, 20.0, 0.3],
	184: ['TL4', 0.50, 1568.0, 25.0, 0.4],
	185: ['TL5', 1.15, 1713.0, 25.0, 0.6],
	186: ['TL6', 2.40, 1936.0, 25.0, 0.3],
	187: ['TL7', 0.30, 1229.0, 25.0, 0.4],
	188: ['TL8', 5.80, 1770.0, 35.0, 0.3],
	189: ['TL9', 6.65, 1733.0, 35.0, 0.6],
	201: ['SB1', 1.50, 1653.0, 25.0, 1.0],
	202: ['SB2', 4.50, 1884.0, 25.0, 1.0],
	203: ['SB3', 5.50, 1935.0, 25.0, 1.2],
	204: ['SB4', 5.25, 1907.0, 25.0, 2.7],
}

# convert parameters to metric
for key in SB40_PARAMS.keys():
	load = SB40_PARAMS[key][1]
	savs = SB40_PARAMS[key][2]
	ht = SB40_PARAMS[key][4]
	if load != 0:
		SB40_PARAMS[key][1] = load*0.22417#(load*ht)*0.735468 # tons/ac-ft to kg/m^3
		SB40_PARAMS[key][2] = (1.0/savs)*0.3048 # inverse feet to meters
		SB40_PARAMS[key][4] = ht*0.3048 # feet to meters

SPECIES_GROUP_COLORS = {

}

# load fia plot cn indices
with open('../../data/fia_ids.txt', 'rb') as f:
    FIA_IDX = json.load(f)

# load treelist layer and surface fuel models
FIA_TREELIST = gdal.Open('../../data/treelist.tif')
SB40_FUELS = gdal.Open('../../data/sb40.tif')
ELEV_30 = gdal.Open('../../data/elevation.tif')

geo_matrix = SB40_FUELS.GetGeoTransform()
ULX, ULY = geo_matrix[0], geo_matrix[3]



class Tree:

    _counter = 0

    def __init__(self, sp, dia, ht, cr, cd=0.5):

        Tree._counter += 1
        self.id = Tree._counter

        self.sp = sp
        self.sp_grp = str(SPECIES_REF[sp][0])
        self.dia = dia
        self.ht = ht
        self.cr = cr
        self.cd = cd
        self.ch = ht*cr
        self.cbh = ht-self.ch

        self.weight = self.get_biomass()

        self.init_beta_canopy()
        self.volume = self.get_volume()

        # currently multiplying foliage weight by crown ratio to resolve
        # high bulk densities when crown length is small
        self.bulk_density = (self.cr*self.weight.foliage)/(self.volume*self.cd)

        self.n_cells = int(self.cd*self.volume + 0.5)

    def get_biomass(self):
        """
        Returns a class with biomass components stored as class
        attributes
        """

        # no biomass cutoff
        if self.dia < 2.54:
            return -1

        # sapling adj factor and hardwood/softwood code
        sap_adj = float(SPECIES_REF[self.sp][1])
        hs = str(SPECIES_REF[self.sp][2])

        # get species group and model coefficents
        b0, b1 = JENKINS_AGB_PARAMS[self.sp_grp]

        # predict biomass
        agb = math.exp(b0 + b1*math.log(self.dia))

        # sapling adjustment factor cutoff
        if self.dia < 12.7:
            agb *= sap_adj

        # foliage
        b0, b1 = JENKINS_COMPONENT_PARAMS[hs]['foliage']
        foliage = agb*math.exp(b0 + (b1/self.dia))

        # roots
        b0, b1 = JENKINS_COMPONENT_PARAMS[hs]['roots']
        roots = agb*math.exp(b0 + (b1/self.dia))

        # bark
        b0, b1 = JENKINS_COMPONENT_PARAMS[hs]['bark']
        bark = agb*math.exp(b0 + (b1/self.dia))

        # wood
        b0, b1 = JENKINS_COMPONENT_PARAMS[hs]['wood']
        wood = agb*math.exp(b0 + (b1/self.dia))

        # dynamically create a class to store components
        components = type('Components', (), {})

        # set attributes with components
        components.total = agb
        components.foliage = foliage
        components.roots = roots
        components.bark = bark
        components.wood = wood
        components.branch = agb - foliage - bark - wood

        return components

    def init_beta_canopy(self):

        self.beta_a, self.beta_b, self.beta_c = BETA_CANOPY_PARAMS[self.sp_grp]

        self.beta_norm = np.sqrt(
            2*np.pi)*((self.beta_a**(self.beta_a-0.5)*self.beta_b**(
                self.beta_b-0.5))/((self.beta_a+self.beta_b)**(
                self.beta_a+self.beta_b-0.5)))

    def get_beta_radius(self, z):

        return self.beta_c*(((z**(self.beta_a-1))*((1-z)**(
            self.beta_b-1))))/self.beta_norm

    def get_volume(self):
        """
        Integrate beta canopy to calculate volume
        """

        z = np.linspace(0, 1, 100)

        dz = self.ch/100.0

        volume = 0

        for z_i in z:
            x = self.get_beta_radius(z_i)
            x_scaled = x*self.ch
            volume += np.pi*(x_scaled**2)*dz

        return volume


class Forest:

    def __init__(self, x, y, w, h):

        proj = projections.AlbersEqualArea()
        x,y = proj.project(x, y)
        x = int(np.abs(x - ULX)/30.0 + 0.5)
        y = int(np.abs(y - ULY)/30.0 + 0.5)

        self.w = int(w/30.0)
        self.h = int(h/30.0)

        self.res = [2,2,1]

        self.treelist_nodata = FIA_TREELIST.GetRasterBand(1).GetNoDataValue()
        self.surface_nodata = SB40_FUELS.GetRasterBand(1).GetNoDataValue()

        self.treelist_array = FIA_TREELIST.ReadAsArray(x, y, self.w, self.h)
        self.surface_array = SB40_FUELS.ReadAsArray(x, y, self.w, self.h)

        dem_array = ELEV_30.ReadAsArray(x, y, self.w, self.h)
        self.dem_array = self.interpolate(dem_array)

        max_elev = np.max(self.dem_array)
        self.dim = [self.w*30, self.h*30, 100]# + max_elev]

        # init parameter arrays
        self.tree_id = domain.ParameterArray(self.dim, self.res)
        self.tree_sp = domain.ParameterArray(self.dim, self.res)
        self.bulk_density = domain.ParameterArray(self.dim, self.res)
        self.moisture = domain.ParameterArray(self.dim, self.res)
        self.sav = domain.ParameterArray(self.dim, self.res)
        self.fuel_depth = domain.ParameterArray(self.dim, self.res)

    def interpolate(self, array):

        return cv2.resize(array, (self.w*(30//self.res[0]),
            self.h*(30//self.res[1])), cv2.INTER_LINEAR)

    def sample_plots(self):

        # samples are stored as a list of Tree objects paired with a key
        # denoting the (i,j) of the treelist_array grid referencing the FIA
        # plot sequence number from which the samples are drawn
        samples = {}

        # loop over each column and row in the treelist_array
        for i in range(self.h):
            print(i)
            for j in range(self.w):

                # initialize the Tree object list at the grid index
                samples[(i,j)] = []

                # if the queried fia idx equals the grid nodata value then
                # skip what's below and return to top of loop
                fia_idx = self.treelist_array[i,j]
                if fia_idx == self.treelist_nodata:
                    continue

                # get the plot sequence number and query the FIA database
                plt_cn = FIA_IDX[str(fia_idx)]
                plot = qfia.query(plt_cn)

                # loop down to the trees in subplots
                for subplot in plot:
                    for tree in subplot:

                        # The tree must have a numeric value for diameter and
                        # a status code of 1 indicating a living tree
                        if (tree.dia) and (tree.statuscd == 1):

                            # Fixed radius plot cutoffs; different expansion
                            # factors apply
                            if tree.dia >= 5.0:
                                n = self._n_subplot_trees()
                            else:
                                n = self._n_microplot_trees()

                            # Extract measurements and add to list
                            # n is the number of tree instances after expansion
                            for k in range(n):
                                spcd = str(int(tree.spcd))
                                dia = tree.dia*2.54
                                ht = tree.actualht*0.3048
                                cr = tree.cr/100
                                samples[(i,j)].append(Tree(spcd, dia, ht, cr,
                                    cd=0.5))

        self.samples = samples

    def _n_subplot_trees(self):

        r = np.random.random()
        if r < 0.338:
            n = 2
        else:
            n = 1

        return n

    def _n_microplot_trees(self):

        r = np.random.random()
        if r < 0.672:
            n = 17
        else:
            n = 16

        return n

    def distribute_canopy(self):

        for i in range(self.h):
            print(i)
            for j in range(self.w):
                for tree in self.samples[(i,j)]:

                    placed = False
                    while not placed:
                        self.tree_id.reset_track()
                        theta = np.random.random()*2*np.pi
                        dist = np.random.random()*30
                        x = np.cos(theta)*dist + 15
                        y = np.sin(theta)*dist + 15
                        x += i*30
                        y += j*30
                        if not ((x > 0) and (x < self.w*30) and (y > 0) and (y < self.h*30)):
                            continue
                        z_m = self.dem_array[int(y/2), int(x/2)]
                        fact = self.res[0]*self.res[1]*self.res[2]
                        n = int(tree.cd*(tree.volume/fact) + 0.5)

                        attemps = 0
                        while n > 0:
                            z = np.random.beta(3,1.5)
                            r = np.random.beta(2.5,1.5)
                            rad = tree.get_beta_radius(z)*r*tree.ch
                            c_theta = np.random.random()*2*np.pi

                            c_x = x + np.cos(c_theta)*rad
                            c_y = y + np.sin(c_theta)*rad
                            c_z = tree.cbh + z*tree.ch #+ z_m

                            if self.tree_id.insert([c_x, c_y, c_z], int(tree.id)):
                                n -= 1
                                attemps = 0
                            else:
                                attemps += 1

                            if attemps > 20:
                                self.tree_id.rewind_insert_track()
                                break
                        placed = True

                    track = self.tree_id.track
                    self.tree_sp.trace_track(track, tree.sp)
                    self.bulk_density.trace_track(track, tree.bulk_density)
                    self.sav.trace_track(track, 1/2000.0*0.348)
                    self.moisture.trace_track(track, 0.68)

    def distribute_surface(self):

        fact = 30 // self.res[0]

        bd_array = np.zeros((self.h*fact, self.w*fact))
        depth_array = np.zeros_like(bd_array)
        sav_array = np.zeros_like(bd_array)
        moisture_array = np.zeros_like(bd_array)

        for i in range(self.h*fact):
            print(i)
            for j in range(self.w*fact):
                theta = np.random.random()*2*np.pi
                radius = np.random.random()*6
                x = np.cos(theta)*radius
                y = np.sin(theta)*radius
                int_x = j//fact + int(x + 0.5)
                int_y = i//fact + int(y + 0.5)
                while (int_x < 0) or (int_x > self.w-1) or (int_y < 0) or (int_y > self.h-1):
                    theta = np.random.random()*2*np.pi
                    radius = np.random.random()*6
                    x = np.cos(theta)*radius
                    y = np.sin(theta)*radius
                    int_x = j//fact + int(x + 0.5)
                    int_y = i//fact + int(y + 0.5)
                bd_array[i,j] = SB40_PARAMS[self.surface_array[int_y, int_x]][1]
                depth_array[i,j] = SB40_PARAMS[self.surface_array[int_y, int_x]][4]
                sav_array[i,j] = SB40_PARAMS[self.surface_array[int_y, int_x]][2]
                moisture_array[i,j] = SB40_PARAMS[self.surface_array[int_y, int_x]][3]*0.25/100

        self.bulk_density.domain[:,:,0] = ndimage.gaussian_filter(bd_array, sigma=5)
        self.fuel_depth.domain[:,:,0] = ndimage.gaussian_filter(depth_array, sigma=5)
        self.sav.domain[:,:,0] = ndimage.gaussian_filter(sav_array, sigma=5)
        self.moisture.domain[:,:,0] = ndimage.gaussian_filter(moisture_array, sigma=5)

import matplotlib.pyplot as plt
import matplotlib

if __name__ == '__main__':

    """
    ele = elevation.DEM(2200, 2200, 30, 30)
    dem = ele.interpolate()
    print(dem)
    plt.imshow(dem)
    plt.show()
    x = np.arange(0, dem.shape[1], 1.0)
    y = np.arange(0, dem.shape[0], 1.0)
    x,y = np.meshgrid(x,y)
    dem = dem - np.min(dem)
    terrain = pv.StructuredGrid(y,x,dem)
    """

    forest = Forest(38.98531294306404, -120.64584148255688, 1000, 1000)

    print('sampling FIA plots')
    forest.sample_plots()

    print('placing trees')
    forest.distribute_canopy()

    print('surface fuels')
    forest.distribute_surface()

    print('plotting')
    pv.set_plot_theme('night')
    p = pv.Plotter()
    viewer = domain.View(forest.moisture)
    grid = viewer.show3d()
    grid.plot()
    #p.add_mesh(grid)
    #p.add_mesh(terrain, color='#bfa300')
    #p.show(show_axes=True)

    from scipy.io import FortranFile
    def save(data, filename):
        data = data.astype(np.float32)
        data = data.T
        f = FortranFile(filename, 'w', 'uint32')
        f.write_record(data)

    save(forest.bulk_density.domain, '../../output/bulk_density.dat')
    save(forest.moisture.domain, '../../output/moisture.dat')
    save(forest.sav.domain, '../../output/sav.dat')
    save(forest.fuel_depth.domain, '../../output/depth.dat')
    save(forest.dem_array, '../../output/topo.dat')


"""
Generator for canopy fuels
"""

# built-in{fil
import time
import json
import sys
import math

# external
import cupy as cp
import numpy as np
import gdal
sys.path.append('../../../qFIA/src/')
import qfia

# package-level
import domain

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

SPECIES_GROUP_COLORS = {

}

# load fia plot cn indices
with open('../../data/fia_ids.txt', 'rb') as f:
    FIA_IDX = json.load(f)

# load treelist layer
FIA_TREELIST = gdal.Open('../../data/treelist.tif')
SB40_FUELS = gdal.Open('../../data/sb40.tif')


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

        self.bulk_density = self.weight.foliage/(self.volume*self.cd)

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

        self.treelist_nodata = FIA_TREELIST.GetRasterBand(1).GetNoDataValue()
        self.treelist_array = FIA_TREELIST.ReadAsArray(x, y, w, h)

        self.w = w
        self.h = h

    def sample_plots(self):

        sampled_cns = []

        samples = {}

        xy = []

        for i in range(self.h):
            for j in range(self.w):
                fia_idx = self.treelist_array[i,j]
                if fia_idx == self.treelist_nodata:
                    continue
                plt_cn = FIA_IDX[str(fia_idx)]
                #if plt_cn in sampled_cns:
                    #continue
                #sampled_cns.append(plt_cn)
                samples[plt_cn] = []
                plot = qfia.query(plt_cn)

                for subplot in plot:
                    for tree in subplot:

                        if (tree.dia) and (tree.statuscd == 1):
                            if tree.dia >= 5.0:
                                n = self.n_subplot_trees()
                            else:
                                n = self.n_microplot_trees()

                            for k in range(n):
                                spcd = str(int(tree.spcd))
                                dia = tree.dia*2.54
                                ht = tree.actualht*0.3048
                                cr = tree.cr
                                theta = np.random.random()*2*np.pi
                                d = np.random.random()*30
                                x = np.cos(theta)*d+15
                                y = np.sin(theta)*d+15
                                x += j*30
                                y += i*30
                                xy.append([x,y,dia/100])
                                samples[plt_cn].append(Tree(spcd, dia, ht, cr))

        return xy

    def n_subplot_trees(self):

        r = np.random.random()
        if r < 0.338:
            n = 2
        else:
            n = 1

        return n

    def n_microplot_trees(self):

        r = np.random.random()
        if r < 0.672:
            n = 17
        else:
            n = 16

        return n


import matplotlib.pyplot as plt
import matplotlib
if __name__ == '__main__':

    for i in range(100):
        tree1 = Tree('202', 50.0, 30, 0.5)
        print(tree1.volume)
        print(tree1.id)

    forest = Forest(500,1200,30,30)
    xy = forest.sample_plots()

    patches = [plt.Circle((i[0],i[1]),i[2]) for i in xy]

    fig, ax = plt.subplots()

    coll = matplotlib.collections.PatchCollection(patches, facecolors='black')
    #ax.imshow(forest.treelist_array, alpha=0.2)
    ax.add_collection(coll)
    ax.set_aspect('equal')
    ax.set_ylim(0,900)
    ax.set_xlim(0,900)
    plt.show()

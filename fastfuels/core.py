"""
The fastfuels module loads a .fio directory-in-file resource on which a user can
perform spatial queries, view fuels data in 3D and export to QuicFire.
"""

__author__     = "Holtz Forestry LLC"
__date__       = "16 November 2020"
__version__    = "0.4.1"
__maintainer__ = "Lucas Wells"
__email__      = "lucas@holtzforestry.com"
__status__     = "Prototype"

# external imports
import colorcet # pip3 install colorcet
import numpy as np # pip3 install numpy
import pyvista as pv # pip3 install pyvista
from scipy.io import FortranFile #pip3 install scipy
import zarr # pip3 install zarr

# new imports in version 0.4.1
import geopandas as gpd
from rasterio.features import geometry_mask
from shapely.geometry import Polygon

# --------------
# USERSPACE DEFS
# --------------
def open(fname):
    """
    Helper function for opening a .fio file. Additional user actions on
    fuels data can be accessed through the class instance returned by this
    function.

    Args:
        fname (str): path and filename of the .fio resource

    Returns:
        An instance of FuelsIO
    """

    return FuelsIO(fname)


# --------------
# HELPER CLASSES
# --------------
class AlbersEqualAreaConic:
    """
    Implements forward and inverse projection on Albers Equal Area Conic

    Attributes:
        lambda_0 (float): origin of longitude
        a (float): semi-major axis of earth (GRS80)
        e2 (float): eccentricity squared
        e (float): eccentricity
        n (float): stored map constant
        C (float): stored map constant
        rho_0 (float): stored map constant
    """

    def __init__(self, phi_1=29.5, phi_2=45.5, phi_0=23.0, lambda_0=-96.0):
        """
        Creates an instance of AlbersEqualAreaConic, initializes attributes and
        computes map constants

        Note:
            Geographic constants are based on the GRS 1980 ellipsoid

        Args:
            phi_1 (float): first standard parallel
            phi_2 (float): second standard parallel
            phi_0 (float): origin of latitude
            lambda_0 (float): origin of longitude
        """

        # convert map params to radians
        phi_0 = np.radians(phi_0)
        phi_1 = np.radians(phi_1)
        phi_2 = np.radians(phi_2)
        self.lambda_0 = np.radians(lambda_0)

        # GRS 1980 REFERENCE ELLIPSIOD CONSTANTS
        # geographic constants
        self.a = 6378137.0
        # derived geometrical constants
        f = 1.0/298.2572221010042 # flattening
        self.e2 = 2*f - f**2 # eccentricity squared
        self.e = np.sqrt(self.e2) # eccentricity

        # preliminaries
        m_1 = self._m(phi_1)
        m_2 = self._m(phi_2)
        q_0 = self._q(phi_0)
        q_1 = self._q(phi_1)
        q_2 = self._q(phi_2)

        # map constants
        self.n = (m_1**2 - m_2**2)/(q_2 - q_1)
        self.C = m_1**2 + self.n*q_1
        self.rho_0 = self.a*(self.C - self.n*q_0)**0.5/self.n

    def _m(self, phi):
        """Private member
        Convenience method for computing map constants
        """

        return np.cos(phi)/np.sqrt(1 - self.e2*(np.sin(phi))**2)

    def _q(self, phi):
        """Private member
        Another convenience method for computing map constants
        """

        return (1 - self.e2)*(np.sin(phi)/(1 - self.e2*(
            np.sin(phi))**2) - (1.0/(2*self.e))*np.log((1-self.e*np.sin(
            phi))/(1 + self.e*np.sin(phi))))

    def forward(self, lat, lon):
        """
        Performs forward projection from geodetic coordinates to projected
        coordinates

        Args:
            lat (float): latitude
            lon (float): longitude

        Returns:
            (x,y) coordinate projected in Albers Equal Area Conic
        """

        # convert to radians for numpy trig functions
        lat = np.radians(lat)
        lon = np.radians(lon)

        # preliminaries
        q = self._q(lat)
        rho = self.a*(self.C - self.n*q)**0.5/self.n
        theta = self.n*(lon - self.lambda_0)

        # retrieve the projected coordinates
        x = rho*np.sin(theta)
        y = self.rho_0 - rho*np.cos(theta)

        return x,y

    def inverse(self, x, y):
        """
        Performs inverse projection from Albers to geodetic coordinates

        Args:
            x (float): x projected in Albers
            y (float): y projected in Albers

        Returns:
            lat and lon in geodetic coordinates
        """

        # preliminaries
        p = np.sqrt(x*x + (self.rho_0 - y)**2)
        theta = np.arctan2(x, self.rho_0 - y)
        q = (self.C - ((p*p)*(self.n**2))/(self.a**2))/self.n

        # convergence criteria
        epsilon = 1e-6

        # iterate latitude calculation until convergence
        phi = np.sin(q/2)
        next_phi = self._inverse_iteration(phi, q)
        while (np.abs(np.degrees(phi) - np.degrees(next_phi)) > epsilon):
            phi = next_phi
            next_phi = self._inverse_iteration(phi, q)

        return np.degrees(phi), np.degrees(self.lambda_0 + theta/self.n)

    def _inverse_iteration(self, phi, q):
        """Private member
        Formula to iterator until convergence for inverse projection of latitude
        """

        return np.radians(np.degrees(phi) + np.degrees((1 -
            self.e2*(np.sin(phi)**2)**2)/(2*np.cos(phi))*(
            (q/(1-self.e2)) - (np.sin(phi)/(1-self.e2*(np.sin(phi)**2))) +
            (1/(2*self.e))*np.log((1 - self.e*np.sin(phi))/(1 +
            self.e*np.sin(phi))))))


class Viewer:
    """
    This class wraps pyvista which wraps VTK. Provides simple convenience
    functions for plotting UniformGrids in 3D.

    Attributes:
        data (dict): dictionary of 3D fuel parameter arrays
        plotter: Pyvista plotter instance

    """

    def __init__(self, data):
        """
        Configures the Pyvista plotter instance

        Args:
            data (dict): dictionary of 3D fuel parameter arrays
        """

        # set pv theme
        pv.set_plot_theme('document')
        self.data = data
        self.plotter = pv.Plotter(title='FastFuels')

    def add(self, property):
        """
        Adds a 3D array to the plotter instance

        Args:
            property (str): fuel parameter to show, parameter string must be in
                the key list of the data dictionary

        Note:
            Use the `get_properties()` method to get available fuel properties.
        """

        # extract the property array from the data dictionary
        fp = self.data[property]

        # move zeros to -1 for thresholding
        fp[fp == 0] = -1

        # convert the 3D array to a Pyvista UniformGrid
        grid = pv.UniformGrid()
        grid.dimensions = np.array(fp.shape) + 1
        grid.spacing = [1,1,1]
        grid.cell_arrays['values'] = fp.flatten(order='F')
        grid = grid.threshold(0)

        # conditionally select the colormap based on property
        if property == 'species_group':
            cm = colorcet.glasbey[:10]
        else:
            cm = 'rainbow'

        # add the array to the plotter instance
        self.plotter.add_mesh(grid, cmap=cm)

    def show(self):
        """
        This will halt the program and display an interactive window populated
        will a 3D array
        """

        self.plotter.show()


class FuelsIO:
    """
    This class handles the .fio resource, extracts metadata and performs spatial
    queries.

    Attributes:
        fio_file: Zarr file manager object
        albers (AlbersEqualAreaConic object): Instance of the
            `AlbersEqualAreaConic` class
        cache_limit (float): Memory cache limit in bytes for spatial queries
        fmwriter (FuelModleWriter object): Instance of the `FuelModleWriter`
            class
    """

    def __init__(self, fio_fname, roads_fname, water_fname, fire_fname):
        """
        Opens a connection to the fio resource, extracts metadata and
        initializes helper classes

        Args:
            fname (str): path and filename of the fio resource.
        """

        # No longer using the remote demo fio file hosted on GCP, moved data
        # over to UCSD server

        # open connection to fio resource
        #if fname == 'remote':
        #    print('connecting to remote FIO server...')
        #    gcs = gcsfs.GCSFileSystem()
        #    self.fio_file = zarr.open(gcs.get_mapper('gs://ca-11-2020/demo.fio'), 'r')
        #else:

        self.fio_file = zarr.open(fname, 'r')

        # get metadata and datasets
        self.extract_meta_data()
        self.parse_contents()

        # instantiate helper classes
        self.albers = AlbersEqualAreaConic()
        self.fmwriter = FireModelWriter()

        # default cache limit
        self.cache_limit = 1e9

        # class attributes for the new stuff from Daniel R-G
        self.roads_fname = roads_fname
        self.water_fname = water_fname
        self.fire_fname = fire_fname

    def extract_meta_data(self):
        """
        Gets metadata from fio resource attributes
        """

        self.extent_x1, self.extent_y1, self.extent_x2, self.extent_y2 = self.fio_file.attrs['extent']

        # new fio version changed extent format key from "extent_format" to
        # "extent_fmt"
        self.extent_fmt = self.fio_file.attrs['extent_fmt']

        self.n_cols = self.extent_x2 - self.extent_x1
        self.n_rows = self.extent_y1 - self.extent_y2

        self.proj = self.fio_file.attrs['proj']
        self.res = self.fio_file.attrs['resolution']
        self.units = self.fio_file.attrs['units']

        # old attributes (these are removed in new version of fio files)
        #self.dim = self.fio_file.attrs['dimensions']
        #self.dim_fmt = self.fio_file.attrs['dim_format']

    def parse_contents(self):
        """
        Get references to fuel array datasets

        Note:
            This method will not load the arrays to memory
        """

        # extract surface fuel properties and elevation
        surface_group = self.fio_file['surface']
        # change key from "loading" to "load" in new fio
        self.surface_loading = surface_group['load']
        # change from "fuel_depth" to "depth" in new fio
        self.surface_fuel_depth = surface_group['depth']
        self.surface_sav = surface_group['sav']
        # change key from "elevation" to "dem" in new fio
        self.elevation = surface_group['dem']

        # removed emc parameter in new version of fio
        #self.surface_emc = surface_group['emc']

        # extract canopy fuel properties
        canopy_group = self.fio_file['canopy']
        # changed key from "bulk_density" to "bd" in new fio
        self.canopy_bulk_density = canopy_group['bd']
        self.canopy_sav = canopy_group['sav']

        # removed species_group parameter in new fio
        #self.canopy_sp_group = canopy_group['species_group']

    def get_extent(self, mode='projected'):
        """
        Get the geographic extent in projected or geographic coordinates

        Args:
            mode (str, default='projected'): mode of extent coordinate system;
                projected or geographic
        """

        if mode == 'projected':
            return self.extent_x1, self.extent_y1, self.extent_x2, self.extent_y2
        elif mode == 'geographic':
            lat1, lon1 = self.albers.inverse(self.extent_x1, self.extent_y1)
            lat2, lon2 = self.albers.inverse(self.extent_x2, self.extent_y2)
            return lon1, lat1, lon2, lat2

    def query(self, lon, lat, radius):
        """
        Performs a geographic spatial query on the fio resource. Extent of the
        query is defined by the point (lat, lon) and a square bounding an inscribed
        circle defined by the radius parameter

        Args:
            lon (float): longitude
            lat (float): latitude
            radius (int): radius of extent circle in meters
        """

        return self.query_geographic(lon, lat, radius)

        # changing the way we query fuels in versino 0.3.1
        """
        if mode == 'relative':
            return self.query_relative(a, b)
        elif mode == 'projected':
            return self.query_projected(a, b)
        elif mode == 'geographic':
            return self.query_geographic(a, b)
        """

    def check_bounds(self, x1, y1, x2, y2):
        """
        Returns true if in bounds and false if out of bounds
        """

        if x1 < 0:
            return 0
        elif y1 < 0:
            return 0
        elif x2 > self.n_cols:
            return 0
        elif y2 > self.n_rows:
            return 0
        else:
            return 1

    def check_cache(self, x1, y1, x2, y2):
        """
        Returns true if query size is less than cache limit and false if
        query size of greater than cache limit
        """

        if 1700*((x2 - x1)*(y2 - y1)) > self.cache_limit:
            return 0
        return 1

    def query_relative(self, a, b):

        x1, y1 = a
        x2, y2 = b

        if self.check_bounds(x1, y1, x2, y2):
            if self.check_cache(x1, y1, x2, y2):
                return self.slice_and_merge(y1, y2, x1, x2)
            else:
                print('ERROR: area too large')
        else:
            print('ERROR: bounding box query out of bounds')
            return -1

    def query_projected(self, a, b):

        x1, y1 = a
        x2, y2 = b

        x1_rel = int(x1 - self.extent_x1)
        y1_rel = int(self.extent_y1 - y1)
        x2_rel = int(x2 - self.extent_x1)
        y2_rel = int(self.extent_y1 - y2)

        return self.query_relative((x1_rel, y1_rel), (x2_rel, y2_rel))

    def query_geographic(self, lon, lat, radius):

        x1, y1 = self.albers.forward(lat, lon)
        x1 -= radius
        y1 += radius

        x2 = x1 + radius*2
        y2 = y1 - radius*2

        self.bbox = (x1, y1, x2, y2)

        return self.query_projected((x1, y1), (x2, y2))

    def slice_and_merge(self, y1, y2, x1, x2):
        """
        Queries the fuel array datasets and loads to memory

        Returns:
            An instance of the FuelsROI class
        """

        data_dict = {}

        canopy_bulk_density = self.canopy_bulk_density[y1:y2, x1:x2, :].astype(np.float32)
        canopy_bulk_density = (canopy_bulk_density/255)*2.0

        canopy_moisture = np.zeros_like(canopy_bulk_density)
        canopy_moisture[canopy_bulk_density != 0] = 1.0 # hardcoded for now

        surface_loading = self.surface_loading[y1:y2, x1:x2].astype(np.float32)
        surface_loading = (surface_loading/255)*3.0

        canopy_bulk_density[:,:,0] = surface_loading
        data_dict['bulk_density'] = canopy_bulk_density

        canopy_sav = self.canopy_sav[y1:y2, x1:x2, :].astype(np.float32)
        canopy_sav = (canopy_sav/255)*8000.0

        surface_sav = self.surface_sav[y1:y2, x1:x2].astype(np.float32)
        surface_sav = (surface_sav/255)*8000.0

        canopy_sav[:,:,0] = surface_sav
        data_dict['sav'] = canopy_sav

        # surface_emc removed in new fio
        #surface_moisture = self.surface_emc[y1:y2, x1:x2]
        canopy_moisture[:,:,0] = 0.2 # hardcoded for new
        data_dict['moisture'] = canopy_moisture

        fuel_depth = np.zeros_like(canopy_sav)
        fuel_depth[:,:,0] = self.surface_fuel_depth[y1:y2, x1:x2].astype(np.float32)
        fuel_depth = (fuel_depth/255)*2.0
        data_dict['fuel_depth'] = fuel_depth

        # removed species group parameter in new fio
        #data_dict['species_group'] = self.canopy_sp_group[y1:y2, x1:x2, :]

        data_dict['elevation'] = self.elevation[y1:y2, x1:x2]

        data_dict = self.mask_features(data_dict)

        return FuelsROI(data_dict)

    # -------------------------------------------------------------
    # New methods from Daniel R-G (Not tested, last minute request)
    # -------------------------------------------------------------
    def roads_buffer(self):
        '''
        Uses MAF/TIGER Feature Class Code Definitions to assign particular widths to
        roads and returns a polygon shapefile containing said buffers
        '''

        sd_roads  = gpd.read_file(self.roads_fname, self.bbox)
        if (sd_roads.shape[0]>0):

            # S1100: Primary Road, S1200: Secondary Road, S1400: Local Road, S1500: Vehicular Trail
            # S1630: Ramp, S1640: Service Drive, S1710: Walkway, S1720: Stairway, S1730: Alley,
            # S1740: Private Road, S1750: Internal USCB Use, S1780: Parking Lot Road, S1820: Bike Path
            # S1830: Bridle Path
            road_categories = ['S1100','S1200','S1400','S1500','S1630','S1640',
                               'S1710','S1720','S1730','S1740','S1750','S1780',
                               'S1820','S1830']
            road_buffers = [15., 10., 5., 3., 8., 5., 0., 0., 2., 3., 3., 2., 0., 0. ]
            roads_geo = gpd.GeoSeries()
            roads_geo.set_crs("ESRI:102039")
            polygon = Polygon([(self.bbox[0], self.bbox[1]),
                               (self.bbox[0], self.bbox[3]),
                               (self.bbox[2], self.bbox[3]),
                               (self.bbox[2], self.bbox[1]),
                               (self.bbox[0], self.bbox[1])])
            poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon], crs=sd_roads.crs)
            rds = gpd.clip(sd_roads,poly_gdf)
            #rds.plot(color='k',ax=ax)

            for i, elem in enumerate(road_categories):
                if (road_buffers[i] > 0):
                    roads_geo = roads_geo.append(rds.query("MTFCC == '{}'".format(elem)).buffer(road_buffers[i]))

            self.remove_shapefile(rds)

        else:
            pass

    def water_buffer(self):
        '''
        Uses NHD Hydrography Features to assign particular widths to the vector shapefile
        containing stream information and returns a polygon shapefile with the stream buffers
        '''

        sd_water = gpd.read_file(self.water_fname, self.bbox)
        if (sd_water.shape[0]>0):

            # 336: Canal/Ditch, 460: Stream/River, 558: Artificial Path,
            # 334: Connector, 428: Pipeline, 420: Underground Conduit
            water_categories = [460, 558,336, 334, 428, 420]
            water_buffers = [3.,6.,2.,0.,0.,0.]
            water_geo = gpd.GeoSeries()
            water_geo.set_crs("ESRI:102039")
            polygon = Polygon([(self.bbox[0], self.bbox[1]),
                               (self.bbox[0], self.bbox[3]),
                               (self.bbox[2], self.bbox[3]),
                               (self.bbox[2], self.bbox[1]),
                               (self.bbox[0], self.bbox[1])])
            poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon], crs=sd_water.crs)
            wtr = gpd.clip(sd_water,poly_gdf)

            for i, elem in enumerate(water_categories):
                if (water_buffers[i] > 0):
                    water_geo = water_geo.append(wtr.query("FType == {}".format(elem)).buffer(water_buffers[i]))

            self.remove_shapefile(water_geo)

        else:
            pass

    def fire_buffer(self):
        '''
        Uses year since last burn to return a polygon shapefile with just the burns
        over the required years
        '''
        sd_fire = gpd.read_file(self.fire_fname, self.bbox)
        if (sd_fire.shape[0]>0):
            fire_geo = gpd.GeoSeries()
            fire_geo.set_crs("ESRI:102039")
            fire_geo = fire_geo.append(sd_fire.query("year > 2015").buffer(0))

            polygon = Polygon([(self.bbox[0], self.bbox[1]),
                               (self.bbox[0], self.bbox[3]),
                               (self.bbox[2], self.bbox[3]),
                               (self.bbox[2], self.bbox[1]),
                               (self.bbox[0], self.bbox[1])])
            poly_gdf = gpd.GeoDataFrame([1], geometry=[polygon], crs=sd_fire.crs)

            fire_clip = gpd.clip(fire_geo, poly_gdf)
            self.remove_shapefile(fire_clip)

        else:
            pass

    def remove_shapefile(self, shapefile):
        '''
        Uses shapefile to mask out the fuels in the array provided
        '''

        total_bounds = shapefile.total_bounds
        geoT = (1.0,0.0,total_bounds[0],
               0.0,1.0,total_bounds[1])
        shp_mask = geometry_mask(shapefile.geometry, np.shape(self.feature_mask), geoT)
        self.feature_mask *= shp_mask

    def apply_mask(self, arr):
        '''
        Applies mask to specified array
        '''

        if (np.shape(self.feature_mask) == np.shape(arr[:,:,0])):
            masked_array = np.copy(arr)
            (ny,nx,nz) = np.shape(arr)
            for l in range(nz):
                masked_array[:,:,l] *= self.feature_mask
            return masked_array
        else:
            print('Mismatch in mask and array size')

    def mask_features(self, data_dict):

        # init mask
        (xmin, ymax, xmax, ymin) = self.bbox
        x_ext, y_ext = [int(xmax-xmin),int(ymax-ymin)]
        self.feature_mask = np.ones([y_ext,x_ext],dtype=bool)

        # construct mask
        self.fire_buffer()
        self.water_buffer()
        self.roads_buffer()
        self.feature_mask = np.flipud(self.feature_mask)

        for array_key in data_dict.keys():
            data_dict[array_key] = self.apply_mask(self.feature_mask, data_dict[array_key])

        return data_dict


class FuelsROI:
    """
    Handles viewing and writing of 3D fuel arrays following a FuelsIO
    spatial query.

    Attributes:
        data_dict (dictionary): Keys are fuel properties and values are 3D
            arrays
    """

    def __init__(self, data_dict):
        """
        initializes attributes and instantiates Viewer and FuelModelWriter
        objects.
        """

        self.data_dict = data_dict
        self.viewer = Viewer(data_dict)
        self.writer = FireModelWriter()

    def get_properties(self):
        """
        Get the available fuel properties

        Returns:
            list of fuel properties
        """

        return list(self.data_dict.keys())

    def view(self, property):
        """
        Display 3D array of given fuel property

        Args:
            property (str): fuel property
        """

        self.viewer.add(property)
        self.viewer.show()

    def write(self, path, model='quicfire', res_xyz=[1,1,1]):
        """
        Writes fuel arrays to a fire model. Currently only implements QUICFire
        """

        if model == 'quicfire':
            print('writing QUICFire input files...')
            self.writer.write_to_quicfire(self.data_dict['bulk_density'],
                path + '/' + 'rhof.dat', res_xyz)
            self.writer.write_to_quicfire(self.data_dict['sav'],
                path + '/' + 'sav.dat', res_xyz)
            self.writer.write_to_quicfire(self.data_dict['moisture'],
                path + '/' + 'moisture.dat', res_xyz)
            self.writer.write_to_quicfire(self.data_dict['fuel_depth'],
                path + '/' + 'fueldepth.dat', res_xyz)
            self.writer.write_to_quicfire(self.data_dict['elevation'],
                path + '/' + 'elevation.dat', res_xyz)
            print('complete')
        elif model == 'wfds':
            print('wfds writer not implemented')


class FireModelWriter:
    """
    Writes fuel arrays to disk as fire model input files
    """

    def write_to_quicfire(self, data, fname, res_xyz):
        """
        Write fuel array as QF input file
        """

        print(f'Writing data to {fname}...')
        rx, ry, rz = res_xyz

        if res_xyz == [1,1,1]:
            pass
        elif res_xyz == [2,2,1]:

            if len(data.shape) == 3:
                h,w,d = data.shape
                # average pool subsampling
                data = data.reshape((h//2, 2, w//2, 2, d)).mean(3).mean(1)
            elif len(data.shape) == 2:
                h,w = data.shape
                data = data.reshape((h//2, 2, w//2, 2)).mean(3).mean(1)
        else:
            print(f'Resolution can be either [1,1,1] or [2,2,1], not {res_xyz}. ' +
                'Defaulting to [1,1,1] resolution')

        print(f'Output resolution is x: {rx}, y: {ry}, z: {rz}\n')

        data = data.astype(np.float32)
        data = data.T
        f = FortranFile(fname, 'w', 'uint32')
        f.write_record(data)

    def write_to_wfds(self, data, fname):
        """
        TODO: write to FDS input file
        """
        pass

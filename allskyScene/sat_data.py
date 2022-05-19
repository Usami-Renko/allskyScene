'''
Description: class of satellite data for all-sky assimilation overview
Author: Hejun Xie
Date: 2022-05-19 23:36:34
LastEditors: Hejun Xie
LastEditTime: 2022-05-20 01:08:59
'''

# Global imports
import numpy as np
# satpy (for FY4A-MWRI)
from satpy.scene import Scene
# h5py (for FY3D)
import h5py

# Local imports
# cython speedup (for FY3D)
from .monotonic import monotonicalize 

class AGRI_data(object):
    def __init__(self, filenames):
        self.filenames = filenames
        self.scn = None
        self.tb_ir = None
        self.area_def = None
        print(self.filenames)

    def read_agri_ndisk(self, ir_channel):
        self.scn = Scene(self.filenames, reader='agri_l1')
        self.scn.load([ir_channel], generate=False, calibration='brightness_temperature')
        self.tb_ir       = self.scn[ir_channel].values - 273.15 # [K] --> [Celsius Degree]
        self.area_def    = self.scn[ir_channel].attrs['area']


class MWRI_data(object):
    def __init__(self, filenames):
        self.filenames = filenames
        self.BTs = list()
        self.lats = list()
        self.lons = list()
        self.divide_meridians = list()
        self.valid_meridian = [-720.,-360.,0.,360.,720.]
    
    def read_mwri_swaths(self, mwri_channel):
        for filename in self.filenames:
            print('read {}'.format(filename))
            result = self.read_mwri_swath(filename, mwri_channel)
            self.BTs.append(result[0])
            self.lats.append(result[1])
            self.lons.append(result[2])
        
        self.calc_divide_meridian()

    def read_mwri_swath(self, filename, mwri_channel):
        file_handle = h5py.File(filename, 'r')
        lat_handler = file_handle['Geolocation']['Latitude']
        lon_handler = file_handle['Geolocation']['Longitude']
        InvalidFlagGeoLoc = (lat_handler[:] < lat_handler.attrs['valid_range'][0]) |\
                            (lat_handler[:] > lat_handler.attrs['valid_range'][1]) |\
                            (lat_handler[:] == lat_handler.attrs['FillValue']) |\
                            (lon_handler[:] < lon_handler.attrs['valid_range'][0]) |\
                            (lon_handler[:] > lon_handler.attrs['valid_range'][1]) |\
                            (lon_handler[:] == lon_handler.attrs['FillValue'])
        lat = lat_handler[:]
        lon = lon_handler[:]
        # print('Invalid GeoLoc number = {} '.format(np.sum(InvalidFlagGeoLoc[:])))

        BT_handler = file_handle['Calibration']['EARTH_OBSERVE_BT_10_to_89GHz']
        InvalidFlag =   (BT_handler[:] < BT_handler.attrs['valid_range'][0]) |\
                        (BT_handler[:] > BT_handler.attrs['valid_range'][1]) |\
                        (BT_handler[:] == BT_handler.attrs['FillValue'])
        BT = BT_handler[:] * BT_handler.attrs['Slope'] + BT_handler.attrs['Intercept']
        (nchannels, nlines, nscans) = BT.shape

        sea = (file_handle['Calibration']['LandSeaMask'][:] == 3)

        BT = np.where(InvalidFlag, np.nan, BT)
        BT = np.where(InvalidFlagGeoLoc, np.nan, BT)
        BT = np.where(~sea, np.nan, BT)
        # BT = np.where((lat>50.0) | (lat<-50.0), np.nan, BT)
        
        # swap byte from big endian to small endian for cython
        lon = lon.newbyteorder().byteswap(inplace=True)
        lat = lat.newbyteorder().byteswap(inplace=True)
        # monotonicalize the longitude
        monotonicalize(lon, lat, -180., 180., 65535.)
        
        return BT[mwri_channel,...], lat, lon
    
    def calc_divide_meridian(self):
        # Get divide meridians 
        for iswath, lon in enumerate(self.lons):
            # 1. case of longitude never cross over [-720.,-260.,0.,360.,720.]
            self.divide_meridians.append(0.0)
            for meridian in self.valid_meridian:
                # 2. case of longitude cross over [-720.,-260.,0.,360.,720.]
                if np.nanmin(lon) < meridian <= np.nanmax(lon):
                    self.divide_meridians[iswath] = meridian
                    print('divide_meridian= {}'.format(meridian)); break

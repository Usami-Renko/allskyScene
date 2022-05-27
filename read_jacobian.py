'''
Description: read jacobian data
Author: Hejun Xie
Date: 2022-05-26 23:29:03
LastEditors: Hejun Xie
LastEditTime: 2022-05-27 15:54:34
'''

import numpy as np

def read1Line(fin):
    while True:
        line = fin.readline()
        if line not in ['\n', '']:
            return line.strip('\n')

class rttov_profile(object):
    def __init__(self):
        self.metaInfo = dict()
        self.viewingConditions = dict()
        self.blackBodyCloud = dict()
        self.skin = dict()
        self.s2m = dict()
        self.gas_units = None
        self.p = None
        self.T = None
        self.q = None
        self.cld_profile = None

    def readfromfile(self, filehandle):
        headLine = read1Line(filehandle) # start with an empty line
        idLine = read1Line(filehandle)
        self.metaInfo['date'] = read1Line(filehandle).split()[-1]
        self.metaInfo['time'] = read1Line(filehandle).split()[-1]
        self.metaInfo['nlevels'] = int(read1Line(filehandle).split()[-1])
        self.metaInfo['nlayers'] = int(read1Line(filehandle).split()[-1])
        
        ViewingConditionLine = read1Line(filehandle)
        self.viewingConditions['zenangle'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['azangle'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['sunzenangle'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['sunazangle'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['elevation'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['latitude'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['longitude'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['EMF'] = float(read1Line(filehandle).split()[-1])
        self.viewingConditions['cosEMF'] = float(read1Line(filehandle).split()[-1])

        BlackBodyCloudLine = read1Line(filehandle)
        self.blackBodyCloud['CTP'] = float(read1Line(filehandle).split()[-1]) # [hPa]
        self.blackBodyCloud['CF'] = float(read1Line(filehandle).split()[-1]) # [0-1]

        SkinLine = read1Line(filehandle)
        self.skin['surf-type'] = int(read1Line(filehandle).split()[-1]) # 0=land, 1=sea, 2=sea-ice
        self.skin['water-type'] = int(read1Line(filehandle).split()[-1]) # 0=fresh water, 1=ocean water
        self.skin['skinT'] = float(read1Line(filehandle).split()[-1]) # [K]
        self.skin['snowFraction'] = float(read1Line(filehandle).split()[-1]) # [0-1]
        self.skin['soilMoisture'] = float(read1Line(filehandle).split()[-1]) # [m3/m3]
        self.skin['salinity'] = float(read1Line(filehandle).split()[-1]) # [%o]
        self.skin['foamFraction'] = float(read1Line(filehandle).split()[-1]) # [0-1]
        FASTEMLine = read1Line(filehandle)
        self.skin['FATSTEMParams'] = np.array([float(param) for param in read1Line(filehandle).split()])

        s2mLine = read1Line(filehandle)
        self.s2m['p2m'] = float(read1Line(filehandle).split()[-1]) # [hPa]
        self.s2m['t2m'] = float(read1Line(filehandle).split()[-1]) # [K]
        self.s2m['q2m'] = float(read1Line(filehandle).split()[-1]) # [ppmv] / [kg/kg]
        self.s2m['o2m'] = float(read1Line(filehandle).split()[-1]) # [ppmv] / [kg/kg]
        self.s2m['u10m'] = float(read1Line(filehandle).split()[-1]) # [m/s]
        self.s2m['v10m'] = float(read1Line(filehandle).split()[-1]) # [m/s]
        self.s2m['windFetch'] = float(read1Line(filehandle).split()[-1]) # [m]

        GasUnitLine = read1Line(filehandle)
        ProfileHeadLine = read1Line(filehandle)
        self.p = np.zeros((self.metaInfo['nlevels']), dtype='float32')
        self.T = np.zeros((self.metaInfo['nlevels']), dtype='float32')
        self.q = np.zeros((self.metaInfo['nlevels']), dtype='float32')

        for ilevel in range(self.metaInfo['nlevels']):
            segs = read1Line(filehandle).split()
            self.p[ilevel] = float(segs[1]) # [hPa]
            self.T[ilevel] = float(segs[2]) # [K]
            self.q[ilevel] = float(segs[3]) # [ppmv] / [kg/kg]
        
class rttov_cld_profile(object):
    def __init__(self):
        self.nlevels = None
        self.userCldFra = None
        self.phTop = None
        self.phBot = None
        self.cc = None
        self.clw = None
        self.rain = None
        self.ciw = None
        self.snow = None

    def readfromfile(self, filehandle):
        headLine = read1Line(filehandle)
        self.nlevels = int(read1Line(filehandle).split()[-1])
        self.userCldFra = float(read1Line(filehandle).split()[-1]) # [0-1]

        ProfileHeadLine = read1Line(filehandle)
        self.phTop = np.zeros((self.nlevels), dtype='float32')
        self.phBot = np.zeros((self.nlevels), dtype='float32')
        self.cc = np.zeros((self.nlevels), dtype='float32')
        self.clw = np.zeros((self.nlevels), dtype='float32')
        self.rain = np.zeros((self.nlevels), dtype='float32')
        self.ciw = np.zeros((self.nlevels), dtype='float32')
        self.snow = np.zeros((self.nlevels), dtype='float32')
        
        for ilevel in range(self.nlevels):
            segs = read1Line(filehandle).split()
            self.phTop[ilevel] = float(segs[1]) # [hPa]
            self.phBot[ilevel] = float(segs[2]) # [hPa]
            self.cc[ilevel] = float(segs[3]) # [0-1]
            self.clw[ilevel] = float(segs[4]) # [kg/kg]
            self.rain[ilevel] = float(segs[5]) # [kg/kg]
            self.ciw[ilevel] = float(segs[6]) # [kg/kg]
            self.snow[ilevel] = float(segs[7]) # [kg/kg]
            
if __name__ == '__main__':
    filename = './jacobian.dat'
    with open(filename, 'r') as fin:
        base_profile = rttov_profile()
        base_profile.readfromfile(fin)
        print(base_profile.metaInfo)
        print(base_profile.viewingConditions)
        print(base_profile.blackBodyCloud)
        print(base_profile.skin)
        print(base_profile.s2m)
        print(base_profile.p)
        print(base_profile.T)
        print(base_profile.q)
        
        base_cld_profile = rttov_cld_profile()
        base_cld_profile.readfromfile(fin) 

        print(base_cld_profile.phTop)
        print(base_cld_profile.phBot)
        print(base_cld_profile.cc)
        print(base_cld_profile.clw)
        print(base_cld_profile.rain)
        print(base_cld_profile.ciw)
        print(base_cld_profile.snow)
        
        base_profile.cld_profile = base_cld_profile
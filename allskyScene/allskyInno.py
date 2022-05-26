'''
Description: store all-sky innovation and plot
Author: Hejun Xie
Date: 2022-05-11 11:12:49
LastEditors: Hejun Xie
LastEditTime: 2022-05-24 13:31:37
'''

# Global imports
import numpy as np
import os
from utils import get_channel_str, chkey

# Local imports
from .read_allskycheck import read_data, get_channel_key, get_domain_key
from .read_innotb import read_data as read_all_data

class AllSkyInno(object):
    def __init__(self, instrument, allsky_innofiles, valid_channels=[3,5,7], total_channels=10):
        self.instrument = instrument
        self.allsky_innofiles = allsky_innofiles
        self.valid_channels = valid_channels
        self.total_channels = total_channels

        # data table keys
        self.datakeys = ['line', 'pos', 'lat', 'lon', 
                        'O', 'BMO', 'BMO(BC)', 'BC',
                        'surf-H', 'surf-T', 
                        'C37obs', 'C37mdl', 'C37sym', 'cldraderr']

        self.raw_data = dict()
        self.domain_merged_data = dict()
        self.timeslots = []
        self.read_allskyInno()
        self.merge_domain()

        # for reference data
        self.expr_dir = None
        self.scatt_data = dict()
        self.direct_data = dict()
        self.singleobs_timeslots = list()
        self.singleobs_channel   = list()
        self.singleobs_latlon    = list()

    def read_allskyInno(self):
        for allsky_innofile in self.allsky_innofiles:
            datetime = os.path.basename(allsky_innofile).split('.')[0][-12:]
            print(datetime)
            self.timeslots.append(datetime)
            self.raw_data[datetime] = read_data(allsky_innofile, len(self.valid_channels))

    def read_ref_data(self, expr_dir):
        self.expr_dir = expr_dir
        for allsky_innofile in self.allsky_innofiles:
            datetime = os.path.basename(allsky_innofile).split('.')[0][-12:]
            direct_innofile = os.path.join(expr_dir, os.path.basename(allsky_innofile).strip('allsky') + 'direct')
            scatt_innofile = os.path.join(expr_dir, os.path.basename(allsky_innofile).strip('allsky') + 'scatt')
            self.scatt_data[datetime] = read_all_data(scatt_innofile, self.total_channels)
            self.direct_data[datetime] = read_all_data(direct_innofile, self.total_channels)
    
    def print_ref_data(self, timeslot, channel, idx):
        channelDataScatt = self.scatt_data[timeslot][chkey(channel-1)]
        channelDataDirect = self.direct_data[timeslot][chkey(channel-1)]
        print('Found reference data in {}'.format(self.expr_dir))

        print('RTTOV-SCATT')
        print('line:{} pos:{}'.format(
            channelDataScatt['line'][idx], channelDataScatt['pos'][idx],
            ))

        print('lat:{:>.2f} lon:{:>.2f} surf-Height:{} surf-Type:{}'.format(
            channelDataScatt['lat'][idx], channelDataScatt['lon'][idx],
            channelDataScatt['surf-height'][idx], channelDataScatt['surf-type'][idx],
            ))

        print('O:{:>.2f} B:{:>.2f} OMB:{:>.2f} OMB(BC):{:>.2f}'.format(
            channelDataScatt['O'][idx], channelDataScatt['B'][idx],
            channelDataScatt['OMB_RAW'][idx], channelDataScatt['OMB_BC'][idx],
            ))
        
        print('C37obs:{:>.3f} C37mdl:{:>.3f} C37sym:{:>.3f} cldraderr:{:>.2f}'.format(
            channelDataScatt['c37obs'][idx], channelDataScatt['c37mdl'][idx], 
            channelDataScatt['c37sym'][idx], channelDataScatt['cldraderr'][idx],
            ))
        
        print('RTTOV-DIRECT')
        print('line:{} pos:{}'.format(
            channelDataDirect['line'][idx], channelDataDirect['pos'][idx],
            ))

        print('lat:{:>.2f} lon:{:>.2f} surf-Height:{} surf-Type:{}'.format(
            channelDataDirect['lat'][idx], channelDataDirect['lon'][idx],
            channelDataDirect['surf-height'][idx], channelDataDirect['surf-type'][idx],
            ))

        print('O:{:>.2f} B:{:>.2f} OMB:{:>.2f} OMB(BC):{:>.2f}'.format(
            channelDataDirect['O'][idx], channelDataDirect['B'][idx],
            channelDataDirect['OMB_RAW'][idx], channelDataDirect['OMB_BC'][idx],
            ))
        
        print('C37obs:{:>.3f} C37mdl:{:>.3f} C37sym:{:>.3f} cldraderr:{:>.2f}'.format(
            channelDataDirect['c37obs'][idx], channelDataDirect['c37mdl'][idx], 
            channelDataDirect['c37sym'][idx], channelDataDirect['cldraderr'][idx],
            ))

    def compare_ref_data(self):
        for timeslot, channel, latlon in zip(self.singleobs_timeslots, 
            self.singleobs_channel, self.singleobs_latlon):
            lat, lon = latlon
            channelData = self.scatt_data[timeslot][chkey(channel)]
            idx = np.argwhere(
                (np.abs(channelData['lat']-lat)<1e-6) &\
                (np.abs(channelData['lon']-lon)<1e-6))[0][0]
            self.print_ref_data(timeslot, channel, idx)    

    def merge_domain(self):
        '''
        Description:
            Merge data in all domains 
        '''

        for timeslot in self.timeslots:
            rawData = self.raw_data[timeslot]

            timeslotData = dict()
            TempPieceData = dict()
            for channel in self.valid_channels:
                timeslotData[get_channel_key(channel)] = dict() # 'datakey'
                TempPieceData[get_channel_key(channel)] = list() # 'piece of data divided by domain'
            self.domain_merged_data[timeslot] = timeslotData

            # Decascading for domains, get TempPieceData
            for domainKey in rawData.keys():
                for channel in self.valid_channels:
                    try:
                        channelDataPiece = rawData[domainKey][get_channel_key(channel)]
                        TempPieceData[get_channel_key(channel)].append(channelDataPiece)
                    except:
                        pass # case if that domain has no valid FOVs for that channel
            
            # merge
            for channel in self.valid_channels:
                TempPieces = TempPieceData[get_channel_key(channel)]
                for TempPiece in TempPieces:
                    for datakey in self.datakeys:
                        value = TempPiece['Data'][datakey]
                        if datakey in timeslotData[get_channel_key(channel)].keys():
                            timeslotData[get_channel_key(channel)][datakey] = \
                                np.hstack((timeslotData[get_channel_key(channel)][datakey], value))
                        else:
                            timeslotData[get_channel_key(channel)][datakey] = value
    
    def print_single_obs(self, timeslot, channel, idx):
        '''
        Print single obs data
        '''
        channelData = self.domain_merged_data[timeslot][get_channel_key(channel)]

        self.singleobs_timeslots.append(timeslot)
        self.singleobs_channel.append(channel)
        self.singleobs_latlon.append([channelData['lat'][idx], channelData['lon'][idx]])

        print('Found active obs in {} for {}!'.format(timeslot, get_channel_key(channel)))

        print('line:{} pos:{}'.format(
            channelData['line'][idx], channelData['pos'][idx],
            ))

        print('lat:{:>.2f} lon:{:>.2f} surf-Height:{} surf-Type:{}'.format(
            channelData['lat'][idx], channelData['lon'][idx],
            channelData['surf-H'][idx], channelData['surf-T'][idx],
            ))

        print('O:{:>.2f} B:{:>.2f} OMB:{:>.2f} OMB(BC):{:>.2f}'.format(
            channelData['O'][idx], channelData['BMO'][idx] + channelData['O'][idx],
            -channelData['BMO'][idx], -channelData['BMO(BC)'][idx],
            ))
        
        print('C37obs:{:>.3f} C37mdl:{:>.3f} C37sym:{:>.3f} cldraderr:{:>.2f}'.format(
            channelData['C37obs'][idx], channelData['C37mdl'][idx], 
            channelData['C37sym'][idx], channelData['cldraderr'][idx],
            ))

    def find_single_obs(self, find_llbox):
        '''
        find llbox = [lonW, lonE, latS, latN]
        '''
        lonW, lonE, latS, latN = find_llbox
        for timeslot in self.timeslots:
            for channel in self.valid_channels:
                channelData = self.domain_merged_data[timeslot][get_channel_key(channel)]
                try:
                    lat = channelData['lat']
                    lon = channelData['lon']
                except:
                    continue # case where it is an empty timeslot
                result = np.argwhere(\
                    (lon > lonW) & (lon < lonE) & \
                    (lat > latS) & (lat < latN))        
                if len(result) != 0:
                    for iresult in range(len(result[0])):
                        idx = result[0][iresult]
                        self.print_single_obs(timeslot, channel, idx)

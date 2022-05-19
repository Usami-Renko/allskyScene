'''
Description: store all-sky innovation and plot
Author: Hejun Xie
Date: 2022-05-11 11:12:49
LastEditors: Hejun Xie
LastEditTime: 2022-05-20 00:47:25
'''

# Global imports
import numpy as np
import os
from utils import get_channel_str

# Local imports
from .read_allskycheck import read_data, get_channel_key, get_domain_key

class AllSkyInno(object):
    def __init__(self, instrument, allsky_innofiles, valid_channels=[3,5,7]):
        self.instrument = instrument
        self.allsky_innofiles = allsky_innofiles
        self.valid_channels = valid_channels

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

    def read_allskyInno(self):
        for allsky_innofile in self.allsky_innofiles:
            datetime = os.path.basename(allsky_innofile).split('.')[0][-12:]
            print(datetime)
            self.timeslots.append(datetime)
            self.raw_data[datetime] = read_data(allsky_innofile)

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

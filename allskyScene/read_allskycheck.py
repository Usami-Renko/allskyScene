'''
Description: read all-sky checked innovations 
Author: Hejun Xie
Date: 2022-05-10 22:14:23
LastEditors: Hejun Xie
LastEditTime: 2022-05-22 19:42:01
'''

import pandas as pd
import numpy as np
from utils import chkey
import re

NUMBER = '[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?'
DomainLine1FMT='PROCESSOR DOMAIN #\s*(\d+)\s* Total FOVs:\s*(\d+)\s*'
QCRound1FMT='\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)'
QCRound2FMT='\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)'
ChannelLine1FMT='CHANNEL #\s*(\d+)\s* Used FOVs:\s*(\d+)\s*'
QCRound3FMT='\s*(\d+)\s+(\d+)\s+(\d+)\s+(\d+)\s+(\d+)'
TabelFMT='\s*({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})\s+({})'.format(\
    NUMBER,NUMBER,NUMBER,NUMBER,NUMBER,NUMBER,NUMBER,\
    NUMBER,NUMBER,NUMBER,NUMBER,NUMBER,NUMBER,NUMBER)

def get_domain_key(PID):
    return 'Domain#{}'.format(PID)

def get_channel_key(CID):
    return 'Channel#{}'.format(CID)

def read_processor_domain_head(fin):
    domainData = dict()
    DomainLine1 = fin.readline()
    if DomainLine1 == '':
        return None
    else:
        p = re.compile(DomainLine1FMT)
        m = p.search(DomainLine1)
        domainData['PID'] = int(m.group(1))
        domainData['TotalFOVs'] = int(m.group(2))
        
        QCRound1Head = fin.readline()
        QCRound1 = fin.readline()
        p = re.compile(QCRound1FMT)
        m = p.search(QCRound1)
        domainData['Round1st'] = dict()
        domainData['Round1st']['Total'] = int(m.group(1))
        domainData['Round1st']['Mixsurface'] = int(m.group(2))
        domainData['Round1st']['Surface'] = int(m.group(3))
        domainData['Round1st']['WindSpeed'] = int(m.group(4))
        domainData['Round1st']['Latitude'] = int(m.group(5))

        QCRound2Head = fin.readline()
        QCRound2 = fin.readline()
        p = re.compile(QCRound2FMT)
        m = p.search(QCRound2)
        domainData['Round2nd'] = dict()
        domainData['Round2nd']['Total'] = int(m.group(1))
        domainData['Round2nd']['c37'] = int(m.group(2))
        domainData['Round2nd']['c37sym'] = int(m.group(3))
        domainData['Round2nd']['c37diff'] = int(m.group(4))
        domainData['RestFOVs'] = domainData['TotalFOVs'] - \
                domainData['Round1st']['Total'] - domainData['Round2nd']['Total']

        return domainData

def read_channel_head(fin):
    channelData = dict()
    ChannelLine1 = fin.readline()
    p = re.compile(ChannelLine1FMT)
    m = p.search(ChannelLine1)
    channelData['CID'] = int(m.group(1))
    channelData['UsedFOVs'] = int(m.group(2))

    QCRound3Head = fin.readline()
    QCRound3 = fin.readline()
    p = re.compile(QCRound3FMT)
    m = p.search(QCRound3)
    channelData['Round3rd'] = dict()
    channelData['Round3rd']['Total'] = int(m.group(1))
    channelData['Round3rd']['RFI'] = int(m.group(2))
    channelData['Round3rd']['Abnormal'] = int(m.group(3))
    channelData['Round3rd']['BgQC_ABS'] = int(m.group(4))
    channelData['Round3rd']['BgQC_STD'] = int(m.group(5))

    return channelData

def read_table(fin, channelData):
    TabelHead = fin.readline()
    nFOVs = channelData['UsedFOVs']

    table = dict()
    channelData['Data'] = table

    table['line']       = np.zeros((nFOVs), dtype='int')
    table['pos']        = np.zeros((nFOVs), dtype='int')
    table['lat']        = np.zeros((nFOVs), dtype='float')
    table['lon']        = np.zeros((nFOVs), dtype='float')
    table['O']          = np.zeros((nFOVs), dtype='float')
    table['BMO']        = np.zeros((nFOVs), dtype='float')
    table['BMO(BC)']    = np.zeros((nFOVs), dtype='float')
    table['BC']         = np.zeros((nFOVs), dtype='float')
    table['surf-H']     = np.zeros((nFOVs), dtype='int')
    table['surf-T']     = np.zeros((nFOVs), dtype='int')
    table['C37obs']     = np.zeros((nFOVs), dtype='float')
    table['C37mdl']     = np.zeros((nFOVs), dtype='float')
    table['C37sym']     = np.zeros((nFOVs), dtype='float')
    table['cldraderr']  = np.zeros((nFOVs), dtype='float')

    for iFOV in range(nFOVs):
        TabelLine = fin.readline()
        p = re.compile(TabelFMT)
        m = p.search(TabelLine)

        table['line'][iFOV]       = int(m.group(1))
        table['pos'][iFOV]        = int(m.group(2))
        table['lat'][iFOV]        = float(m.group(3))
        table['lon'][iFOV]        = float(m.group(4))
        table['O'][iFOV]          = float(m.group(5))
        table['BMO'][iFOV]        = float(m.group(6))
        table['BMO(BC)'][iFOV]    = float(m.group(7))
        table['BC'][iFOV]         = float(m.group(8))
        table['surf-H'][iFOV]     = int(m.group(9))
        table['surf-T'][iFOV]     = int(m.group(10))
        table['C37obs'][iFOV]     = float(m.group(11))
        table['C37mdl'][iFOV]     = float(m.group(12))
        table['C37sym'][iFOV]     = float(m.group(13))
        table['cldraderr'][iFOV]  = float(m.group(14))


def read_data(file, nactivechannel):

    data = dict()

    with open(file, 'r') as fin:

        while True:
            domainData = read_processor_domain_head(fin)
            if domainData is not None:
                data[get_domain_key(domainData['PID'])] = domainData
                if domainData['RestFOVs'] == 0: # No channel info then
                    continue
                for _ in range(nactivechannel): # channel 3, 5, 7
                    channelData = read_channel_head(fin)
                    read_table(fin,channelData)
                    domainData[get_channel_key(channelData['CID'])] = channelData
            else:
                break


    
    return data

if __name__ == '__main__':

    TESTFILE = './singleInst/checkinno/innotbinno_fy3dmwi201906040300.dat_allsky'
    data = read_data(TESTFILE)

    print(data[get_domain_key(163)][get_channel_key(3)])
    print(data[get_domain_key(163)][get_channel_key(5)])
    print(data[get_domain_key(163)][get_channel_key(7)])


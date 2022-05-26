'''
Description: 
Author: Hejun Xie
Date: 2022-05-23 20:35:07
LastEditors: Hejun Xie
LastEditTime: 2022-05-24 13:26:44
'''
'''
Description: Utilities for Innovation Brightness Temperature
Author: Hejun Xie
Date: 2022-04-07 19:03:32
LastEditors: Hejun Xie
LastEditTime: 2022-04-09 20:02:24
'''

import pandas as pd
import numpy as np

def read_domain_head(fin):
    line = fin.readline()
    if line == '':
        return None
    else:
        segs = line.split()
        num_rad = segs[1]
        num_chn = segs[3]
        return (int(num_rad), int(num_chn))

def read_channel_head(fin):
    fin.readline()
    fin.readline()

def read_table(fin, num_rad):
    table = {}

    table['line'] = np.zeros((num_rad), dtype='int')
    table['pos']  = np.zeros((num_rad), dtype='int')
    table['lat']  = np.zeros((num_rad), dtype='float')
    table['lon']  = np.zeros((num_rad), dtype='float')
    table['Tb']  = np.zeros((num_rad), dtype='float')
    table['Tb(xb)-Tb']  = np.zeros((num_rad), dtype='float')
    table['Tb(xb)-Tb-Bias']  = np.zeros((num_rad), dtype='float')
    table['Bias']  = np.zeros((num_rad), dtype='float')
    table['surf-height'] = np.zeros((num_rad), dtype='int')
    table['surf-type']  = np.zeros((num_rad), dtype='int')
    table['effecfrac']  = np.zeros((num_rad), dtype='float')
    table['c37obs'] = np.zeros((num_rad), dtype='float32')
    table['c37mdl'] = np.zeros((num_rad), dtype='float32')
    table['c37sym'] = np.zeros((num_rad), dtype='float32')
    table['cldraderr'] = np.zeros((num_rad), dtype='float32')

    for irad in range(num_rad):
        line = fin.readline()
        segs = line.split()

        try:
            table['line'][irad]           = int(segs[0])
            table['pos'][irad]            = int(segs[1])
            table['lat'][irad]            = float(segs[2])
            table['lon'][irad]            = float(segs[3])
            table['Tb'][irad]             = float(segs[4])
            table['Tb(xb)-Tb'][irad]      = float(segs[5])
            table['Tb(xb)-Tb-Bias'][irad] = float(segs[6])
            table['Bias'][irad]           = float(segs[7])
            table['surf-height'][irad]    = int(segs[8])
            table['surf-type'][irad]      = int(segs[9])
            table['effecfrac'][irad]      = float(segs[10])
            table['c37obs'][irad]         = float(segs[11])
            table['c37mdl'][irad]         = float(segs[12])
            table['c37sym'][irad]         = float(segs[13])
            table['cldraderr'][irad]      = float(segs[14])

        except ValueError:
            table['Tb'][irad]             = np.nan
            table['Tb(xb)-Tb'][irad]      = np.nan
            table['Tb(xb)-Tb-Bias'][irad] = np.nan
            table['Bias'][irad]           = np.nan
            table['effecfrac'][irad]      = np.nan
            table['c37obs'][irad]         = np.nan
            table['c37mdl'][irad]         = np.nan
            table['c37sym'][irad]         = np.nan
            table['cldraderr'][irad]      = np.nan
            continue

        if table['Tb'][irad] < 0.0:
            table['Tb'][irad] = np.nan
            table['Tb(xb)-Tb'][irad] = np.nan
            table['Tb(xb)-Tb-Bias'][irad] = np.nan
            table['Bias'][irad] = np.nan
        
    table['B'] = table['Tb(xb)-Tb'] + table['Tb']
    table['O'] = table['Tb']
    table['OMB_RAW'] = - table['Tb(xb)-Tb']
    table['OMB_BC'] = - table['Tb(xb)-Tb-Bias']
    table['BC'] = - table['Bias']

    return table

def read_data(file, num_chn):

    data = {}
    for i_chn in range(num_chn):
        data['ch{}'.format(i_chn+1)] = {}

    with open(file, 'r') as fin:
        while True:
            pack = read_domain_head(fin)
            if pack is not None:
                num_rad, num_chn = pack
                for i_chn in range(num_chn):
                    ch_str = 'ch{}'.format(i_chn+1)
                    read_channel_head(fin)
                    table = read_table(fin, num_rad)
                    for key, value in table.items():
                        if key in data[ch_str].keys():
                            data[ch_str][key] = np.hstack((data[ch_str][key], value))
                        else:
                            data[ch_str][key] = value

            else:
                break
    
    return data
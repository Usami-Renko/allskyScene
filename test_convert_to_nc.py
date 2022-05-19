'''
Description: test transfer from GRAPES modelvar to nc
Author: Hejun Xie
Date: 2020-10-28 16:33:40
LastEditors: Hejun Xie
LastEditTime: 2022-05-10 19:12:28
'''

# unit test import
import sys
sys.path.append('/home/xhj/wkspcs/Radar-Operator/ZJU_AERO/')

from ZJU_AERO.nwp.grapes import convert_to_nc_general_ctl
from ZJU_AERO.nwp.grapes import convert_to_nc_specific_ctl
import glob
import os

type_var = {
    'dxa': ['pi', 'th', 'u', 'v', 'w', 'q', 'qr', 't', 'p', 'rh'],
    'xa':  ['pi', 'th', 'u', 'v', 'w', 'q', 'qr', 't', 'p', 'rh'],
    'grapes_input': ['pi', 'th', 'u', 'v', 'w', 'qv', 'qc', 'qr', 'qi', 'qs', 'qg', 'pcf'],
}

def convert_ctls(ctl_matches, datatype):

    CTLFILES = glob.glob(ctl_matches)

    NCS = []
    for CTLFILE in CTLFILES:
        ctl_basename = os.path.basename(CTLFILE)
        ctl_dirname = os.path.dirname(CTLFILE)

        nc_basename = '{}.nc'.format(os.path.splitext(ctl_basename)[0])
        NC = os.path.join(ctl_dirname, nc_basename)
        NCS.append(NC)

    for CTLFILE, NC in zip(CTLFILES, NCS):
        print(NC)
        convert_to_nc_specific_ctl(CTLFILE, NC, var=type_var[datatype], v_grid_type='gi')


if __name__ == "__main__":

    ctlfiles = sys.argv[1]
    datatype = sys.argv[2]
    print(ctlfiles)
    print(datatype)

    convert_ctls(ctlfiles, datatype)

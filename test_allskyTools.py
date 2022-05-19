'''
Description: test module AllSkyOverview
Author: Hejun Xie
Date: 2022-05-20 00:14:16
LastEditors: Hejun Xie
LastEditTime: 2022-05-20 01:06:58
'''

# Global imports
import os
import glob
import numpy as np
import xarray as xr

from allskyScene import AllSkyInno, \
    AGRI_data, MWRI_data, \
    plot_allsky_level, interp2plevel, \
    get_dxa_dataset, get_xb_dataset, get_grapesinput_dataset 

FORCE_INTERP = False

# plot settings
# PLOT_LEVELS = [1000.,850.,700.,500.,300.,200.]
PLOT_LEVELS = [850.]
# PLOT_LEVELS = ['bottom']

'''
1. READ INNOVATION
'''
expr_dir = './singleInstYanhua'
allsky_innofiles = glob.glob('{}/checkinno/innotbinno_fy3dmwi????????????.dat_allsky'.format(expr_dir))
asi = AllSkyInno('mwri', allsky_innofiles)

'''
2. Read satellite data
'''
AGRI_NDISK_files = glob.glob('./AGRI/FY4A-_AGRI--_N_DISK_1047E_L1-_FDI-_MULT_NOM_20210719180000_??????????????_4000M_V0001.HDF')
agri = AGRI_data(AGRI_NDISK_files)
agri.read_agri_ndisk(ir_channel='C12')

MWRI_SWATH_files = glob.glob('./MWRI/2021071915-21/FY3D_MWRI?_GBAL_L1_20210719_????_010KM_MS.HDF')
mwri = MWRI_data(MWRI_SWATH_files)
mwri.read_mwri_swaths(mwri_channel=6)

'''
3. PLOT SINGLE CHANNEL FILE
'''
ch_nos = [5]
for ch_no in ch_nos:
    expr_dir = './singleChannel{}Yanhua'.format(ch_no)
    print(expr_dir)
    dxa_file = '{}/binary/dxa.nc'.format(expr_dir)
    xb_file =  '{}/binary/xb.nc'.format(expr_dir)
    grapesinput_file = '{}/binary/grapes_input.in.nc'.format(expr_dir)
    dxa_plevel_file = "{}/binary/dxa_plevel.nc".format(expr_dir)
    xb_plevel_file = "{}/binary/xb_plevel.nc".format(expr_dir)
    
    ds_dxa = get_dxa_dataset(dxa_file)
    ds_xb = get_xb_dataset(xb_file)
    ds_grapesinput = get_grapesinput_dataset(grapesinput_file)
    
    if FORCE_INTERP:
        ds_dxa_plevel = interp2plevel(ds_dxa, ds_dxa, ds_xb.data_vars['Pressure'])
        ds_dxa_plevel.to_netcdf(dxa_plevel_file)
        ds_xb_plevel = interp2plevel(ds_xb, ds_dxa, ds_xb.data_vars['Pressure'])
        ds_xb_plevel.to_netcdf(xb_plevel_file)
    else:
        ds_dxa_plevel = xr.open_dataset(dxa_plevel_file)
        ds_xb_plevel = xr.open_dataset(xb_plevel_file)
    
    for level in PLOT_LEVELS:
        plot_allsky_level(asi, ch_no, 
            ds_dxa, ds_dxa_plevel, 
            ds_xb, ds_xb_plevel, ds_grapesinput, 
            agri, mwri,
            level=level)
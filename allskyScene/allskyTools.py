'''
Description: plot dxa
Author: Hejun Xie
Date: 2022-05-10 16:53:55
LastEditors: Hejun Xie
LastEditTime: 2022-05-20 20:49:14
'''

'''
Global imports
'''
import os
import sys
import glob
import copy

from itertools import product

import warnings
warnings.simplefilter("ignore", UserWarning)

import numpy as np
import netCDF4 as nc
import xarray as xr
from scipy import interpolate

# matplotlib
import matplotlib.pyplot as plt
from matplotlib import colors
import matplotlib as mpl

# proplot
rcParams_copy = copy.copy(mpl.rcParams)  
import proplot as pplt
pplt.rc.update(**rcParams_copy)

# cartopy
import cartopy.crs as ccrs

# pyUtils
from utils import get_channel_str, get_vminvmax

'''
Local imports
'''
from .read_allskycheck import get_channel_key

'''
User Settings
'''
# Coordinates (x,y) --> (lon, lat) 
# x = [-180,180.], x = lon 
data_projection = ccrs.PlateCarree(central_longitude=0.0) 
plot_projection = ccrs.PlateCarree(central_longitude=180.0)   

# Manual levels for contour reserved
Manual_RH_LEVELS = False
RH_LEVELS_CF = np.linspace(-8.0,8.0,21)
RH_LEVELS_CS = np.linspace(-8.0,8.0,11)
RH_TICKS = [-8.0, -6.0, -4.0, -2.0, 0.0, 2.0, 4.0, 6.0, 8.0]

Manual_T_LEVELS = False
T_LEVELS_CF = np.linspace(-0.20,0.20,21)
T_LEVELS_CS = np.linspace(-0.20,0.20,11)
T_TICKS = [-0.2, -0.10, 0.0, 0.10, 0.20]

# Hydrometeor Manual levels for contour
Q_LEVELS = [1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2]
Q_LEVEL_LABELS = [r'$1.0e^{-7}$', r'$1.0e^{-6}$', r'$1.0e^{-5}$', 
                r'$1.0e^{-4}$', r'$1.0e^{-3}$', r'$1.0e^{-2}$']

# interped plevels
PLEVEL = np.array([1000.,850.,700.,500.,300.,200.])

LAYOUT_DEBUG = False
BBOX_CLIP = True

class AllSkyOverview(object):
    '''
    Description:
        A class for drawing the overview of all-sky assimilation
    '''
    def __init__(self, asi, ch_no,
                ds_dxa, ds_dxa_plevel,
                ds_xb, ds_xb_plevel, ds_grapesinput, 
                agri, mwri):
        self.hydro_overlay_jobs = ['AGRI_IR']
        self.analy_incre_jobs = ['RH']
        self.level_jobs = [850.]
        self.region_jobs = ['Global']
        self.layout_settings = dict()

        # save data
        self.asi = asi
        self.ch_no = ch_no
        self.ds_dxa = ds_dxa
        self.ds_dxa_plevel = ds_dxa_plevel
        self.ds_xb = ds_xb
        self.ds_xb_plevel = ds_xb_plevel
        self.agri = agri
        self.mwri = mwri

    def update_layout_settings(self, region):
        '''
        Description : layout setting for matplotlib plot for different region
        '''    
        if region == 'Global':
            # self.layout_settings['llbox'] = [0.,360.,-90.,90.]
            self.layout_settings['llbox'] = [0.,360.,-60.,60.]
            self.layout_settings['gridlineSpacing'] = 30.
            self.layout_settings['TwoColorBarfigsize'] = (22,7.5)
            self.layout_settings['colorbarRightPad'] = 0.01
            self.layout_settings['colorbarLeftPad']  = 0.03
            self.layout_settings['colorbarShrink'] = 0.67
            self.layout_settings['innoMakerSize'] = 8
            self.layout_settings['innoMakerEdge'] = 1/24
            self.layout_settings['contourLineWidth'] = 1.0
            self.layout_settings['topoHatchDensity'] = 4
            self.layout_settings['streamLineWidth'] = 3.0
            self.layout_settings['streamLineDensity'] = 3.0
            self.layout_settings['coastlineWidth'] = 1.5
        elif region in ['EastAsia', 'NorthIndianOcean']:
            if region == 'EastAsia':
                self.layout_settings['llbox'] = [80., 160.,  0.,  50.]
            elif region == 'NorthIndianOcean':
                self.layout_settings['llbox'] = [40., 100.,-10., 27.5]
            self.layout_settings['gridlineSpacing'] = 10.
            self.layout_settings['TwoColorBarfigsize'] = (17,10)
            self.layout_settings['colorbarRightPad'] = 0.01
            self.layout_settings['colorbarLeftPad']  = 0.03
            self.layout_settings['colorbarShrink'] = 0.73
            self.layout_settings['innoMakerSize'] = 30
            self.layout_settings['innoMakerEdge'] = 1/10
            self.layout_settings['contourLineWidth'] = 1.5
            self.layout_settings['topoHatchDensity'] = 2
            self.layout_settings['streamLineWidth'] = 4.0
            self.layout_settings['streamLineDensity'] = 2.0
            self.layout_settings['coastlineWidth'] = 2.0
        else:
            raise ValueError('Invalid region {}'.format(region))
    
    def assign_jobs(self,
        region_jobs,
        level_jobs,
        analy_incre_jobs,
        hydro_overlay_jobs,
        ):
        '''
        Assign jobs here
        '''
        self.region_jobs = region_jobs
        self.level_jobs = level_jobs
        self.analy_incre_jobs = analy_incre_jobs
        self.hydro_overlay_jobs = hydro_overlay_jobs
        
    def do_jobs(self):
        '''
        Do plot jobs here
        '''
        for region in self.region_jobs:
            self.update_layout_settings(region)
            for level, analy_incre, hydro_overlay in \
                product(self.level_jobs, self.analy_incre_jobs, self.hydro_overlay_jobs):
                self.plot_allsky_level(region, level, analy_incre, hydro_overlay)

    def plot_allsky_level(self, region, level, analy_incre, hydro_overlay):
        '''
        plot allsky assimilation overview on a given level
        1) type(level) == int: model layer specification
        2) type(level) == float: pressure layer specification
        3) type(level) == 'bottom' bottom model layer
        '''
        if level == 'bottom':
            level = 0
        
        LayerTag = '{}hPa'.format(level) if isinstance(level, float) \
            else 'Level{}'.format(level)
        FIGNAME = 'd{}_{}_ch{}_{}_{}.png'.format(analy_incre, LayerTag, self.ch_no, hydro_overlay, region)
        print(FIGNAME)
        
        fig = plt.figure(figsize=self.layout_settings['TwoColorBarfigsize'])
        ax = self.draw_cartopy_map(fig, 111)
        plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0, wspace=0.0, hspace=0.0)

        if not LAYOUT_DEBUG:
            '''
            a). Plot analysis increment contour [zorder = 5] 
            '''
            if analy_incre == 'RH':
                self.plot_incre_RH(ax, fig, level, CBlocation='left', zorder=5)
            elif analy_incre == 'T':
                self.plot_incre_T(ax, fig, level, CBlocation='left', zorder=5)
            else:
                raise ValueError('Invalid analy_incre {}'.format(analy_incre))

            '''
            b). Plot cloud proxy [zorder = 10]
            '''
            # CHANNEL 12 of FY4A-AGRI (IR): 10.8um Window channel
            if hydro_overlay == 'AGRI_IR':
                self.plot_agri_ir(ax, fig, CBlocation='right', zorder=10)
            # CHANNEL 7 of FY3D-MWRI: 37GHz V pol. Window channel 
            elif hydro_overlay == 'MWRI':
                self.plot_mwri(ax, fig, CBlocation='right', zorder=10)
            # MODEL hydrometeor plot: Qr column Maximum [kg/kg]
            elif hydro_overlay == 'HYDRO':
                self.plot_hydro(ax, fig, CBlocation='right', zorder=10)
            else:
                raise ValueError('Invalid hydro_overlay {}'.format(hydro_overlay))

            '''
            c). Plot stream line [zorder = 12]
            '''
            self.plot_streamline(ax, fig, level, zorder=12)

            '''
            d). Plot active observation [zorder = 15]
            '''
            self.plot_active_obs(ax, fig, CBlocation=None, zorder=15)
            
        '''
        e). Title, Annotation and Output
        '''
        LayerText = r'$\textbf{' + LayerTag + '}$'
        InstTimeText = r'$\textbf{'+get_channel_str(self.asi.instrument, self.ch_no-1)+ \
            ' [{}-{}]'.format(self.asi.timeslots[0], self.asi.timeslots[-1])+'}$'
        ax.text(0.0, 1.0, LayerText, 
            ha='left', va='bottom', size=24, transform=ax.transAxes)
        ax.text(1.0, 1.0, InstTimeText, 
            ha='right', va='bottom', size=24, transform=ax.transAxes)
        plt.savefig(FIGNAME, bbox_inches='tight') if BBOX_CLIP \
            else plt.savefig(FIGNAME)
        plt.close()

    def _plot_colorbar(self, plotFig, plotAx, plotStuff, label, location=None, ticks='auto'):
        '''
        Description:
            Plot colorbar for a given (plot plotFig, plotAx, plotStuff)
        Params: 
            1. label (a string)
            2. location: None, left or right
            3. ticks: a tick list or 'auto'
        '''
        if location == None:
            return

        if location == 'left':
            CB = plotFig.colorbar(plotStuff, ax=[plotAx], shrink=self.layout_settings['colorbarShrink'], 
                pad=self.layout_settings['colorbarLeftPad'], location='left')
        elif location == 'right':
            CB = plotFig.colorbar(plotStuff, ax=[plotAx], shrink=self.layout_settings['colorbarShrink'], 
                pad=self.layout_settings['colorbarRightPad'], location='right')
        else:
            raise ValueError('Invalid location {}'.format(location))
        CB.ax.set_ylabel(label, fontsize=20)

        if isinstance(ticks, list) or isinstance(ticks, np.ndarray):
            CB.set_ticks(ticks)
        elif isinstance(ticks, str) and ticks == 'auto':
            pass
        else:
            raise ValueError('Invalid ticks type:{}'.format(type(ticks)))

    def _plot_nan_hatch(self, ax, x, y, data, transform=None, 
        color='brown', density=4, zorder=5):
        '''
        Description: Fill nan of data with hatches of given color & given transform
        '''
        if transform == None:
            transform = ax.transData
        # there must be two levels
        hatch = ax.contourf(x, y, np.isnan(data), 
            levels=[.5,1.5], colors='none',
            hatches=[density*'/', density*'/'], transform=transform,
            zorder=zorder)

        hatch.collections[0].set_edgecolor(color)
        hatch.collections[0].set_linewidth(0.)

    def _select_layer(self, ds, ds_plevel, level, label):
        '''
        Params:
            type(level) == int: model layer specification
            type(level) == float: pressure layer specification
        '''
        if isinstance(level, int):
            return ds.data_vars[label].isel(level=level) if 'level' in ds.data_vars[label].dims \
                else ds.data_vars[label].isel(level_uv=level)
        elif isinstance(level, float):
            return ds_plevel.data_vars[label].sel(Pressure=level, method='nearest')
        else:
            raise ValueError('Invalid level specification: {}'.format(level))

    def _convert_lon(self, lon):
        '''
        Description:
            Put logitude to [-180, 180]
            to avoid the bug of gridlines & fix stride over bug
            can be viewed as the reverse processof monotonicalize
        '''
        return (lon + 180) % 360 - 180

    def draw_cartopy_map(self, fig, ax_no):
        import cartopy.crs as ccrs
        import matplotlib.ticker as mticker
        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter

        Delt=1e-6
        dashdotdotted = (0, (3, 5, 1, 5, 1, 5))

        # locatations for GridLine and ticks 
        longridlocs = np.arange(self.layout_settings['llbox'][0], self.layout_settings['llbox'][1]+Delt, 
            self.layout_settings['gridlineSpacing'])
        latgridlocs = np.arange(self.layout_settings['llbox'][2], self.layout_settings['llbox'][3]+Delt, 
            self.layout_settings['gridlineSpacing'])

        # Attention: extent cannnot stride over the antilongitude of central longtitude!
        extent = [self.layout_settings['llbox'][0], self.layout_settings['llbox'][1]-Delt, 
            self.layout_settings['llbox'][2], self.layout_settings['llbox'][3]] 

        ax = fig.add_subplot(ax_no, projection=plot_projection)
        ax.set_xticks(longridlocs, crs=data_projection)
        ax.set_yticks(latgridlocs, crs=data_projection)
        ax.xaxis.set_major_formatter(LongitudeFormatter())
        ax.yaxis.set_major_formatter(LatitudeFormatter())
        ax.minorticks_off()
        ax.set_extent(extent, crs=data_projection)
        
        # Attention: gridline parameter xlocs only accept xlocs with [-180., 180.]
        gl = ax.gridlines(draw_labels=False, 
            xlocs=self._convert_lon(longridlocs), ylocs=latgridlocs,
            linewidth=1.5, linestyle=dashdotdotted, color='k', alpha=0.7)

        ax.coastlines(resolution='50m', color='k', alpha=0.8, 
            lw=self.layout_settings['coastlineWidth'])

        return ax

    def plot_active_obs(self, ax, fig, CBlocation=None, zorder=15):
        '''
        plot active observation
        '''
        for itimeslot,timeslot in enumerate(self.asi.timeslots):
            data = self.asi.domain_merged_data[timeslot][get_channel_key(self.ch_no)]
            try:
                lonTag, latTag = data['lon'], data['lat']
                OMB = (- data['BMO(BC)']) / data['cldraderr']
            except KeyError:
                print('Empty timeslot {}'.format(timeslot))
                continue
            
            im = ax.scatter(lonTag, latTag, marker='o', edgecolor='k', 
                s=self.layout_settings['innoMakerSize'], c=OMB, 
                linewidths=itimeslot*self.layout_settings['innoMakerEdge']+0.2,
                cmap='rainbow', vmin = -2.0, vmax = 2.0, transform=data_projection, zorder=zorder)

        self._plot_colorbar(fig, ax, im, 
            label=r'Normalized FG departure [K]',
            location=CBlocation)

    def plot_hydro(self, ax, fig, CBlocation=None, zorder=10):
        '''
        PLOT HYDRO OF GRAPESINPUT
        now implemented as QR column maximum 
        '''
        lon, lat = self.ds_grapesinput.coords['lon'], self.ds_grapesinput.coords['lat']
        
        hydro = np.max(self.ds_grapesinput.data_vars['QR'].data, axis=0)
        
        plot_hydro = np.ma.masked_array(hydro, hydro<Q_LEVELS[0])
        norm = colors.LogNorm(vmin=Q_LEVELS[0], vmax=Q_LEVELS[-1])
        cmap = pplt.Colormap('rainbow', alpha=(0.2, 1.0))
        
        im = ax.imshow(plot_hydro, cmap=cmap, norm=norm, 
            origin='lower', extent=[np.min(lon), np.max(lon), np.min(lat), np.max(lat)],
            transform=data_projection, zorder=zorder)
        
        self._plot_colorbar(fig, ax, im, 
            label=r'QR Column Maximum [kg/kg]',
            location=CBlocation)

    def plot_agri_ir(self, ax, fig, CBlocation=None, zorder=10):
        '''
        PLOT AGRI IR FULL DISK
        '''
        agri = self.agri
        # get agri data crs
        agri_crs    = agri.area_def.to_cartopy_crs()
    
        # spectral [-90, -40] Greys [-40, 30] [Celsius Degree]
        rbtop = pplt.Colormap('spectral', 'Greys_r', ratios=(5, 7), name='rbtop')

        # make clear sky transparent (do not modify original data)
        tb_ir_plot = np.where(agri.tb_ir>0, np.nan, agri.tb_ir) 
        
        im = plt.imshow(tb_ir_plot, origin='upper', 
            cmap=rbtop, vmin=-90, vmax=30, 
            extent=agri_crs.bounds, transform=agri_crs, 
            alpha=0.5, zorder=zorder)

        self._plot_colorbar(fig, ax, im, 
            label=r'Brightness Temperature [$^{\circ}C$]',
            location=CBlocation)

    def plot_mwri(self, ax, fig, CBlocation=None, zorder=10):
        '''
        PLOT MWRI SWATHS
        '''    
        mwri = self.mwri
        # Plot MWRI Swaths
        safe_pad = 1.0
        cmap = pplt.Colormap('rainbow', alpha=(0.0, 1.0))
        for divide_meridian, filename, BT, lat, lon in \
            zip(mwri.divide_meridians, mwri.filenames, mwri.BTs, mwri.lats, mwri.lons):
            print('plot {}'.format(filename))
            BT_left = np.where(lon<divide_meridian-safe_pad, BT, np.nan)
            BT_right = np.where(lon>divide_meridian+safe_pad, BT, np.nan)
            pm = plt.pcolormesh(lon, lat, BT_left, cmap=cmap, vmin=190., vmax=270., shading='auto', 
                transform=data_projection, zorder=zorder)
            pm = plt.pcolormesh(lon, lat, BT_right, cmap=cmap, vmin=190., vmax=270., shading='auto', 
                transform=data_projection, zorder=zorder)

        self._plot_colorbar(fig, ax, pm, 
            label=r'Brightness Temperature [K]',
            location=CBlocation)

    def plot_streamline(self, ax, fig, level, zorder=12):
        '''
        Plot stream line from model
        '''
        lon, lat = self.ds_xb_plevel.coords['lon'], self.ds_xb_plevel.coords['lat']
        U = self._select_layer(self.ds_xb, self.ds_xb_plevel,
                level=level, label='U wind')
        V = self._select_layer(self.ds_xb, self.ds_xb_plevel,
                level=level, label='V wind')
        
        speed = np.sqrt(U.data**2+V.data**2)
        lw = speed / np.nanmax(speed)
        sp = ax.streamplot(lon, lat, U.data, V.data, 
            color='r', linewidth=self.layout_settings['streamLineWidth']*lw, 
            density=self.layout_settings['streamLineDensity'],
            transform=data_projection,
            zorder=zorder)

    def plot_incre_RH(self, ax, fig, level, CBlocation=None, zorder=5):
        lon, lat = self.ds_dxa.coords['lon'], self.ds_dxa.coords['lat']
        dRH = self._select_layer(self.ds_dxa, self.ds_dxa_plevel,
                level=level, label='Increment Relative Humidity')
        
        vstage = 2.0

        if not Manual_RH_LEVELS:
            dRHmin, dRHmax = get_vminvmax(dRH.data.flatten(), vstage=vstage, ign_perc=0.5, lsymmetric=True)
            print('dRHmin={:>.1f}, dRHmax={:>.1f}'.format(dRHmin, dRHmax))
            RH_LEVELS_CF = np.linspace(dRHmin, dRHmax, 21)
            RH_LEVELS_CS = np.linspace(dRHmin, dRHmax, 11)
            RH_TICKS = np.arange(dRHmin, dRHmax+1e-6, 2.0)

        norm = colors.Normalize(vmin=RH_TICKS[0], vmax=RH_TICKS[-1])
        CF = ax.contourf(lon, lat, dRH, RH_LEVELS_CF, cmap='BrBG', norm=norm, extend='both', transform=data_projection)    
        CS = ax.contour(lon, lat, dRH, RH_LEVELS_CS, colors=('k',), 
            linewidths=(self.layout_settings['contourLineWidth'],), transform=data_projection)    
        ax.clabel(CS, fmt='%2.1f', colors='k', fontsize=10)

        self._plot_nan_hatch(ax, lon, lat, dRH, transform=data_projection, 
            color='brown', density=self.layout_settings['topoHatchDensity'], zorder=zorder)

        self._plot_colorbar(fig, ax, CF, 
            label=r'Analysis Increment RH [\%]',
            location='left', 
            ticks=RH_TICKS)

    def plot_incre_T(self, ax, fig, level, CBlocation=None, zorder=5):
        lon, lat = self.ds_dxa.coords['lon'], self.ds_dxa.coords['lat']
        dT = self._select_layer(self.ds_dxa, self.ds_dxa_plevel,
                level=level, label='Increment Temperature')

        vstage = 0.2

        if not Manual_T_LEVELS:
            dTmin, dTmax = get_vminvmax(dT.data.flatten(), vstage=vstage, ign_perc=0.0, lsymmetric=True)
            print('dTmin={:>.1f}, dTmax={:>.1f}'.format(dTmin, dTmax))
            T_LEVELS_CF = np.linspace(dTmin, dTmax, 21)
            T_LEVELS_CS = np.linspace(dTmin, dTmax, 11)
            T_TICKS = np.arange(dTmin, dTmax+1e-6, 0.2)

        norm = colors.Normalize(vmin=T_TICKS[0], vmax=T_TICKS[-1])
        CF = ax.contourf(lon, lat, dT, T_LEVELS_CF, cmap='RdBu_r', norm=norm, extend='both', transform=data_projection)    
        CS = ax.contour(lon, lat, dT, T_LEVELS_CS, colors=('k',), 
            linewidths=(self.layout_settings['contourLineWidth'],), transform=data_projection)
        ax.clabel(CS, fmt='%2.1f', colors='k', fontsize=10)

        self._plot_nan_hatch(ax, lon, lat, dT, transform=data_projection, 
            color='brown', density=self.layout_settings['topoHatchDensity'], zorder=zorder)

        self._plot_colorbar(fig, ax, CF, 
            label=r'Analysis Increment Temperature [K]',
            location='left', 
            ticks=T_TICKS)


def get_dxa_dataset(DXAFILE):
    ds = xr.open_dataset(DXAFILE)
    new_ds = xr.Dataset(
        {
            "Increment Temperature": (["level", "y", "x"], ds.data_vars['t'].data[0,...]),
            "Increment Relative Humidity": (["level", "y", "x"], ds.data_vars['rh'].data[0,...]),
            "Increment Pressure" : (["level", "y", "x"], ds.data_vars['p'].data[0,...]),
        },
        coords={
            "lon" : (["x"], ds.data_vars['longitudes'].data),
            "lat" : (["y"], ds.data_vars['latitudes'].data),
        }
    )
    ds.close()
    return new_ds

def get_xb_dataset(XBFILE):
    ds = xr.open_dataset(XBFILE)

    def interp_uv(raw_uv, p):
        interped_uv = np.zeros(p.shape, dtype='float32')
        interped_uv[1:-1,...]   = (raw_uv[1:,...] + raw_uv[:-1,...]) / 2.0
        interped_uv[0,...]      = raw_uv[0,...]
        interped_uv[-1,...]     = raw_uv[-1,...]
        return interped_uv

    interped_u = interp_uv(ds.data_vars['u'].data[0,...], ds.data_vars['p'].data[0,...])
    interped_v = interp_uv(ds.data_vars['v'].data[0,...], ds.data_vars['p'].data[0,...])
    
    new_ds = xr.Dataset(
        {
            "Temperature": (["level", "y", "x"], ds.data_vars['t'].data[0,...]),
            "Relative Humidity": (["level", "y", "x"], ds.data_vars['rh'].data[0,...]),
            "Pressure" : (["level", "y", "x"], ds.data_vars['p'].data[0,...]),
            'U wind': (["level", "y", "x"], interped_u),
            'V wind': (["level", "y", "x"], interped_v),
        },
        coords={
            "lon" : (["x"], ds.data_vars['longitudes'].data),
            "lat" : (["y"], ds.data_vars['latitudes'].data),
        }
    )
    ds.close()
    return new_ds

def get_grapesinput_dataset(GRAPESINPUT):
    ds = xr.open_dataset(GRAPESINPUT)
    new_ds = xr.Dataset(
        {
            "QR": (["level", "y", "x"], ds.data_vars['qr'].data[0,...]),
        },
        coords={
            "lon" : (["x"], ds.data_vars['longitudes'].data),
            "lat" : (["y"], ds.data_vars['latitudes'].data),
        }
    )
    ds.close()
    return new_ds

def horiz_nearest_interp(ds_out, var_in):
    '''
    Make sure that ds_out and var_in has the same vertical levels
    ds_out is just a template for var_out
    '''
    var_out = np.zeros((len(ds_out.coords['level'].data), 
        len(ds_out.coords['lat'].data), len(ds_out.coords['lon'].data)), dtype='float32')
    for ilon, lon in enumerate(ds_out.coords['lon'].data):
        jlon = np.argmin(np.abs(var_in.coords['lon'].data - lon))
        for ilat, lat in enumerate(ds_out.coords['lat'].data):
            jlat = np.argmin(np.abs(var_in.coords['lat'].data - lat))
            var_out[:,ilat,ilon] = var_in.data[:,jlat,jlon]
    return var_out

def interp2plevel(ds, ds_out, pressure):
    '''
    Description:
        interpolate the data of ds to horizontal grid of ds_grid and pressure
        ds_out is just a template to make horizontal nearest interpolation
    '''
    '''
    1. A nearest interpolation xa grids to dxa grids
    '''
    pressure_on_ds_out_grids = horiz_nearest_interp(ds_out, pressure)
    
    '''
    2. Generate new dataset
    '''
    new_data = dict()
    for var in ds.data_vars.keys():
        if var != 'Pressure':
            new_data[var] = (["Pressure", 'y', 'x'], 
                np.full((len(PLEVEL), len(ds_out.coords['lat'].data), len(ds_out.coords['lon'].data)), np.nan, dtype='float32'))
    new_coords = copy.deepcopy(ds_out.coords)
    new_coords['Pressure'] = PLEVEL 
    new_ds = xr.Dataset(new_data, coords=new_coords)

    '''
    3. Interp ds to ds_out grids in horizontal, then interpolate to pressure levels
    '''   
    for var in new_ds.data_vars.keys():
        data_on_ds_out_grids = horiz_nearest_interp(ds_out, ds.data_vars[var])
        for ilon in range(len(ds_out.coords['lon'].data)):
            for ilat in range(len(ds_out.coords['lat'].data)):
                x = pressure_on_ds_out_grids[::-1,ilat,ilon]    # from descending to ascending
                y = data_on_ds_out_grids[::-1,ilat,ilon]  # from descending to ascending
                f = interpolate.interp1d(x, y, kind='linear', 
                    bounds_error=False, fill_value=np.nan,  # no extrapolating
                    assume_sorted=True)                     # speed-up by canceling ascending check
                new_ds.data_vars[var][:,ilat,ilon] = f(PLEVEL)    
    return new_ds


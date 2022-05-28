'''
Description: plot jacobian of single observation
Author: Hejun Xie
Date: 2022-05-27 15:56:31
LastEditors: Hejun Xie
LastEditTime: 2022-05-28 21:34:35
'''

# Global imports
import numpy as np
import proplot as pplt

from read_jacobian import rttov_profile, rttov_cld_profile
from utils import get_channel_str

directory = './Jacobian/profile3/'
JACOBIANFILE_allsky = '{}/jacobian_allsky.dat'.format(directory)
JACOBIANFILE_clearsky = '{}/jacobian_clearsky.dat'.format(directory)
nchannels = 10
PLOT_CHANNELS = [0,1,2,3,4,5,6,7,8,9]
# PLOT_CHANNELS = [4] # 23.8GHz V pol.

c0 = 'k'
c1 = 'red8'
c2 = 'blue8'
c3 = 'yellow8'
c4 = 'green8'
c5 = 'cyan8'

def plot_jacobian_profile(baseProf, JacobianProfs, JacobianClearProfs, pic, ch_no):
    '''
    Plot Jacobian
    '''
    JacobianProf = JacobianProfs[ch_no]

    # T, q, cc, clw/rain, ciw/snow 
    pplt.rc.update({'meta.width': 1, 'label.weight': 'bold', 'tick.labelweight': 'bold', 
    'fontname': 'Source Sans Pro'})

    handles_clearsky = []
    handles_allsky = []
    
    array = [
        [1,2,0],
        [3,4,5],
    ]
    
    fig, axes = pplt.subplots(array, 
        sharex=False, spanx=False, 
        sharey=True,  spany=True,
        figwidth=8., figheight=6.,
        wspace=3.3)
    
    axes.format(
        abc=r'\textbf{[A].}', abcloc='ul',
        xtickdir='in', ytickdir='in',
        ylim=(1000., 200.), yreverse=True, yscale='log', 
        yticks=(1000., 850., 700., 500., 300., 200.), ytickminor=False,
        xgridminor=True, xticklabelpad=5, xtickminor=True,
        )
    
    '''
    1. T & JacobianT
    '''
    axes[0].linex(baseProf.p, baseProf.T, color=c0)

    ox0 = axes[0].altx(color=c1, label=r'Jacobian Temperature $ \partial{T_b} / \partial{T} $ [K/K]', lw=1)
    art_allsky = ox0.linex(baseProf.p, JacobianProf.T, color=c1, ls='-',
        label=r'$ \partial{T_b} / \partial{T} $')
    art_clearsky = ox0.linex(baseProf.p, JacobianClearProf.T, color=c1, ls='--',
        label=r'$ \partial{T_b} / \partial{T} $')
    handles_allsky.extend(art_allsky)
    handles_clearsky.extend(art_clearsky)

    axes[0].format(
        urtitle=r'\textbf{T}',
        xlabel=r'Temperature [K]',
        xlim=[200, 300],
        yformatter='{x:.0f}hPa',
        xticks=(200, 220, 240., 260., 280., 300), 
    )
    
    ox0.format(
        xlim=(-2e-1, 2e-1), 
        xformatter='sci', xscale='symlog', xscale_kw=dict(linthreshx=1e-2), 
        xtickminor=True, xgrid=True, xgridcolor=c1, xtickdir='in',
        rc_kw={'grid.alpha':0.4},
    )
    
    
    '''
    2. q & Jacobianq
    '''
    axes[1].linex(baseProf.p, baseProf.q, color=c0)

    ox1 = axes[1].altx(color=c2, label=r'Jacobian Humidity $ \partial{T_b} / \partial{q} $  [K/ppmv]', lw=1)
    art_allsky = ox1.linex(baseProf.p, JacobianProf.q, color=c2, ls='-',
        label=r'$ \partial{T_b} / \partial{q} $')
    art_clearsky = ox1.linex(baseProf.p, JacobianClearProf.q, color=c2, ls='--',
        label=r'$ \partial{T_b} / \partial{q} $')
    handles_allsky.extend(art_allsky)
    handles_clearsky.extend(art_clearsky)

    axes[1].format(
        urtitle=r'\textbf{q}',
        xlabel=r'Water Vapour [ppmv]',
        xscale='log', xformatter='sci',
        yformatter='{x:.0f}hPa',
        xlim=(1e1, 1e4), xticks=(1e1, 1e2, 1e3, 1e4), 
    )

    ox1.format(
        xlim=(-1e-3, 1e-3), 
        xformatter='sci', xscale='symlog', xscale_kw=dict(linthreshx=1e-4), 
        xtickminor=True, xgrid=True, xgridcolor=c2, xtickdir='in',
        rc_kw={'grid.alpha':0.4},
    )


    '''
    3. cc & Jacobiancc
    '''
    axes[2].linex(baseProf.p, baseProf.cld_profile.cc, color=c0)

    ox2 = axes[2].altx(color=c3, label=r'Jacobian cc $ \partial{T_b} / \partial{cc} $  [K/-]', lw=1)
    ox2.linex(baseProf.p, JacobianProf.cld_profile.cc, color=c3, label=r'$ \partial{T_b} / \partial{cc} $')
    handles_allsky.extend(ox2.get_legend_handles_labels()[0])

    axes[2].format(
        urtitle=r'\textbf{CC}',
        xlabel=r'Cloud Fraction [0-1]',
        xlim=(0., 1.), xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
        yformatter='{x:.0f}hPa',
    )

    ox2.format(
        xlim=(-5.0, 5.0), xscale='symlog', xscale_kw=dict(linthreshx=1e-1),
        xtickminor=True, xgrid=True, xgridcolor=c3, xtickdir='in',
        rc_kw={'grid.alpha':0.4}, 
    )    

    '''
    4. clw/rain & Jacobian clw/rain
    '''
    axes[3].linex(baseProf.p, baseProf.cld_profile.clw, color=c0, ls='--')
    axes[3].linex(baseProf.p, baseProf.cld_profile.rain, color=c0, ls='-')

    ox3 = axes[3].altx(color=c4, 
        label=r'Jacobian CLW/RAIN $ \partial{T_b} / \partial{q_{liq}} $  [K/(kg/kg)]', lw=1)
    ox3.linex(baseProf.p, JacobianProf.cld_profile.clw, color=c4, ls='--', 
        label=r'$ \partial{T_b} / \partial{q_{CLW}} $')
    ox3.linex(baseProf.p, JacobianProf.cld_profile.rain, color=c4, ls='-', 
        label=r'$ \partial{T_b} / \partial{q_{RAIN}} $')
    handles_allsky.extend(ox3.get_legend_handles_labels()[0])

    axes[3].format(
        urtitle=r'\textbf{CLW/RAIN}',
        xlabel=r'CLW/RAIN Mxing Ratio [kg/kg]',
        xscale='log', xformatter='sci',
        yformatter='{x:.0f}hPa',
        xlim=(1e-6, 1e-3), xticks=(1e-6, 1e-5, 1e-4, 1e-3),
    )

    ox3.format(
        xlim=(-1e5, 1e5), 
        xformatter='sci', xscale='symlog', xscale_kw=dict(linthreshx=1e4),
        xtickminor=True, xgrid=True, xgridcolor=c4, xtickdir='in',
        rc_kw={'grid.alpha':0.4}, 
    )
    

    '''
    5. ciw/snow & Jacobian ciw/snow
    '''

    axes[4].linex(baseProf.p, baseProf.cld_profile.ciw, color=c0, ls='--')
    axes[4].linex(baseProf.p, baseProf.cld_profile.snow, color=c0, ls='-')

    ox4 = axes[4].altx(color=c5, 
        label=r'Jacobian CIW/SNOW $ \partial{T_b} / \partial{q_{sol}} $ [K/(kg/kg)]', lw=1)
    ox4.linex(baseProf.p, JacobianProf.cld_profile.ciw, color=c5, ls='--',
        label=r'$ \partial{T_b} / \partial{q_{CIW}} $')
    ox4.linex(baseProf.p, JacobianProf.cld_profile.snow, color=c5, ls='-',
        label=r'$ \partial{T_b} / \partial{q_{SNOW}} $')
    handles_allsky.extend(ox4.get_legend_handles_labels()[0])
    
    axes[4].format(
        urtitle=r'\textbf{CIW/SNOW}',
        xlabel=r'CIW/SNOW Mxing Ratio [kg/kg]',
        xscale='log', xformatter='sci',
        yformatter='{x:.0f}hPa',
        xlim=(1e-6, 1e-3), xticks=(1e-6, 1e-5, 1e-4, 1e-3),
    )

    ox4.format(
        xlim=(-1e5, 1e5), 
        xformatter='sci', xscale='symlog', xscale_kw=dict(linthreshx=1e4),
        xtickminor=True, xgrid=True, xgridcolor=c5, xtickdir='in',
        rc_kw={'grid.alpha':0.4},  
    )
    
    '''
    6. auxiliary axes
    '''
    with pplt.rc.context({'axes.spines.left': False, 'axes.spines.right': False,
        'axes.spines.bottom': False, 'axes.spines.left': False}):
        ax_aux = fig.add_subplot(233) # 2rows, 3cols, the 3rd axes

        ax_aux.format(
            yformatter='{x:.0f}hPa',
        )

        ax_aux.get_xaxis().set_visible(False)
        ax_aux.get_yaxis().set_visible(False)

        ax_aux.legend(title=r'\textbf{Clear-Sky}', handles=handles_clearsky, loc='ul', align='left',
            frame=True, ncols=1, order='F')
        ax_aux.legend(title=r'\textbf{All-Sky}', handles=handles_allsky, loc='ur', align='left',
            frame=True, ncols=1, order='F')
        
        InstText = r'$\textbf{'+get_channel_str('mwri', ch_no)+'}$'
        ax_aux.text(0.0, 0.2, InstText, transform='axes', size='large', ha='left', va='top')

        from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
        lon_formatter, lat_formatter = LongitudeFormatter(), LatitudeFormatter() 
        LatLonText = r'$\textbf{'+ 'Profile at {}, {}'.format( lat_formatter(baseProf.viewingConditions['latitude']), 
            lon_formatter(baseProf.viewingConditions['longitude']) ) +'}$'
        ax_aux.text(0.0, 0.05, LatLonText, transform='axes', size='large', ha='left', va='top')

    
    fig.save(pic)

if __name__ == '__main__':
    with open(JACOBIANFILE_allsky, 'r') as fin:
        baseProf = rttov_profile()
        baseProf.readfromfile(fin)
        baseCldProf = rttov_cld_profile()
        baseCldProf.readfromfile(fin)
        baseProf.cld_profile = baseCldProf

        JacobianProfs = []
        for ichannel in range(nchannels):
            JacobianProf = rttov_profile()
            JacobianProf.readfromfile(fin)
            JacobianCldProf = rttov_cld_profile()
            JacobianCldProf.readfromfile(fin)
            JacobianProf.cld_profile = JacobianCldProf
            JacobianProfs.append(JacobianProf)
    
    with open(JACOBIANFILE_clearsky, 'r') as fin:
        baseProfClear = rttov_profile()
        baseProfClear.readfromfile(fin)
        
        JacobianClearProfs = []
        for ichannel in range(nchannels):
            JacobianClearProf = rttov_profile()
            JacobianClearProf.readfromfile(fin)
            JacobianClearProfs.append(JacobianClearProf)
    
    for ch_no in PLOT_CHANNELS:
        pic = '{}Jacobian_ch{}.png'.format(directory, ch_no+1)
        print(pic)
        plot_jacobian_profile(baseProf, JacobianProfs, JacobianClearProfs, pic, ch_no)

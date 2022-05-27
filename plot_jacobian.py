'''
Description: plot jacobian of single observation
Author: Hejun Xie
Date: 2022-05-27 15:56:31
LastEditors: Hejun Xie
LastEditTime: 2022-05-28 00:47:00
'''

# Global imports
import numpy as np
import proplot as pplt

from read_jacobian import rttov_profile, rttov_cld_profile

JACOBIANFILE = './jacobian.dat'
nchannels = 10

c0 = 'k'
c1 = 'red8'
c2 = 'blue8'
c3 = 'yellow8'
c4 = 'green8'
c5 = 'cyan8'

def plot_jacobian_profile(baseProf, JacobianProf, pic):
    # T, q, cc, clw/rain, ciw/snow 
    pplt.rc.update({'meta.width': 1, 'label.weight': 'bold', 'tick.labelweight': 'bold', 
    'fontname': 'Source Sans Pro'})
    
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
        # ylabel=r'Pressure [hPa]',
        abc=r'\textbf{[A].}', abcloc='ul',
        xtickdir='in', ytickdir='in',
        ylim=(1000., 200.), yreverse=True, yscale='log', 
        yticks=(200., 300., 500., 700., 850., 1000.), ytickminor=False,
        yticklabels=(r'200hPa', r'300hPa', r'500hPa', r'700hPa', r'850hPa', r'1000hPa'),
        xgridminor=True, xticklabelpad=5, xtickminor=True,
        )
    '''
    1. T & JacobianT
    '''
    axes[0].format(
        urtitle=r'\textbf{T}',
        xlabel=r'Temperature [K]',
        xlim=[200, 300],
        xticks=(200, 220, 240., 260., 280., 300), 
    )

    axes[0].linex(baseProf.p, baseProf.T, color=c0)

    ox0 = axes[0].altx(color=c1, label=r'Jacobian Temperature [-]', lw=1)
    ox0.format(
        xlim=(-1e-1, 1e-1), 
        xformatter='sci', # xscale='symlog', xticks=(-1e-1, -5e-2, 0., 5e-2, 1e-1), 
        xtickminor=True, xgrid=True, xgridcolor=c1,
        yloc='zero',
    )
    ox0.linex(baseProf.p, JacobianProf.T, color=c1)

    
    '''
    2. q & Jacobianq
    '''
    axes[1].format(
        urtitle=r'\textbf{q}',
        xlabel=r'Water Vapour [ppmv]',
        xscale='log', xformatter='sci',
        xlim=(1e1, 1e4), xticks=(1e1, 1e2, 1e3, 1e4), 
    )

    axes[1].linex(baseProf.p, baseProf.q, color=c0)

    ox1 = axes[1].altx(color=c2, label=r'Jacobian Humidity [K/ppmv]', lw=1)
    ox1.format(
        xlim=(-1e-3, 1e-3), 
        xformatter='sci', # xscale='symlog', xticks=(-1e-1, -5e-2, 0., 5e-2, 1e-1)
        xtickminor=True, xgrid=True, xgridcolor=c2,
        yloc='zero',
    )

    ox1.linex(baseProf.p, JacobianProf.q, color=c2)

    '''
    3. cc & Jacobiancc
    '''
    axes[2].format(
        urtitle=r'\textbf{CC}',
        xlabel=r'Cloud Fraction [0-1]',
        xlim=(0., 1.), xticks=(0.0, 0.2, 0.4, 0.6, 0.8, 1.0),
    )

    axes[2].linex(baseProf.p, baseProf.cld_profile.cc, color=c0)

    ox2 = axes[2].altx(color=c3, label=r'Jacobian cc [K/-]', lw=1)
    ox2.format(
        xlim=(-5.0, 5.0), xscale='symlog', xscale_kw=dict(linthreshx=1e-1),
        # xscale='symlog', xticks=(-1e-1, -5e-2, 0., 5e-2, 1e-1)
        xtickminor=True, xgrid=True, xgridcolor=c3,
        yloc='zero',
    )
    ox2.linex(baseProf.p, JacobianProf.cld_profile.cc, color=c3)

    '''
    4. clw/rain & Jacobian clw/rain
    '''
    axes[3].format(
        urtitle=r'\textbf{CLW/RAIN}',
        xlabel=r'CLW/RAIN Mxing Ratio [kg/kg]',
        xscale='log', xformatter='sci',
        xlim=(1e-6, 1e-3), xticks=(1e-6, 1e-5, 1e-4, 1e-3),
    )

    axes[3].linex(baseProf.p, baseProf.cld_profile.clw, color=c0, ls='--')
    axes[3].linex(baseProf.p, baseProf.cld_profile.rain, color=c0, ls='-')

    ox3 = axes[3].altx(color=c4, label=r'Jacobian CLW/RAIN [K/(kg/kg)]', lw=1)
    ox3.format(
        xlim=(-1e5, 1e5), 
        xformatter='sci', xscale='symlog', xscale_kw=dict(linthreshx=1e4),
        xtickminor=True, xgrid=True, xgridcolor=c4,
        yloc='zero',
    )
    ox3.linex(baseProf.p, JacobianProf.cld_profile.clw, color=c4, ls='--')
    ox3.linex(baseProf.p, JacobianProf.cld_profile.rain, color=c4, ls='-')

    '''
    5. ciw/snow & Jacobian ciw/snow
    '''
    
    axes[4].format(
        urtitle=r'\textbf{CIW/SNOW}',
        xlabel=r'CIW/SNOW Mxing Ratio [kg/kg]',
        xscale='log', xformatter='sci',
        xlim=(1e-6, 1e-3), xticks=(1e-6, 1e-5, 1e-4, 1e-3),
    )

    axes[4].linex(baseProf.p, baseProf.cld_profile.ciw, color=c0, ls='--')
    axes[4].linex(baseProf.p, baseProf.cld_profile.snow, color=c0, ls='-')

    ox4 = axes[4].altx(color=c5, label=r'Jacobian CIW/SNOW [K/(kg/kg)]', lw=1)
    ox4.format(
        xlim=(-1e5, 1e5), 
        xformatter='sci', xscale='symlog', xscale_kw=dict(linthreshx=1e4),
        xtickminor=True, xgrid=True, xgridcolor=c5,
        yloc='zero',
    )
    ox4.linex(baseProf.p, JacobianProf.cld_profile.ciw, color=c5, ls='--')
    ox4.linex(baseProf.p, JacobianProf.cld_profile.snow, color=c5, ls='-')
    
    fig.save(pic)

if __name__ == '__main__':
    with open(JACOBIANFILE, 'r') as fin:
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
    
    for ichannel in range(nchannels):
        ichannel = 4
        pic = 'Jacobian_ch{}.png'.format(ichannel+1)
        plot_jacobian_profile(baseProf, JacobianProfs[ichannel], pic) # Channel 5 23.8GHz V pol.
        exit()
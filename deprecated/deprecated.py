'''
Description: save deprecated code here
Author: Hejun Xie
Date: 2022-05-19 21:33:48
LastEditors: Hejun Xie
LastEditTime: 2022-05-19 21:35:37
'''


def draw_worldmap(ax):
    from mpl_toolkits.basemap import Basemap
    
    map = Basemap(projection='cyl', llcrnrlat=slat, urcrnrlat=elat, llcrnrlon=slon, urcrnrlon=elon, resolution='i', ax=ax)
    map.drawparallels(np.arange(slat, elat+5, 30.0), linewidth=1, dashes=[4, 3], labels=[1, 0, 0, 0])
    map.drawmeridians(np.arange(slon, elon+5, 30.0), linewidth=1, dashes=[4, 3], labels=[0, 0, 0, 1])
    map.drawcoastlines()
    return map


def plot_inno_scatter(map, fig, ax, asi, ch_no):
    '''
    plot innovation scatter
    '''
    for itimeslot,timeslot in enumerate(asi.timeslots):
        data = asi.domain_merged_data[timeslot][get_channel_key(ch_no)]
        try:
            lonTag, latTag = data['lon'], data['lat']
            OMB = (- data['BMO(BC)']) / data['cldraderr']
        except KeyError:
            print('Empty timeslot {}'.format(timeslot))
            continue
        
        x, y = map(lonTag, latTag)
        im = map.scatter(x, y, marker='o', edgecolor='k', s=8, c=OMB, linewidths=itimeslot/24+0.2,cmap='rainbow', vmin = -2.0, vmax = 2.0)

    position = fig.add_axes([0.972,0.245,0.012,0.68])
    cb  = fig.colorbar(im, cax=position, shrink=0.80, pad=0.01)
    cb.ax.set_ylabel('Normalized FG departure [K]')


def plot_hydro(map, ax, TLON, TLAT, hydro):
    '''
    PLOT HYDRO
    '''
    plot_hydro = np.ma.masked_array(hydro, hydro<Q_LEVELS[0])
    x, y = map(TLON.T, TLAT.T)
    norm = colors.LogNorm(vmin=Q_LEVELS[0], vmax=Q_LEVELS[-1])
    im = map.pcolormesh(x, y, plot_hydro.T, cmap='Greys', norm=norm, shading='auto', alpha=0.5)

def plot_level(ds_dxa, ds_dxa_plevel, ds_xb, ds_xb_plevel, ds_grapesinput, asi, ch_no, level='bottom'):
    '''
    plot allsky innovations
    '''

    TLON, TLAT = np.meshgrid(ds_dxa.coords['lon'], ds_dxa.coords['lat'])
    TLON_GI, TLAT_GI = np.meshgrid(ds_grapesinput.coords['lon'], ds_grapesinput.coords['lat'])

    if level == 'bottom':
        dRH = ds_dxa.data_vars['Increment Relative Humidity'].data[0,...]
        dT = ds_dxa.data_vars['Increment Temperature'].data[0,...]
    else:
        dRH = ds_dxa_plevel.data_vars['Increment Relative Humidity'].sel(Pressure=level, method='nearest')
        dT= ds_dxa_plevel.data_vars['Increment Temperature'].sel(Pressure=level, method='nearest')

    qr_column = np.max(ds_grapesinput.data_vars['QR'].data, axis=0)

    figsize=(17,7.5)

    '''
    1. RH
    '''
    fig, ax = plt.subplots(figsize=figsize)
    map = draw_worldmap(ax)

    x, y = map(TLON.T, TLAT.T)

    dRHmin, dRHmax = get_vminvmax(dRH.data.flatten(), vstage=2.0, ign_perc=0.5, lsymmetric=True)
    print('dRHmin={:>.1f}, dRHmax={:>.1f}'.format(dRHmin, dRHmax))
    RH_LEVELS_CF = np.linspace(dRHmin, dRHmax, 21)
    RH_LEVELS_CS = np.linspace(dRHmin, dRHmax, 11)
    RH_TICKS = np.linspace(dRHmin, dRHmax, 11)

    norm = colors.Normalize(vmin=RH_TICKS[0], vmax=RH_TICKS[-1])

    CF = map.contourf(x, y, dRH.T, RH_LEVELS_CF, cmap='BrBG', norm=norm, extend='both')    
    CS = map.contour(x, y, dRH.T, RH_LEVELS_CS, colors=('k',), linewidths=(1.,))    
    ax.clabel(CS, fmt='%2.1f', colors='k', fontsize=10)

    CB = fig.colorbar(CF, ax=ax, shrink=0.80, pad=0.08, aspect=60, orientation='horizontal')
    CB.ax.set_xlabel(r'Analysis Increment Relative Humidity [\%]', fontsize=25)
    CB.set_ticks(RH_TICKS)

    # plot_hydro(map, ax, TLON_GI, TLAT_GI, qr_column)
    # plot_inno_scatter(map, fig, ax, asi, ch_no)

    ax.text(0.0, 1.0, r'$\textbf{'+'{}hPa'.format(level)+'}$', 
        ha='left', va='bottom', size=24, transform=ax.transAxes)
    ax.text(1.0, 1.0, r'$\textbf{'+get_channel_str(asi.instrument, ch_no-1)+' [{}-{}]'.format(asi.timeslots[0], asi.timeslots[-1])+'}$', 
        ha='right', va='bottom', size=24, transform=ax.transAxes)
    plt.savefig('dRH_level-{}_ch{}.png'.format(level, ch_no), bbox_inches='tight')
    plt.close()

    exit()

    '''
    2. T
    '''
    fig, ax = plt.subplots(figsize=figsize)
    map = draw_worldmap(ax)   
    x, y = map(TLON.T, TLAT.T)

    dTmin, dTmax = get_vminvmax(dT.data.flatten(), vstage=0.2, ign_perc=0.0, lsymmetric=True)
    print('dTmin={:>.1f}, dTmax={:>.1f}'.format(dTmin, dTmax))
    T_LEVELS_CF = np.linspace(dTmin, dTmax, 21)
    T_LEVELS_CS = np.linspace(dTmin, dTmax, 11)
    T_TICKS = np.linspace(dTmin, dTmax, 11)

    norm = colors.Normalize(vmin=T_TICKS[0], vmax=T_TICKS[-1])

    CF = map.contourf(x, y, dT.T, T_LEVELS_CF, cmap='RdBu_r', norm=norm, extend='both')    
    CS = map.contour(x, y, dT.T, T_LEVELS_CS, colors=('k',), linewidths=(1.,))
        
    ax.clabel(CS, fmt='%2.1f', colors='k', fontsize=10)

    CB = fig.colorbar(CF, ax=ax, shrink=0.80, pad=0.08, aspect=60, orientation='horizontal')
    CB.ax.set_xlabel(r'Analysis Increment Temperature [K]', fontsize=25)
    CB.set_ticks(T_TICKS)

    plot_hydro(map, ax, TLON_GI, TLAT_GI, qr_column)
    plot_inno_scatter(map, fig, ax, asi, ch_no)

    ax.text(0.0, 1.0, r'$\textbf{'+'{}hPa'.format(level)+'}$', 
        ha='left', va='bottom', size=24, transform=ax.transAxes)
    ax.text(1.0, 1.0, r'$\textbf{'+get_channel_str(asi.instrument, ch_no-1)+' [{}-{}]'.format(asi.timeslots[0], asi.timeslots[-1])+'}$', 
        ha='right', va='bottom', size=24, transform=ax.transAxes)
    plt.savefig('dT_level-{}_ch{}.png'.format(level, ch_no), bbox_inches='tight')
    plt.close()

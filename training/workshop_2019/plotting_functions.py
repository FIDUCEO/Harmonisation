import numpy as np
import cartopy
import cartopy.crs as ccrs
from cartopy.io import shapereader
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
import matplotlib.pylab as plt; plt.close('all')
import matplotlib.pyplot as plt; plt.close("all")
import matplotlib.dates as mdates
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy.ma as ma



def plot_lisbon_var(var, lon, lat, label, cmap='viridis', projection='platecarree'):

    x = lon[::4,::4]
    y = lat[::4,::4]
    z = var[::4,::4]

    fig  = plt.figure()
    p = ccrs.PlateCarree(central_longitude=0)
    threshold = 0

    ax = plt.axes(projection=p)
    ax.stock_img()
    ax.coastlines(resolution='50m')

    rivers_50m = cfeature.NaturalEarthFeature('physical', 'rivers_lake_centerlines', '50m')    
    ax.add_feature(rivers_50m, facecolor='None', edgecolor='b', linewidth=0.5)
    ax.add_feature(cfeature.BORDERS.with_scale('50m'))
    ax.add_feature(cfeature.LAKES.with_scale('50m'))
    
    g = ccrs.Geodetic()
    trans = ax.projection.transform_points(g, x, y)
    x0 = trans[:,:,0]
    x1 = trans[:,:,1]

    ax.set_extent([-12,-1,35,45])
    gl = ax.gridlines(crs=p, draw_labels=True, linewidth=1, color='gray', alpha=0.5, linestyle='-')
    gl.xlabels_top = False
    gl.ylabels_right = False
    gl.xlines = True
    gl.ylines = True
    gl.xlocator = mticker.FixedLocator([-180,-120,-60,0,60,120,180])
    gl.ylocator = mticker.FixedLocator([-90,-60,-30,0,30,60,90])
    gl.xformatter = LONGITUDE_FORMATTER
    gl.yformatter = LATITUDE_FORMATTER

    im = ax.pcolor(x, y, z, transform=ax.projection, cmap='gnuplot2')

    im.set_clim(np.nanmin(var),np.nanmax(var))

    # Lisbon: 38.736946, -9.142685
    ax.plot(-9.142685, 38.736946, 'ro', markersize=3, transform=p)
    ax.text(-10, 38, 'Lisbon', transform=p)

    cb = plt.colorbar(im, orientation="horizontal", extend='both', label=label)
    plt.show()
    plt.close('all')

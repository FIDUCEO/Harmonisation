""" Test code for reading netCDF files with xarray package. """

import xarray as xr
import numpy as np
import pandas as pd

# Folder where the data is
dataDir = "D:\Projects\FIDUCEO\Data\Simulated"

# Name of the netCDF file to read
dF = "matchup_sim_srf.n16.nc"

temp = 15 + 8 * np.random.randn(2, 2, 3)
precip = 10 * np.random.rand(2, 2, 3)
lon = [[-99.83, -99.32], [-99.79, -99.23]]
lat = [[42.25, 42.21], [42.63, 42.59]]

# for real use cases, its good practice to supply array attributes such as
# units, but we won't bother here for the sake of brevity
ds = xr.Dataset({'temperature': (['x', 'y', 'time'],  temp),
                'precipitation': (['x', 'y', 'time'], precip)},
                coords={'lon': (['x', 'y'], lon),
                'lat': (['x', 'y'], lat),
                'time': pd.date_range('2014-09-06', periods=3),
                'reference_time': pd.Timestamp('2014-09-05')})
ds

# Read netCDF file to an xarray.Dataset
# data = xr.open_dataset(dF, decode_cf=False)
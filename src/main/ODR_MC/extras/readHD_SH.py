""" Functions to read netCDF file of a sensor pair.
Use netCDF4 package to read netCDF data, 
and pandas with data frame for storage and use of netCDF data. """


'''---Python Modules---'''
import numpy as np
from netCDF4 import Dataset
import pandas as pd
from datetime import datetime as dt


class HData:
    def __init__(self, inputPath):
        """Return harmonisation data for input path"""

        # Open data for:
        # H - Harmonisation matchup data
        # Ur - Random matchup uncertainties
        # Us - Systematic matchup uncertainties
        # timedata - Matchup times
        # sensor1 - Sensor 1 ID
        # sensor2 - Sensor 2 ID
        # n_mu - Number of matchups

        self.H, self.Ur, self.Us, self.timedata, \
            self.sensor1, self.sensor2, self.n_mu = self.openData_pd(inputPath)

    def openData_pd(self, inputPath):
        """Return harmonisation data from input NetCDF file as pandas data frame"""

        print('Opening: ', inputPath)

        # Open matchup NetCDF file
        rootgrp = Dataset(inputPath, 'r')

        # Open auxiliary data
        timedata = self.conv2date(rootgrp.variables['CorrIndexArray'][:])
        sensor1 = rootgrp.variables['lm'][0, 0]   # sensor 1 number
        sensor2 = rootgrp.variables['lm'][0, 1]   # sensor 2 number
        n_mu = rootgrp.variables['lm'][0, 2]      # number of match ups

        #  Open harmonisation data as pandas data frame

        # Select data required to open
        rspair = False
        if sensor1 == -1:
            rspair = True

        didx, colID = self.select_variables(rspair)

        # 1. Harmonisation data
        H = pd.DataFrame(data=rootgrp.variables['H'][:, didx], columns=colID)
        H['sAdj'] = rootgrp.variables['K'][:]                # add K values

        # 2. Random  uncertainties
        Ur = pd.DataFrame(data=rootgrp.variables['Ur'][:, didx], columns=colID)
        Ur['sAdj'] = rootgrp.variables['Kr'][:]    # add K random uncertainty

        # 3. Systematic Uncertainties
        Us = pd.DataFrame(data=rootgrp.variables['Us'][:, didx], columns=colID)
        Us['sAdj'] = rootgrp.variables['Ks'][:]    # add K random uncertainty

        # Close NetCDF file
        rootgrp.close()

        return H, Ur, Us, timedata, sensor1, sensor2, n_mu

    def opendata_xarray(self, inputPath):
        """Return harmonisation data from input NetCDF file as xarray"""

        # TODO

        return 0

    def select_variables(self, rspair):
        """Return index of NetCDF variables required and column ids for data frame"""

        # Select data variables required:
        #  > RS Pair - for reference-sensor pair only Lref is needed
        #  > SS pair - for sensor-sensor pair all data required
        if rspair:
            didx = [0, 5, 6, 7, 8, 9]
            colID =['Lref', 'Cspace2', 'CICT2', 'CEarth2', 'LICT2', 'Torb2']
        else:
            didx = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
            colID = ['Lref', 'Cspace', 'CICT', 'CEarth', 'LICT', 'Torb'
                     'Cspace2', 'CICT2', 'CEarth2', 'LICT2', 'Torb2']

        return didx, colID

    def conv2date(self, inTime):
        """Return unix time in seconds (from 1970) from AVHRR time in seconds (from 1975)"""

        # Calculate difference from AVHRR start time and unix start time in seconds
        start_time_AVHRR = dt(1975, 1, 1)
        start_time_unix = dt(1970, 1, 1)
        time_diff = (start_time_AVHRR - start_time_unix).total_seconds()

        # Convert to time from 1975 to date
        outTime = [dt.fromtimestamp(time+time_diff) for time in inTime]

        return outTime

if __name__ == '__main__':
    def main():
        inputPath = r"C:\Users\seh2\data\Harmonisation\m02_n15.nc"

        hdata = HData(inputPath)

        return 0

    main()


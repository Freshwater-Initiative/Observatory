
import pandas as pd

import xarray as xr

import os as os

homedir = 'D:/UW_PhD/Python/xarray/'
os.chdir(homedir)

#open pnnl lat long el file
pnnlxy = xr.open_dataset('data_LatLonGht.nc')

#range of dates to create new files for
start_date = '2006-11-01'
end_date = '2006-11-20'

dates = [x.strftime('%Y-%m-%d') for x in pd.date_range(start=start_date, end=end_date, freq='D')]
#for each day (file) open pnnl netcdf file as xarray dataset, add coordinate data, save as new netcdf
for ind, ymd in enumerate(dates):
    date = ymd
    pnnl = xr.open_dataset('data.' + date +'.nc')


    #merge pnnl lat long data with climate data
    pnnlnew = xr.merge([pnnl, pnnlxy], compat='no_conflicts')

    #convert variables LAT, LON and Z to coordinates
    pnnlnewc = pnnlnew.set_coords({'LAT','LON','Z'}, inplace=False)

    #create series of dates to add to dataset
    time = pd.date_range(start=date, periods=24, freq='H')

    #add coordinates using series of dates
    pnnlnewc.update({'time': ('time', time)})


    #problem with xarray dataframe.to_netcdf method: does not work if attributes contain attribute titled "coordinates"; 
    #for variables that have coordinates in attributes, remove "coordinate" from attributes
    pnnlnewc.Q2.attrs = ([('stagger', ''),
                          ('units', 'kg kg-1'),
                          ('description', 'QV at 2 M'),
                          ('MemoryOrder', 'XY '),
                          ('FieldType', 104)])
    
    pnnlnewc.PSFC.attrs = ([('stagger', ''),
                            ('units', 'Pa'),
                            ('description', 'SFC PRESSURE'),
                            ('MemoryOrder', 'XY '),
                            ('FieldType', 104)])
    
    pnnlnewc.GLW.attrs = ([('stagger', ''),
                           ('units', 'W m-2'),
                           ('description', 'DOWNWARD LONG WAVE FLUX AT GROUND SURFACE'),
                           ('MemoryOrder', 'XY '),
                           ('FieldType', 104)])
    
    pnnlnewc.SWDOWN.attrs = ([('stagger', ''),
                              ('units', 'W m-2'),
                              ('description', 'DOWNWARD SHORT WAVE FLUX AT GROUND SURFACE'),
                              ('MemoryOrder', 'XY '),
                              ('FieldType', 104)])
    
    pnnlnewc.PREC_ACC_NC.attrs = ([('stagger', ''),
                                   ('units', 'mm'),
                                   ('description', 'GRID SCALE  PRECIPITATION'),
                                   ('MemoryOrder', 'XY '),
                                   ('FieldType', 104)])
    
    pnnlnewc.SNOW_ACC_NC.attrs = ([('stagger', ''),
                                   ('units', 'mm'),
                                   ('description', 'SNOW WATER EQUIVALENT'),
                                   ('MemoryOrder', 'XY '),
                                   ('FieldType', 104)])
    
    #save new netcdf file
    pnnlnewc.to_netcdf(date + '_coord.nc')
    print(str(ind) + ' : ' + date)
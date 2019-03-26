# import libraries
import os
import numpy as np
import pandas as pd
import dask as da
import xarray as xray
import ftplib, wget, bz2, urllib
from multiprocessing.pool import ThreadPool

from dask.diagnostics import ProgressBar
import ogh


def compile_x_wrfnnrp_raw_Salathe2014_locations(time_increments):
    """
    Compile a list of file URLs for Salathe et al., 2014 raw WRFdata
    
    time_increments: (list) a list of dates that identify each netcdf file
    """
    locations=[]
    domain='http://cses.washington.edu'
    subdomain='/rocinante/WRF/NNRP/vic_16d/PNW_1970_1999/WRF_NNRP_noBC/netcdf_daily'
WRF_NNRP_noBC.197001.nc	
    for ind, yearmo in enumerate(time_increments):
        basename='/WRF_NNRP_noBC.{0}.nc'.format(yearmo)
        url='{0}{1}{2}'.format(domain, subdomain, basename)
        locations.append(url)
    return(locations)

def compile_x_wrfpnnl2018_raw_locations(time_increments):
    """
    Compile a list of file URLs for PNNL 2018 raw WRF NNRP data   
    
    time_increments: (list) a list of dates that identify each netcdf file
    """
    locations=[]
    domain='http://cses.washington.edu'
    subdomain='/rocinante/WRF/PNNL_NARR_6km'

    for ind, ymd in enumerate(time_increments):
        basename='/' + ymd[0:4] + '/data.' + ymd[0:4] + '-' + ymd[4:6] + '-' + ymd[6:8] + '.nc'
        url='{0}{1}{2}'.format(domain, subdomain, basename)
        locations.append(url)
    return(locations)

def wget_x_download_spSubset(fileurl, 
                             spatialbounds, 
                             file_prefix='sp_', 
                             rename_latlong_names={'LAT':'LAT','LON':'LON'}, 
                             replace_file=True):
    """
    Download files from an http domain
    
    fileurl: (str) a urls to request a netcdf file
    spatialbounds: (dict) dictionary providing the minx, miny, maxx, and maxy of the spatial region
    file_prefix: (str) a string to mark the output file as a spatial subset
    rename_latlong_names: (dict) a dictionary to standardize latitude/longitude synonyms to LAT/LON, respectively
    replace_file: (logic) If True, the existing file will be replaced; if False, the file download is skipped
    """
    
    # check if the file path already exists; if so, apply replace_file logic
    basename = os.path.basename(fileurl)
    if os.path.isfile(basename):
        os.remove(basename)
        
    if os.path.isfile(file_prefix+basename) and replace_file:
        os.remove(file_prefix+basename)
    elif os.path.isfile(file_prefix+basename) and not replace_file:
        # replace_file is False; return file path and skip
        return(os.path.join(os.getcwd(), file_prefix+basename))
        
    # try the file connection
    #print('connecting to: '+basename)
    try:
        ping = urllib.request.urlopen(fileurl)
        
        # if the file exists, download it
        if ping.getcode() != 404:
            ping.close()
            wget.download(fileurl)

            # open the parent netcdf file
            ds = xray.open_dataset(basename, engine = 'netcdf4')
            
            # rename latlong if they are not LAT and LON, respectively
            if not isinstance(rename_latlong_names, type(None)):
                ds = ds.rename(rename_latlong_names)
                
            # slice by the bounding box
            spSubset = ds.sel(LON=slice(spatialbounds['minx'], spatialbounds['maxx']),
                              LAT=slice(spatialbounds['miny'], spatialbounds['maxy']))

            # print the spatial subset
            spSubset.to_netcdf(file_prefix+basename)
            #print('downloaded: spatial subset of ' + basename)
            
            # remove the parent
            ds.close()
            os.remove(basename)
            return(os.path.join(os.getcwd(), file_prefix+basename))

        else:
            ping.close()
    except:
        print('File does not exist at this URL: ' + basename)

        
def get_x_dailywrf_Salathe2014(homedir,
                               spatialbounds,
                               subdir='salathe2014/Daily_WRF_1970_1999/noBC',
                               nworkers=4,
                               start_date='1970-01-01',
                               end_date='1989-12-31',
                               rename_timelatlong_names={'LAT':'LAT','LON':'LON'},
                               file_prefix='sp_',
                               replace_file=True):
    """
    get Daily WRF data from Salathe et al. (2014) using xarray on netcdf files
    """
    # check and generate DailyMET livneh 2013 data directory
    filedir=os.path.join(homedir, subdir)
    ogh.ensure_dir(filedir)
    
    # modify each month between start_date and end_date to year-month
    dates = [x.strftime('%Y%m') for x in pd.date_range(start=start_date, end=end_date, freq='M')]

    # initialize parallel workers
    da.set_options(pool=ThreadPool(nworkers))
    ProgressBar().register()
    
    # generate the list of files to download
    filelist = compile_x_wrfnnrp_raw_Salathe2014_locations(dates)

    # download files of interest
    NetCDFs=[]
    for url in filelist:
        NetCDFs.append(da.delayed(wget_x_download_spSubset)(fileurl=url,
                                                            spatialbounds=spatialbounds,
                                                            file_prefix=file_prefix,
                                                            rename_latlong_names=rename_timelatlong_names,
                                                            replace_file=replace_file))

    # run operations
    outputfiles = da.compute(NetCDFs)[0]

    # reset working directory
    os.chdir(homedir)
    return(outputfiles)

def get_x_hourlywrf_PNNL2018(homedir,
                               spatialbounds,
                               subdir='PNNL2018/Hourly_WRF_1981_2015/noBC',
                               nworkers=4,
                               start_date='1980-10-01',
                               end_date='2015-12-31',
                               rename_timelatlong_names={'LAT':'LAT','LON':'LON'},
                               file_prefix='sp_',
                               replace_file=True):
    """
    get hourly WRF data from a 2018 PNNL WRF run using xarray on netcdf files
    """
    # check and generate data directory
    filedir=os.path.join(homedir, subdir)
    ogh.ensure_dir(filedir)
    
    # modify each month between start_date and end_date to year-month
    dates = [x.strftime('%Y%m%d') for x in pd.date_range(start=start_date, end=end_date, freq='D')]
    
    # initialize parallel workers
    da.set_options(pool=ThreadPool(nworkers))
    ProgressBar().register()
    
    # generate the list of files to download
    filelist = compile_x_wrfpnnl2018_raw_locations(dates)
    
    # download files of interest
    NetCDFs=[]
    for url in filelist:
        NetCDFs.append(da.delayed(wget_x_download_spSubset)(fileurl=url,
                                                            spatialbounds=spatialbounds,
                                                            file_prefix=file_prefix,
                                                            rename_latlong_names=rename_timelatlong_names,
                                                            replace_file=replace_file))
    
    # run operations
    outputfiles = da.compute(NetCDFs)[0]
    
    # reset working directory
    os.chdir(homedir)
    return(outputfiles)

def netcdf_to_ascii(homedir, subdir, netcdfs, mappingfile, catalog_label, meta_file):
    # initialize list of dataframe outputs
    outfiledict = {}
    
    # generate destination folder
    filedir=os.path.join(homedir, subdir)
    ogh.ensure_dir(filedir)

    # connect with collection of netcdfs
    ds_mf = xray.open_mfdataset(netcdfs, engine = 'netcdf4')

    # generate list of variables
    ds_vars = [ds_var for ds_var in dict(ds_mf.variables).keys() 
               if ds_var not in ['YEAR','MONTH','DAY','TIME','LAT','LON']]

    # convert netcdfs to pandas.Panel API
    ds_pan = ds_mf.to_dataframe()[ds_vars]

    # read in gridded cells of interest
    maptable, nstation = ogh.mappingfileToDF(mappingfile, colvar=None)

    # at each latlong of interest
    for ind, eachrow in maptable.iterrows():

        # generate ASCII time-series
        ds_df = ds_pan.loc[eachrow['LAT'],eachrow['LONG_'],:].reset_index(drop=True, level=[0,1])

        # create file name
        outfilename = os.path.join(filedir,'data_{0}_{1}'.format(eachrow['LAT'],eachrow['LONG_']))

        # save ds_df
        outfiledict[outfilename] = da.delayed(ds_df.to_csv)(path_or_buf=outfilename, sep='\t', header=False, index=False)

    # compute ASCII time-series files
    ProgressBar().register()
    outfiledict = da.compute(outfiledict)[0]
    
    # update metadata file
    meta_file[catalog_label] = dict(ds_mf.attrs)
    meta_file[catalog_label]['variable_list']=np.array(ds_vars)
    
    # catalog the output files
    ogh.addCatalogToMap(outfilepath=mappingfile, maptable=maptable, folderpath=filedir, catalog_label=catalog_label)
    os.chdir(homedir)
    return(list(outputfiles.keys()))

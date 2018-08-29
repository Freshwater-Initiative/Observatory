
# -*- coding:  utf-8 -*-

#import os
#import json


def meta_file(): 
    tmp = {}
    
    """
    Daily Meteorology -  Livneh et al.,  2013
    """
    tmp['dailymet_livneh2013']={}
    tmp['dailymet_livneh2013']['spatial_resolution'] = '1/16-degree'
    tmp['dailymet_livneh2013']['web_protocol'] = 'ftp'
    tmp['dailymet_livneh2013']['domain'] = 'livnehpublicstorage.colorado.edu'
    tmp['dailymet_livneh2013']['subdomain'] = '/public/Livneh.2013.CONUS.Dataset/Meteorology.asc.v.1.2.1915.2011.bz2/'
    tmp['dailymet_livneh2013']['decision_steps'] = 'files organized by spatial bounding boxes'
    tmp['dailymet_livneh2013']['filename_structure'] = 'Meteorology_Livneh_CONUSExt_v.1.2_2013_{LAT}_{LONG}'
    tmp['dailymet_livneh2013']['file_format'] = 'bz2-compressed ASCII'
    
    tmp['dailymet_livneh2013']['reference']={}
    tmp['dailymet_livneh2013']['reference'][1] = 'Livneh,  B.,  E. A. Rosenberg,  C. Lin,  B. Nijssen,  V. Mishra,  K. M. Andreadis,  E. P. Maurer,  and D. P. Lettenmaier,  2013:  A Long-Term Hydrologically Based Dataset of Land Surface Fluxes and States for the Conterminous United States:  Update and Extensions. J. Climate,  26,  9384-9392.'
    tmp['dailymet_livneh2013']['reference'][2]='ftp://livnehpublicstorage.colorado.edu/public/Livneh.2013.CONUS.Dataset/readme.txt'
    
    tmp['dailymet_livneh2013']['start_date'] = '1915-01-01'
    tmp['dailymet_livneh2013']['end_date'] = '2011-12-31'
    tmp['dailymet_livneh2013']['temporal_resolution'] = 'D'
    tmp['dailymet_livneh2013']['delimiter'] = '\t'
    tmp['dailymet_livneh2013']['variable_list'] = ['PRECIP',  'TMAX',  'TMIN',  'WINDSPD']
    
    tmp['dailymet_livneh2013']['variable_info']={}
    tmp['dailymet_livneh2013']['variable_info']['PRECIP']={}
    tmp['dailymet_livneh2013']['variable_info']['PRECIP']={'desc': 'daily precipitation (mm)', 
                                                           'dtypes': 'float64', 'units': 'mm'}
    tmp['dailymet_livneh2013']['variable_info']['TMAX']={}
    tmp['dailymet_livneh2013']['variable_info']['TMAX']={'desc': 'daily maximum temperature (C)', 
                                                         'dtypes': 'float64', 'units': 'C'}
    tmp['dailymet_livneh2013']['variable_info']['TMIN']={}
    tmp['dailymet_livneh2013']['variable_info']['TMIN']={'desc': 'daily minimum temperature (C)', 
                                                         'dtypes': 'float64', 'units': 'C'}
    tmp['dailymet_livneh2013']['variable_info']['WINDSPD']={}
    tmp['dailymet_livneh2013']['variable_info']['WINDSPD']={'desc': 'daily mean wind speed (m/s)', 
                                                            'dtypes': 'float64', 'units': 'm/s'}
    
    """
    Daily Meteorology - Livneh et al.,  2015
    """
    tmp['dailymet_livneh2015']={}
    tmp['dailymet_livneh2015']['spatial_resolution'] = '1/16-degree'
    tmp['dailymet_livneh2015']['web_protocol'] = 'ftp'
    tmp['dailymet_livneh2015']['domain'] = '192.12.137.7'
    tmp['dailymet_livneh2015']['subdomain'] = '/pub/dcp/archive/OBS/livneh2014.1_16deg/ascii/daily/'
    tmp['dailymet_livneh2015']['decision_steps'] = 'files organized by Latitude'
    tmp['dailymet_livneh2015']['filename_structure'] = 'Meteorology_Livneh_NAmerExt_15Oct2014_{LAT}_{LONG}'
    tmp['dailymet_livneh2015']['file_format'] = 'bz2-compressed ASCII'
    
    tmp['dailymet_livneh2015']['reference']={}
    tmp['dailymet_livneh2015']['reference'][1]= 'Livneh B.,  T.J. Bohn,  D.S. Pierce,  F. Munoz-Ariola,  B. Nijssen,  R. Vose,  D. Cayan,  and L.D. Brekke,  2015:  A spatially comprehensive,  hydrometeorological data set for Mexico,  the U.S.,  and southern Canada 1950-2013,  Nature Scientific Data,  5: 150042,  doi: 10.1038/sdata.2015.42.'
    tmp['dailymet_livneh2015']['reference'][2]='ftp: //livnehpublicstorage.colorado.edu/public/Livneh.2013.CONUS.Dataset/readme.txt'
    
    tmp['dailymet_livneh2015']['start_date'] = '1950-01-01'
    tmp['dailymet_livneh2015']['end_date'] = '2013-12-31'
    tmp['dailymet_livneh2015']['temporal_resolution'] = 'D'
    tmp['dailymet_livneh2015']['delimiter'] = '\\s+'
    tmp['dailymet_livneh2015']['variable_list'] = ['PRECIP',  'TMAX',  'TMIN',  'WINDSPD']
    
    tmp['dailymet_livneh2015']['variable_info']={}
    tmp['dailymet_livneh2015']['variable_info']['PRECIP']={}
    tmp['dailymet_livneh2015']['variable_info']['PRECIP']={'desc': 'daily precipitation (mm)', 
                                                           'dtypes': 'float64', 'units': 'mm'}
    tmp['dailymet_livneh2015']['variable_info']['TMAX']={}
    tmp['dailymet_livneh2015']['variable_info']['TMAX']={'desc': 'daily maximum temperature (C)', 
                                                         'dtypes': 'float64', 'units': 'C'}
    tmp['dailymet_livneh2015']['variable_info']['TMIN']={}
    tmp['dailymet_livneh2015']['variable_info']['TMIN']={'desc': 'daily minimum temperature (C)', 
                                                         'dtypes': 'float64', 'units': 'C'}
    tmp['dailymet_livneh2015']['variable_info']['WINDSPD']={}
    tmp['dailymet_livneh2015']['variable_info']['WINDSPD']={'desc': 'daily mean wind speed (m/s)', 
                                                            'dtypes': 'float64', 'units': 'm/s'}
    
    """
    Daily VIC - Livneh et al.,  2013
    """
    tmp['dailyvic_livneh2013']={}
    tmp['dailyvic_livneh2013']['spatial_resolution'] = '1/16-degree'
    tmp['dailyvic_livneh2013']['web_protocol'] = 'ftp'
    tmp['dailyvic_livneh2013']['domain'] = 'livnehpublicstorage.colorado.edu'
    tmp['dailyvic_livneh2013']['subdomain'] = '/public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2/'
    tmp['dailyvic_livneh2013']['decision_steps'] = 'files organized by spatial bounding boxes'
    tmp['dailyvic_livneh2013']['filename_structure'] = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_{LAT}_{LONG}'
    tmp['dailyvic_livneh2013']['file_format'] = 'bz2-compressed ASCII'
    
    tmp['dailyvic_livneh2013']['reference']={}
    tmp['dailyvic_livneh2013']['reference'][1] = 'Livneh,  B.,  E. A. Rosenberg,  C. Lin,  B. Nijssen,  V. Mishra,  K. M. Andreadis,  E. P. Maurer,  and D. P. Lettenmaier,  2013:  A Long-Term Hydrologically Based Dataset of Land Surface Fluxes and States for the Conterminous United States:  Update and Extensions. J. Climate,  26,  9384-9392.'
    tmp['dailyvic_livneh2013']['reference'][2] = 'ftp: //livnehpublicstorage.colorado.edu/public/Livneh.2013.CONUS.Dataset/readme.txt'
    
    tmp['dailyvic_livneh2013']['start_date'] = '1915-01-01'
    tmp['dailyvic_livneh2013']['end_date'] = '2011-12-31'
    tmp['dailyvic_livneh2013']['temporal_resolution'] = 'D'
    tmp['dailyvic_livneh2013']['delimiter'] = '\t'
    tmp['dailyvic_livneh2013']['variable_list'] = ['YEAR', 'MONTH', 'DAY', 'EVAP', 'RUNOFF', 'BASEFLOW', 'SMTOP', 'SMMID', 'SMBOT', 'SWE', 'WDEW', 'SENSIBLE', 'LATENT', 'GRNDFLUX', 'RNET', 'RADTEMP', 'PREC']
    
    tmp['dailyvic_livneh2013']['variable_info']={}
    tmp['dailyvic_livneh2013']['variable_info']['YEAR']={}
    tmp['dailyvic_livneh2013']['variable_info']['YEAR']={'desc': 'year', 
                                                         'dtypes': 'int8', 'units': 'yr'}
    tmp['dailyvic_livneh2013']['variable_info']['MONTH']={}
    tmp['dailyvic_livneh2013']['variable_info']['MONTH']={'desc': 'month', 
                                                          'dtypes': 'int8', 'units': 'mo'}
    tmp['dailyvic_livneh2013']['variable_info']['DAY']={}
    tmp['dailyvic_livneh2013']['variable_info']['DAY']={'desc': 'day', 
                                                        'dtypes': 'int8', 'units': 'day'}
    tmp['dailyvic_livneh2013']['variable_info']['EVAP']={}
    tmp['dailyvic_livneh2013']['variable_info']['EVAP']={'desc': 'Total ET rate-- includes Canopy,  Sub-canopy Evaporation,  Transpiration,  and Snow Sublimation', 
                                                         'dtypes': 'float64', 'units': 'mm/s'}
    tmp['dailyvic_livneh2013']['variable_info']['RUNOFF']={}
    tmp['dailyvic_livneh2013']['variable_info']['RUNOFF']={'desc': 'Runoff', 
                                                           'dtypes': 'float64', 'units': 'mm/s'}
    tmp['dailyvic_livneh2013']['variable_info']['BASEFLOW']={}
    tmp['dailyvic_livneh2013']['variable_info']['BASEFLOW']={'desc': 'Baseflow', 
                                                             'dtypes': 'float64', 'units': 'mm/s'}
    tmp['dailyvic_livneh2013']['variable_info']['SMTOP']={}
    tmp['dailyvic_livneh2013']['variable_info']['SMTOP']={'desc': 'Soil moisture top layer', 
                                                          'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2013']['variable_info']['SMMID']={}
    tmp['dailyvic_livneh2013']['variable_info']['SMMID']={'desc': 'Soil moisture middle layer', 
                                                          'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2013']['variable_info']['SMBOT']={}
    tmp['dailyvic_livneh2013']['variable_info']['SMBOT']={'desc': 'Soil moisture bottom layer', 
                                                          'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2013']['variable_info']['SWE']={}
    tmp['dailyvic_livneh2013']['variable_info']['SWE']={'desc': 'Snow water equivalent (SWE)', 
                                                        'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2013']['variable_info']['WDEW']={}
    tmp['dailyvic_livneh2013']['variable_info']['WDEW']={'desc': 'Canopy water', 
                                                         'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2013']['variable_info']['SENSIBLE']={}
    tmp['dailyvic_livneh2013']['variable_info']['SENSIBLE']={'desc': 'Net sensible heat flux', 
                                                             'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2013']['variable_info']['LATENT']={}
    tmp['dailyvic_livneh2013']['variable_info']['LATENT']={'desc': 'Net latent heat flux', 
                                                           'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2013']['variable_info']['GRNDFLUX']={}
    tmp['dailyvic_livneh2013']['variable_info']['GRNDFLUX']={'desc': 'Net heat flux into ground', 
                                                             'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2013']['variable_info']['RNET']={}
    tmp['dailyvic_livneh2013']['variable_info']['RNET']={'desc': 'Net downward radiation flux', 
                                                         'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2013']['variable_info']['RADTEMP']={}
    tmp['dailyvic_livneh2013']['variable_info']['RADTEMP']={'desc': 'Mean radiative surface temperature', 
                                                            'dtypes': 'float64', 'units': 'K'}
    tmp['dailyvic_livneh2013']['variable_info']['PREC']={}
    tmp['dailyvic_livneh2013']['variable_info']['PREC']={'desc': 'Incoming precipitation rate', 
                                                         'dtypes': 'float64', 'units': 'mm/s'}

    """
    Daily Meteorology - bias-corrected Livneh et al. 2013
    """
    tmp['dailymet_bclivneh2013']={}
    tmp['dailymet_bclivneh2013']['spatial_resolution'] = '1/16-degree'
    tmp['dailymet_bclivneh2013']['web_protocol'] = 'http'
    tmp['dailymet_bclivneh2013']['domain'] = 'cses.washington.edu'
    tmp['dailymet_bclivneh2013']['subdomain'] = '/rocinante/Livneh/bcLivneh_WWA_2013/forcings_ascii/'
    tmp['dailymet_bclivneh2013']['decision_steps'] = ''
    tmp['dailymet_bclivneh2013']['filename_structure'] = 'data_{LAT}_{LONG}'
    tmp['dailymet_bclivneh2013']['file_format'] = 'ASCII'
    
    tmp['dailymet_bclivneh2013']['reference']={}
    tmp['dailymet_bclivneh2013']['reference'][1] = 'Livneh,  B.,  E. A. Rosenberg,  C. Lin,  B. Nijssen,  V. Mishra,  K. M. Andreadis,  E. P. Maurer,  and D. P. Lettenmaier,  2013:  A Long-Term Hydrologically Based Dataset of Land Surface Fluxes and States for the Conterminous United States:  Update and Extensions. J. Climate,  26,  9384-9392.'
    tmp['dailymet_bclivneh2013']['reference'][2]='ftp: //livnehpublicstorage.colorado.edu/public/Livneh.2013.CONUS.Dataset/readme.txt'
    
    tmp['dailymet_bclivneh2013']['start_date'] = '1915-01-01'
    tmp['dailymet_bclivneh2013']['end_date'] = '2011-12-31'
    tmp['dailymet_bclivneh2013']['temporal_resolution'] = 'D'
    tmp['dailymet_bclivneh2013']['delimiter'] = '\t'
    tmp['dailymet_bclivneh2013']['variable_list'] = ['PRECIP',  'TMAX',  'TMIN',  'WINDSPD']

    tmp['dailymet_bclivneh2013']['variable_info']={}
    tmp['dailymet_bclivneh2013']['variable_info']['PRECIP']={}
    tmp['dailymet_bclivneh2013']['variable_info']['PRECIP']={'desc': 'daily precipitation (mm)', 
                                                             'dtypes': 'float64', 'units': 'mm'}
    tmp['dailymet_bclivneh2013']['variable_info']['TMAX']={}
    tmp['dailymet_bclivneh2013']['variable_info']['TMAX']={'desc': 'daily maximum temperature (C)', 
                                                           'dtypes': 'float64', 'units': 'C'}
    tmp['dailymet_bclivneh2013']['variable_info']['TMIN']={}
    tmp['dailymet_bclivneh2013']['variable_info']['TMIN']={'desc': 'daily minimum temperature (C)', 
                                                           'dtypes': 'float64', 'units': 'C'}
    tmp['dailymet_bclivneh2013']['variable_info']['WINDSPD']={}
    tmp['dailymet_bclivneh2013']['variable_info']['WINDSPD']={'desc': 'daily mean wind speed (m/s)', 
                                                              'dtypes': 'float64', 'units': 'm/s'}
    
    """
    Daily VIC - Livneh et al.,  2015
    """
    tmp['dailyvic_livneh2015']={}
    tmp['dailyvic_livneh2015']['spatial_resolution']='1/16-degree'
    tmp['dailyvic_livneh2015']['web_protocol']='ftp'
    tmp['dailyvic_livneh2015']['domain']='192.12.137.7'
    tmp['dailyvic_livneh2015']['subdomain']='/pub/dcp/archive/OBS/livneh2014.1_16deg/VIC.ASCII/'
    tmp['dailyvic_livneh2015']['decision_steps']='files organized by Latitude'
    tmp['dailyvic_livneh2015']['filename_structure']='Fluxes_Livneh_NAmerExt_15Oct2014_{LAT}_{LONG}'
    tmp['dailyvic_livneh2015']['file_format']='bz2-compressed ASCII'
    
    tmp['dailyvic_livneh2015']['reference']={}
    tmp['dailyvic_livneh2015']['reference'][1]="Livneh B.,  T.J. Bohn,  D.S. Pierce,  F. Munoz-Ariola,  B. Nijssen,  R. Vose,  D. Cayan,  and L.D. Brekke,  2015:  A spatially comprehensive,  hydrometeorological data set for Mexico,  the U.S.,  and southern Canada 1950-2013,  Nature Scientific Data,  5: 150042,  doi: 10.1038/sdata.2015.42."
    tmp['dailyvic_livneh2015']['reference'][2]="ftp: //192.12.137.7/pub/dcp/archive/OBS/livneh2014.1_16deg/README.Livneh.Grids.txt.v3.txt"
    
    tmp['dailyvic_livneh2015']['start_date']='1950-01-01'
    tmp['dailyvic_livneh2015']['end_date']='2013-12-31'
    tmp['dailyvic_livneh2015']['temporal_resolution']='D'
    tmp['dailyvic_livneh2015']['delimiter']='\t'
    tmp['dailyvic_livneh2015']['variable_list']=['YEAR', 'MONTH', 'DAY', 'EVAP', 'RUNOFF', 'BASEFLOW', 'SMTOP', 'SMMID', 'SMBOT', 'SWE', 'WDEW', 'SENSIBLE', 'LATENT', 'GRNDFLUX', 'RNET', 'PETTALL', 'PETSHORT', 'PETNATVEG']
    
    tmp['dailyvic_livneh2015']['variable_info']={}
    tmp['dailyvic_livneh2015']['variable_info']['YEAR']={}
    tmp['dailyvic_livneh2015']['variable_info']['YEAR']={'desc': 'year', 
                                                         'dtypes': 'int8', 'units': 'yr'}
    tmp['dailyvic_livneh2015']['variable_info']['MONTH']={}
    tmp['dailyvic_livneh2015']['variable_info']['MONTH']={'desc': 'month', 
                                                          'dtypes': 'int8', 'units': 'mo'}
    tmp['dailyvic_livneh2015']['variable_info']['DAY']={}
    tmp['dailyvic_livneh2015']['variable_info']['DAY']={'desc': 'day', 
                                                        'dtypes': 'int8', 'units': 'day'}
    tmp['dailyvic_livneh2015']['variable_info']['EVAP']={}
    tmp['dailyvic_livneh2015']['variable_info']['EVAP']={'desc': 'Total ET rate-- includes Canopy,  Sub-canopy Evaporation,  Transpiration,  and Snow Sublimation', 
                                                         'dtypes': 'float64', 'units': 'mm/day'}
    tmp['dailyvic_livneh2015']['variable_info']['RUNOFF']={}
    tmp['dailyvic_livneh2015']['variable_info']['RUNOFF']={'desc': 'Runoff', 
                                                           'dtypes': 'float64', 'units': 'mm/day'}
    tmp['dailyvic_livneh2015']['variable_info']['BASEFLOW']={}
    tmp['dailyvic_livneh2015']['variable_info']['BASEFLOW']={'desc': 'Baseflow', 
                                                             'dtypes': 'float64', 'units': 'mm/day'}
    tmp['dailyvic_livneh2015']['variable_info']['SMTOP']={}
    tmp['dailyvic_livneh2015']['variable_info']['SMTOP']={'desc': 'Soil moisture top layer', 
                                                          'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2015']['variable_info']['SMMID']={}
    tmp['dailyvic_livneh2015']['variable_info']['SMMID']={'desc': 'Soil moisture middle layer', 
                                                          'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2015']['variable_info']['SMBOT']={}
    tmp['dailyvic_livneh2015']['variable_info']['SMBOT']={'desc': 'Soil moisture bottom layer', 
                                                          'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2015']['variable_info']['SWE']={}
    tmp['dailyvic_livneh2015']['variable_info']['SWE']={'desc': 'Snow water equivalent (SWE)', 
                                                        'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2015']['variable_info']['WDEW']={}
    tmp['dailyvic_livneh2015']['variable_info']['WDEW']={'desc': 'Canopy water', 
                                                         'dtypes': 'float64', 'units': 'mm'}
    tmp['dailyvic_livneh2015']['variable_info']['SENSIBLE']={}
    tmp['dailyvic_livneh2015']['variable_info']['SENSIBLE']={'desc': 'Net sensible heat flux', 
                                                             'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2015']['variable_info']['LATENT']={}
    tmp['dailyvic_livneh2015']['variable_info']['LATENT']={'desc': 'Net latent heat flux', 
                                                           'dtypes': 'float64',  'units': 'W/m^2'}
    tmp['dailyvic_livneh2015']['variable_info']['GRNDFLUX']={}
    tmp['dailyvic_livneh2015']['variable_info']['GRNDFLUX']={'desc': 'Net heat flux into ground', 
                                                             'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2015']['variable_info']['RNET']={}
    tmp['dailyvic_livneh2015']['variable_info']['RNET']={'desc': 'Net downward radiation flux', 
                                                         'dtypes': 'float64', 'units': 'W/m^2'}
    tmp['dailyvic_livneh2015']['variable_info']['PETTALL']={}
    tmp['dailyvic_livneh2015']['variable_info']['PETTALL']={'desc': 'Potential Evapotranspiration from tall crop (Alfalfa)', 
                                                            'dtypes': 'float64', 'units': 'mm/day'}
    tmp['dailyvic_livneh2015']['variable_info']['PETSHORT']={}
    tmp['dailyvic_livneh2015']['variable_info']['PETSHORT']={'desc': 'Potential Evapotranspiration from short crop (Grass)', 
                                                             'dtypes': 'float64', 'units': 'mm/day'}
    tmp['dailyvic_livneh2015']['variable_info']['PETNATVEG']={}
    tmp['dailyvic_livneh2015']['variable_info']['PETNATVEG']={'desc': 'Potential Evapotranspiration from current vegetation', 
                                                              'dtypes': 'float64', 'units': 'mm/day'}

    """
    Daily wrf-nnrp - Salathe et al.,  2014
    """
    tmp['dailywrf_salathe2014']={}
    tmp['dailywrf_salathe2014']['spatial_resolution'] = '1/16-degree'
    tmp['dailywrf_salathe2014']['web_protocol'] = 'http'
    tmp['dailywrf_salathe2014']['domain'] = 'cses.washington.edu/'
    tmp['dailywrf_salathe2014']['subdomain'] = '/rocinante/WRF/NNRP/vic_16d/WWA_1950_2010/raw/forcings_ascii/'
    tmp['dailywrf_salathe2014']['decision_steps'] = ''
    tmp['dailywrf_salathe2014']['filename_structure'] = 'data_{LAT}_{LONG}'
    tmp['dailywrf_salathe2014']['file_format'] = 'ASCII'

    tmp['dailywrf_salathe2014']['reference']={}
    tmp['dailywrf_salathe2014']['reference'][1]=u'Salathé Jr EP,  Hamlet AF,  Mass CF,  Lee SY,  Stumbaugh M,  Steed R. Estimates of twenty-first-century flood risk in the Pacific Northwest based on regional climate model simulations. Journal of Hydrometeorology. 2014 Oct;15(5): 1881-99. DOI:  10.1175/JHM-D-13-0137.1'
    tmp['dailywrf_salathe2014']['reference'][2]='http://cses.washington.edu/rocinante/WRF/README'

    tmp['dailywrf_salathe2014']['start_date'] = '1950-01-01'
    tmp['dailywrf_salathe2014']['end_date'] = '2010-12-31'
    tmp['dailywrf_salathe2014']['temporal_resolution'] = 'D'
    tmp['dailywrf_salathe2014']['delimiter'] = '\\s+'
    tmp['dailywrf_salathe2014']['variable_list'] = ['PRECIP',  'TMAX',  'TMIN',  'WINDSPD']

    tmp['dailywrf_salathe2014']['variable_info']={}
    tmp['dailywrf_salathe2014']['variable_info']['PRECIP']={}
    tmp['dailywrf_salathe2014']['variable_info']['PRECIP']={'desc': 'Daily accumulated precipitation', 
                                                            'dtypes': 'float64', 'units': 'mm'}
    tmp['dailywrf_salathe2014']['variable_info']['TMAX']={}
    tmp['dailywrf_salathe2014']['variable_info']['TMAX']={'desc': 'Maximum temperature at 2m', 
                                                          'dtypes': 'float64', 'units': 'C'}
    tmp['dailywrf_salathe2014']['variable_info']['TMIN']={}
    tmp['dailywrf_salathe2014']['variable_info']['TMIN']={'desc': 'Minimum temperature at 2m', 
                                                          'dtypes': 'float64', 'units': 'C'}
    tmp['dailywrf_salathe2014']['variable_info']['WINDSPD']={}
    tmp['dailywrf_salathe2014']['variable_info']['WINDSPD']={'desc': 'Wind Speed', 
                                                             'dtypes': 'float64', 'units': 'm/s'}

    """
    Daily WRF-nnrp - bias-corrected Salathe et al.,  2014
    """
    tmp['dailywrf_bcsalathe2014']={}
    tmp['dailywrf_bcsalathe2014']['spatial_resolution'] = '1/16-degree'
    tmp['dailywrf_bcsalathe2014']['web_protocol'] = 'http'
    tmp['dailywrf_bcsalathe2014']['domain'] = 'cses.washington.edu/'
    tmp['dailywrf_bcsalathe2014']['subdomain'] = '/rocinante/WRF/NNRP/vic_16d/WWA_1950_2010/bc/forcings_ascii/'
    tmp['dailywrf_bcsalathe2014']['decision_steps'] = ''
    tmp['dailywrf_bcsalathe2014']['filename_structure'] = 'data_{LAT}_{LONG}'
    tmp['dailywrf_bcsalathe2014']['file_format'] = 'ASCII'

    tmp['dailywrf_bcsalathe2014']['reference']={}
    tmp['dailywrf_bcsalathe2014']['reference'][1]='Salathé Jr EP,  Hamlet AF,  Mass CF,  Lee SY,  Stumbaugh M,  Steed R. Estimates of twenty-first-century flood risk in the Pacific Northwest based on regional climate model simulations. Journal of Hydrometeorology. 2014 Oct;15(5): 1881-99. DOI:  10.1175/JHM-D-13-0137.1'
    tmp['dailywrf_bcsalathe2014']['reference'][2]='http://cses.washington.edu/rocinante/WRF/README'

    tmp['dailywrf_bcsalathe2014']['start_date'] = '1950-01-01'
    tmp['dailywrf_bcsalathe2014']['end_date'] = '2010-12-31'
    tmp['dailywrf_bcsalathe2014']['temporal_resolution'] = 'D'
    tmp['dailywrf_bcsalathe2014']['delimiter'] = '\\s+'
    tmp['dailywrf_bcsalathe2014']['variable_list'] = ['PRECIP',  'TMAX',  'TMIN',  'WINDSPD']

    tmp['dailywrf_bcsalathe2014']['variable_info']={}
    tmp['dailywrf_bcsalathe2014']['variable_info']['PRECIP']={}
    tmp['dailywrf_bcsalathe2014']['variable_info']['PRECIP']={'desc': 'Daily accumulated precipitation', 
                                                              'dtypes': 'float64', 'units': 'mm'}
    tmp['dailywrf_bcsalathe2014']['variable_info']['TMAX']={}
    tmp['dailywrf_bcsalathe2014']['variable_info']['TMAX']={'desc': 'Maximum temperature at 2m', 
                                                            'dtypes': 'float64', 'units': 'C'}
    tmp['dailywrf_bcsalathe2014']['variable_info']['TMIN']={}
    tmp['dailywrf_bcsalathe2014']['variable_info']['TMIN']={'desc': 'Minimum temperature at 2m', 
                                                            'dtypes': 'float64', 'units': 'C'}
    tmp['dailywrf_bcsalathe2014']['variable_info']['WINDSPD']={}
    tmp['dailywrf_bcsalathe2014']['variable_info']['WINDSPD']={'desc': 'Wind Speed', 
                                                               'dtypes': 'float64', 'units': 'm/s'}

    # output the file
    #json.dump(tmp,  open('ogh_meta.json',  'w'),  ensure_ascii=False)

    # perform a test read
    #json.load(open('ogh_meta.json'))

    return(tmp)
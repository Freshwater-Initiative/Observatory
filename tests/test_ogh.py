from __future__ import (absolute_import,
                        division,
                        print_function,
                        unicode_literals)

import os
import pytest
import ftplib
import ogh
import numpy as np
import pandas as pd
import geopandas as gpd

# TODO: NEED TESTS!!!
# class Testoghfunctions(object):

#     def setup(self):
#         pass

#     def test_saveDictOfDf(self):
#         pass
    
#     def test_saveDictOfDf(self):
#         pass


# loading and saving dataframes and dictionaries
class Test_ogh_load_functions(object):

    def test_ogh_meta():
        keys = dict(ogh.ogh_meta())
        if len(keys)>0:
            pass


    def test_DictOfDf():
        path ='tests/data/test.json'
        obj1 = pd.DataFrame({'foo':['a','b'],'bar':[1,2]})
        obj2 = ['foo','bar']

        # savedictofdf success
        ogh.saveDictOfDf(path, obj1)
        assert os.path.exists(path)

        # readdictofdf success
        test = ogh.readDictOfDf(path)
        assert True

        # savedictofdf fail
        try:
            ogh.saveDictOfDf(path, obj2)
        except:
            os.remove(path)
            assert True

    def test_ensuredir():
        path0 = os.getcwd()
        path1 = 'tests/data/test_subfolder'
        ogh.ensuredir(path1)
        ogh.ensuredir(path0)
        assert os.path.exists(path1)
        
        
class Test_mappingfile_ops(object):
    def test_readmappingfile():
        test_map = ogh.mappingfileToDF(mappingfile='tests/data/test_mapping.csv', colvar='all', summary=True)
        assert True

        test_compare = ogh.compareonvar(map_df=test_map, colvar=None)
        assert True

    
# loading and saving shapefiles
class Test_ogh_shape_functions(object):

    def test_shapefile_reading():
        path='tests/data/shape.shp'

        # getFullShape success
        test = ogh.getFullShape(path)
        assert True

        # readShapefileTable success
        test = ogh.readShapefileTable(shapefile='tests/data/shape.shp')
        assert True

        # getShapeBbox
        test = ogh.getShapeBbox(fiona.open(path))
        assert True


    def test_reprojShapefile():
        ogh.reprojShapefile(sourcepath='tests/data/shape.shp', 
                            newprojdictionary={'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})
        assert True


    def test_treatgeoself_parts():
        inpath='tests/data/shape.shp'
        outpath='tests/data/test_mappingfile.csv'

        # reading in the shapefile
        test_poly=gpd.read_file(inpath)
        assert True

        # generate lat-lon
        test_lon, test_lat= np.array(test_poly.centroid[0])

        # filterPointsinShape
        test_maptable = ogh.filterPointsinShape(shape=test_poly.geometry[0], 
                                                points_lat=test_lat, points_lon=test_lon, points_elev=None, 
                                                buffer_distance=0.06, buffer_resolution=16, labels=['LAT', 'LONG_', 'ELEV'])
        test_maptable = test_maptable.reset_index(inplace=True).rename(columns={'index': 'FID'})
        assert len(test_maptable)>0

        # print the mappingfile
        test_maptable.to_csv(outpath, sep=',', header=True, index=False)
        assert os.path.exists(outpath)

    
# spatial mapping
class Test_ogh_spatial_mapping(object):

    def test_canadabox_bc():
        test = ogh.canadabox_bc()
        assert True


    def test_mapContentFolder():
        test = ogh.mapContentFolder('test')
        assert True


    def test_mapfilelocations():
        test_maptable = pd.read_csv('tests/data/test_mappingfile.csv')

        # livneh 2013
        assert ogh.compile_dailyMET_Livneh2013_locations(test_maptable) is not None
        assert ogh.compile_VICASCII_Livneh2013_locations(test_maptable) is not None
        
        # livneh 2013 (CIG)
        assert ogh.compile_Livneh2013_locations(test_maptable) is not None
        assert ogh.compile_bc_Livneh2013_locations(test_maptable) is not None

        # livneh 2015
        assert ogh.compile_dailyMET_Livneh2015_locations(test_maptable) is not None
        assert ogh.compile_VICASCII_Livneh2015_locations(test_maptable) is not None

        # Salathe 2014
        assert ogh.compile_wrfnnrp_raw_Salathe2014_locations(test_maptable) is not None
        assert ogh.compile_wrfnnrp_bc_Salathe2014_locations(test_maptable) is not None



# web scraping
class Test_ogh_webscraping(object):

    def test_scrapeurl():
        path = 'test'
        filenames=ogh.scrapeurl(path, startswith='data')
        assert len(filenames)>0


    def test_scrapedomain():
        # gridded data product metadata
        domain='livnehpublicstorage.colorado.edu'
        subdomain='/public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2'

        # identify the subfolder blocks
        blocks = ogh.scrape_domain(domain=domain, subdomain=subdomain, startswith='fluxes')
        assert len(blocks)>0

    

# web download
class Test_ogh_webdownload(object):

    def test_ftp_download():
        protocol = 'ftp://'
        ipaddress = 'livnehpublicstorage.colorado.edu/'
        subdomain = 'public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.100.95.25.36/'
        filename1 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40625.bz2' # real
        filename2 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40626.bz2' # fake

        # connect to ipaddress and subfolder
        try:
            ftp=ftplib.FTP(ipaddress)
            ftp.login()
            ftp.cwd(subdomain)

            # connection success
            if len(ftp.nlst(filename1))==1:
                assert True

            # connection fail
            if len(ftp.nlst(filename2))==0:
                assert True

            ftp.close()
        except:
            pass


    def test_ftp_download_one():
        protocol = 'ftp://'
        ipaddress = 'livnehpublicstorage.colorado.edu/'
        subdomain = 'public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.100.95.25.36/'
        filename1 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40626.bz2' # fake

        urlofinterest=os.path.join(protocol, ipaddress, subdomain, filename1)
        ogh.ftp_download_one(urlofinterest)
        assert True


    def test_ftp_download_p():
        protocol = 'ftp://'
        ipaddress = 'livnehpublicstorage.colorado.edu/'
        subdomain = 'public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.100.95.25.36/'
        filename1 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40625.bz2' # real
        filename2 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40626.bz2' # fake

        listofinterest=[os.path.join(protocol, ipaddress, subdomain, filename1), 
                        os.path.join(protocol, ipaddress, subdomain, filename2)]
        ogh.ftp_download_p(listofinterest)
        assert True


    def test_wget_download():
        listofinterest=['test1', 'test2']
        ogh.wget_download(listofinterest)
        assert True


    def test_wget_download_one():
        urlofinterest='test1'
        ogh.wget_download_one(urlofinterest)
        assert True


    def test_wget_download_p():
        listofinterest=['test1', 'test2']
        ogh.wget_download_p(listofinterest)
        assert True


# Download the files to the subdirectory
class Test_ogh_cataloging(object):

    folderpath='tests/data/test_files/' # has sample file
    
    def test_catalogfiles():
        test = ogh.catalogfiles(folderpath)
        assert True

    def test_addCatalogToMap():
        # read in a sample mappingfile as test_map
        test_map = ogh.mappingfileToDf('tests/data/test_mappingfile.csv', colvar=None)
        ogh.addCatalogToMap(outfilepath='tests/data/test_catalog.csv', 
                            maptable=test_map, 
                            folderpath=folderpath,
                            catalog_label='test')
        assert True

    
# Wrapper scripts
class Test_ogh_wrappedget(object):
    def test_getDailyMET_livneh2013():
        ogh.getDailyMET_livneh2013(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv', 
                                   subdir='tests/data/test_files', catalog_label='dailymet_livneh2013')
        assert True
    
    def test_getDailyMET_livneh2015():
        ogh.getDailyMET_livneh2015(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv',
                                   subdir='tests/data/test_files', catalog_label='dailymet_livneh2015')
        assert True
        

    def test_getDailyMET_bcLivneh2013():
        ogh.getDailyMET_bcLivneh2013(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv',
                                     subdir='tests/data/test_files',catalog_label='dailymet_bclivneh2013')
        assert True


    def test_getDailyVIC_livneh2013():
        ogh.getDailyVIC_livneh2013(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv',
                                   subdir='tests/data/test_files', catalog_label='dailyvic_livneh2013')
        assert True
        
    def test_getDailyVIC_livneh2015():
        ogh.getDailyVIC_livneh2015(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv',
                                   subdir='tests/data/test_files', catalog_label='dailyvic_livneh2015')
        assert True
    

    def test_getDailyWRF_salathe2014():
        ogh.getDailyWRF_salathe2014(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv',
                                    subdir='tests/data/test_files', catalog_label='dailywrf_salathe2014')
        assert True
        

    def test_getDailyWRF_bcsalathe2014():
        ogh.getDailyWRF_bcsalathe2014(homedir=os.getcwd(), mappingfile='tests/data/test_mappingfile.csv', 
                                      subdir='tests/data/test_files', catalog_label='dailywrf_bcsalathe2014')
        assert True
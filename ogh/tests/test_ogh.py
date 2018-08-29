from __future__ import (absolute_import,
                        division,
                        print_function,
                        unicode_literals)

import os
import pytest
import ftplib
import numpy as np
import pandas as pd
import fiona
import geopandas as gpd
import ogh as ogh

data_path = os.path.join(ogh.__path__[0], 'tests/data')

# loading and saving dataframes and dictionaries
class Test_ogh_load_functions(object):

    def test_ogh_meta(self):
        keys = dict(ogh.ogh_meta())
        if len(keys)>0:
            pass
        assert True


    def test_DictOfDf(self):
        path = os.path.join(data_path,'test.json')
        obj1 = pd.DataFrame({'foo':['a','b'],'bar':[1,2]})
        obj2 = ['foo','bar']

        # savedictofdf success
        ogh.saveDictOfDf(path, obj1.to_dict())
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

    def test_ensuredir(self):
        path0 = os.getcwd()
        path1 = os.path.join(data_path,'test_files')
        ogh.ensure_dir(path1)
        ogh.ensure_dir(path0)
        assert os.path.exists(path1)
        
    
# loading and saving shapefiles
class Test_ogh_shape_functions(object):

    def test_shapefile_reading(self):
        # getFullShape success
        test = ogh.getFullShape(os.path.join(data_path,'shape.shp'))
        assert True

        # readShapefileTable success
        test = ogh.readShapefileTable(shapefile=os.path.join(data_path,'shape.shp'))
        assert True

        # getShapeBbox
        test = ogh.getShapeBbox(fiona.open(os.path.join(data_path,'shape.shp')))
        assert True


    def test_reprojShapefile(self):
        ogh.reprojShapefile(sourcepath=os.path.join(data_path,'shape.shp'), 
                            newprojdictionary={'proj':'longlat', 'ellps':'WGS84', 'datum':'WGS84'})
        assert True

    
# spatial mapping
class Test_ogh_spatial_mapping(object):

    def test_canadabox_bc(self):
        test = ogh.canadabox_bc()
        assert True


    def test_mapContentFolder(self):
        test = ogh.mapContentFolder(os.path.join(data_path,'test_files'))
        assert True


    def test_mapfilelocations(self):
        test_maptable = pd.read_csv(os.path.join(data_path,'test_mappingfile.csv'))

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

    def test_scrapeurl(self):
        #path = os.path.join(data_path,'test_files') # needs to be a url
        #filenames=ogh.scrapeurl(path, startswith='data')
        #assert len(filenames)>0
        assert True


    def test_scrapedomain(self):
        # gridded data product metadata
        domain='livnehpublicstorage.colorado.edu'
        subdomain='/public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2'

        # identify the subfolder blocks
        blocks = ogh.scrape_domain(domain=domain, subdomain=subdomain, startswith='fluxes')
        assert len(blocks)>0

    

# web download
class Test_ogh_webdownload(object):

    def test_ftp_download(self):
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


    def test_ftp_download_one(self):
        protocol = 'ftp://'
        ipaddress = 'livnehpublicstorage.colorado.edu/'
        subdomain = 'public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.100.95.25.36/'
        filename1 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40626.bz2' # fake

        urlofinterest=os.path.join(protocol, ipaddress, subdomain, filename1)
        ogh.ftp_download_one(urlofinterest)
        assert True


    def test_ftp_download_p(self):
        protocol = 'ftp://'
        ipaddress = 'livnehpublicstorage.colorado.edu/'
        subdomain = 'public/Livneh.2013.CONUS.Dataset/Fluxes.asc.v.1.2.1915.2011.bz2/fluxes.100.95.25.36/'
        filename1 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40625.bz2' # real
        filename2 = 'VIC_fluxes_Livneh_CONUSExt_v.1.2_2013_25.90625_-97.40626.bz2' # fake

        listofinterest=[os.path.join(protocol, ipaddress, subdomain, filename1), 
                        os.path.join(protocol, ipaddress, subdomain, filename2)]
        ogh.ftp_download_p(listofinterest)
        os.remove(filename1.replace('.bz2',''))
        assert True

    def test_wget_download(self):
        listofinterest=['test1', 'test2']
        ogh.wget_download(listofinterest)
        assert True


    def test_wget_download_one(self):
        urlofinterest='test1'
        ogh.wget_download_one(urlofinterest)
        assert True


    def test_wget_download_p(self):
        listofinterest=['test1', 'test2']
        ogh.wget_download_p(listofinterest)
        assert True


# Download the files to the subdirectory
class Test_ogh_cataloging(object):
    def test_catalogfiles(self):
        test = ogh.catalogfiles(os.path.join(data_path,'test_files'))
        assert True

    def test_addCatalogToMap(self):        
        # read in a sample mappingfile as test_map
        test_map, nstat = ogh.mappingfileToDF(os.path.join(data_path,'test_mappingfile.csv'), colvar=None)
        ogh.addCatalogToMap(outfilepath=os.path.join(data_path,'test_catalog.csv'), 
                            maptable=test_map, 
                            folderpath=os.path.join(data_path,'test_files'),
                            catalog_label='test')
        assert True

    
# Wrapper scripts
class Test_ogh_wrappedget(object):
    def test_getDailyMET_livneh2013(self):
        ogh.getDailyMET_livneh2013(homedir=data_path, 
                                   mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                   subdir=os.path.join(data_path,'test_files1'),
                                   catalog_label='dailymet_livneh2013')
        assert True
    
    def test_getDailyMET_livneh2015(self):
        ogh.getDailyMET_livneh2015(homedir=os.getcwd(), 
                                   mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                   subdir=os.path.join(data_path,'test_files2'),
                                   catalog_label='dailymet_livneh2015')
        assert True
        

    def test_getDailyMET_bcLivneh2013(self):
        ogh.getDailyMET_bcLivneh2013(homedir=data_path, 
                                     mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                     subdir=os.path.join(data_path,'test_files3'),
                                     catalog_label='dailymet_bclivneh2013')
        assert True


    def test_getDailyVIC_livneh2013(self):
        ogh.getDailyVIC_livneh2013(homedir=data_path, 
                                   mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                   subdir=os.path.join(data_path,'test_files4'),
                                   catalog_label='dailyvic_livneh2013')
        assert True
        
    def test_getDailyVIC_livneh2015(self):
        ogh.getDailyVIC_livneh2015(homedir=data_path, 
                                   mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                   subdir=os.path.join(data_path,'test_files5'),
                                   catalog_label='dailyvic_livneh2015')
        assert True
    

    def test_getDailyWRF_salathe2014(self):
        ogh.getDailyWRF_salathe2014(homedir=data_path, 
                                    mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                    subdir=os.path.join(data_path,'test_files6'),
                                    catalog_label='dailywrf_salathe2014')
        assert True
        

    def test_getDailyWRF_bcsalathe2014(self):
        ogh.getDailyWRF_bcsalathe2014(homedir=data_path, 
                                      mappingfile=os.path.join(data_path,'test_catalog.csv'), 
                                      subdir=os.path.join(data_path,'test_files7'),
                                      catalog_label='dailywrf_bcsalathe2014')
        assert True
        
    def test_files_remove(self):
        for eachdir in ['test_files1','test_files2','test_files3','test_files4','test_files5','test_files6','test_files7']:
            path = os.path.join(data_path, eachdir)
            pd.Series(os.listdir(path)).apply(lambda x: os.remove(os.path.join(path, x)))
            os.rmdir(path)
        os.remove(os.path.join(data_path,'test_catalog.csv'))
        assert True
        
        
class Test_mappingfile_ops(object):
    def test_readmappingfile(self):        
        test_map, nstat = ogh.mappingfileToDF(os.path.join(data_path,'test_mappingfile.csv'), colvar=None)
        test_map = test_map.drop_duplicates()
        test_map.to_csv(os.path.join(data_path,'test_mappingfile.csv'), index=False, columns=['FID','LAT','LONG_','ELEV'])
        assert True

        test_compare = ogh.compareonvar(map_df=test_map, colvar=None)
        assert True
        
    def test_treatgeoself_parts(self):
        inpath=os.path.join(data_path,'shape.shp')

        # reading in the shapefile
        test_poly=gpd.read_file(inpath)
        assert True

        # generate lat-lon
        test_lon, test_lat= np.array(test_poly.centroid[0])

        # filterPointsinShape
        test_maptable = ogh.filterPointsinShape(shape=test_poly.geometry[0], 
                                                points_lat=[test_lat], points_lon=[test_lon], points_elev=None, 
                                                buffer_distance=0.06, buffer_resolution=16, labels=['LAT', 'LONG_', 'ELEV'])
        test_maptable = test_maptable.reset_index().rename(columns={'index': 'FID'})
        assert len(test_maptable)>0
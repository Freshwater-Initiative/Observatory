# unit test notes

Core unit tests were compiled to test functions and operations related to use-cases 1-4. Additional unit tests are needed for use-cases 5+.

Due to issues interacting with ftp web domains, a number of unit tests related to web-scraping or parallel file requests from ftp domains have been deactivated. These deactivated unit tests have been relate to ftp domains, specifically those coming from the Livneh et al. 2013 ('ftp://livnehpublicstorage.colorado.edu') gridded data products. New unit tests are needed.

Additional unit tests are also needed for the ogh_xarray_landlab class, such as for landlab rastermodelgrid, netcdf input data formats, netcdf-to-ascii and grid/regrid operations.
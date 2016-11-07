#!/bin/python
# pylint: disable=C0103

"""Download and extract file"""

import urllib.request
import zipfile
import os

url = "http://www.ngdc.noaa.gov/mgg/global/relief/ETOPO2/ETOPO2v2-2006/"\
        "ETOPO2v2g/netCDF/ETOPO2v2g_f4_netCDF.zip"
directory = 'data'
fname = os.path.basename(url)
urllib.request.urlretrieve(url, directory + '/' + fname)

zf = zipfile.ZipFile(directory + '/' + fname)
zf.extractall(directory)

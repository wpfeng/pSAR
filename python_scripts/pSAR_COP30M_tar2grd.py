#!/usr/bin/env python
#
#
#
#
import os
import sys
import glob
#
def roi2copdem_tars(lon0,lon1,lat0,lat1):
    #
    tars = []
    #
    dem_prefix = 'Copernicus_DSM_10_'
    #
    #
    s_lat = 'N'
    s_lon = 'E'
    # dem_prefix = 'Copernicus_DSM_30_'
    for clon in range(lon0,lon1,1):
        for clat in range(lat0,lat1,1):
            #
            if clat <0:
                s_lat = 'S'
            if clon <0:
                s_lon = 'W'
            #
            outname = '%s%s%2d_00_%s%03d_00.tar' % (dem_prefix,s_lat,abs(clat),s_lon,abs(clon))
            tars.append(outname)
            #
        #
    #
    return tars
#

#
if len(sys.argv)<2:
    #
    helpstr = \
            '''
            %s <run> -roi <lon0,lon1,lat0,lat1> -outgrd [dem.grd in default]
            +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            To merge and subset copernicus DEM data based on given information


            Three steps are required to reach the final output data:
            1) untar files from target folder <in_dir>
            2) merge all geotiff files into a merge.tif
            3) subset merge.tif based on given projwin <-roi> to a grd format (netCDF)

            -outgrd output name, dem.grd will be in default

            '''
    print(helpstr % sys.argv[0])
    #
    sys.exit(-1)
    #
#
roi = None
outgrd = 'dem.grd'
#
for i,key in enumerate(sys.argv):
    #
    if key == '-roi':
        roi = [float(croi) for croi in sys.argv[i+1].split(',')]
    if key == '-outgrd':
        outgrd = sys.argv[i+1]
#
if True:
    tmp_folder = 'Copernicus_DEM_merge_DIR'
    #
    if not os.path.exists(tmp_folder):
        os.makedirs(tmp_folder)
        #
    #
    copdem_30m_home = os.environ['COPDEM_HOME']
    #
    #
    if roi is None:
      tars = glob.glob(sys.argv[1]+'/Coperni*DSM*.tar')
    else:
      lon0 = int(roi[0]-1)
      lon1 = int(roi[1]+1)
      lat0 = int(roi[2]-1)
      lat1 = int(roi[3]+1)
      tars = roi2copdem_tars(lon0,lon1,lat0,lat1)
      tars = ['%s/%s' % (copdem_30m_home,ctar) for ctar in tars]
    #
    if len(tars)<1:
        print(" Err: no valid DEM was found!!!")
        sys.exit(-1)
        #
    #
    for ctar in tars:
        #
        outname = tmp_folder+'/%s' % os.path.basename(ctar)
        #
        if not os.path.exists(outname):
            #
            os.system('ln -s %s %s' % (os.path.abspath(ctar),outname))
            #
        #
        #
    #
    # now let us untar all files
    #
    topdir = os.getcwd()
    os.chdir(tmp_folder)
    #
    goStr = 'pSAR_COP30M_untar.py'
    print('  running: %s' % goStr)
    os.system(goStr)
    #
    goStr = ' ls Coper*/DEM/Coper*.tif > dem.list'
    print('  running: %s' % goStr)
    os.system(goStr)
    #
    goStr = 'gdal_merge.py -o merge.tif --optfile dem.list'
    print('  running: %s' % goStr)
    os.system(goStr)
    # gdal_translate -of netCDF -projwin -9.35625 31.8529166667 -6.94791666667 29.7529166667  merge.tif dem.grd
    goStr = 'gdal_translate -of netCDF -projwin %f %f %f %f %s %s' % (roi[0],roi[3],roi[1],roi[2],'merge.tif',outgrd)
    print('  running: %s ' % goStr)
    os.system(goStr)
    #
    if os.path.exists(outgrd):
        #
        os.system('ln -s $PWD/%s %s/%s' % (outgrd,topdir,outgrd))
        #
    #
    #
    os.chdir(topdir)
    #




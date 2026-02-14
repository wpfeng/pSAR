#!/usr/bin/env python
#
#

import datetime
import os
import sys
import numpy as np
import xarray as xr
import matplotlib
import quadtree_ds as qds
import utm
from shapely.geometry import Polygon
from pyproj import Transformer, CRS
#
#
# matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt
#
def calculate_polygon_area(lls):
    #
    """
    Calculate the area of a geographic coordinate polygon (unit: square meters)
    coords: List of polygon vertices in format [(lon1, lat1), (lon2, lat2), ..., (lon1, lat1)]
            Note: The polygon must be closed (first and last points are the same)
    """
    # Calculate center coordinates to determine appropriate UTM zone
    lon_center = np.mean(lls[0:4,0]) 
    lat_center = np.mean(lls[0:4,1])
    #
    #
    # Calculate UTM zone from longitude (关键修复：根据经度计算带号)
    utm_zone = int((lon_center + 180) // 6) + 1  # 带号计算公式
    hemisphere = 'north' if lat_center >= 0 else 'south'  # 判断南北半球

    # 正确的UTM投影参数（使用带号和半球）
    utm_crs = CRS.from_string(f"+proj=utm +zone={utm_zone} +{hemisphere} +datum=WGS84")
    #
    # Define UTM projection based on center coordinates
    # utm_crs = CRS.from_proj4(f"+proj=utm +lon_0={lon_center} +lat_0={lat_center} +datum=WGS84")

    # Transform geographic coordinates (WGS84) to UTM projected coordinates
    transformer = Transformer.from_crs(CRS("EPSG:4326"), utm_crs, always_xy=True)
    projected_coords = [transformer.transform(lls[i,0], lls[i,1]) for i in range(lls.shape[0])]

    # Calculate area using shapely's built-in method
    polygon = Polygon(projected_coords)
    return polygon.area/1e6  # Returns area in square meters
#
def get_utm_info(latitude, longitude):
    """
    Calculate UTM zone number, letter, and related coordinate reference system information.

    Args:
        latitude (float): Latitude in decimal degrees.
        longitude (float): Longitude in decimal degrees.

    Returns:
        dict: A dictionary containing UTM zone number, hemisphere, letter,
              EPSG code, and CRS object.
    """
    # Calculate UTM zone number (1-60)
    zone_number = int((longitude + 180) / 6) + 1
    zone_number = max(1, min(zone_number, 60))  # Clamp to valid range

    # Determine hemisphere and EPSG code
    is_northern = latitude >= 0
    hemisphere = 'north' if is_northern else 'south'
    epsg_base = 32600 if is_northern else 32700
    epsg_code = epsg_base + zone_number

    # Calculate UTM letter designator (latitude band)
    # Mapping of latitude ranges to letters (excluding special cases)
    letter_bands = "CDEFGHJKLMNPQRSTUVWXX"  # X is repeated for 84-80°N and 80-84°S
    letter_index = int((latitude + 80) / 8)

    # Handle special cases and clamp to valid index
    if latitude >= 84:
        letter = 'X'
    elif latitude < -80:
        letter = 'C'  # Southernmost band
    else:
        letter_index = max(0, min(letter_index, len(letter_bands) - 1))
        letter = letter_bands[letter_index]


    return {
        'un': zone_number,
        'hemisphere': hemisphere,
        'ul': letter,
    }
#
def grd_write_xr(data,lons,lats,outgrd):
    #
    data = np.flipud(data)
    #
    da = xr.DataArray(data, coords=[lats,lons], dims=['y', 'x'])
    #
    # set attributes
    # to make sure we can have valid v_min and v_max in GMT
    #
    vmin = np.min(data[~np.isnan(data)])
    vmax = np.max(data[~np.isnan(data)])
    da.attrs['actual_range'] = [vmin,vmax]

    # create DataArray as Dataset
    ds = da.to_dataset(name='z')

    # save NetCDF
    ds.to_netcdf(outgrd)
#
def quad_xy2poly(cmll):
    #
    outpoly = np.zeros([5,2])
    outpoly[0,:] = [cmll[0],cmll[2]]
    outpoly[1,:] = [cmll[0],cmll[3]]
    outpoly[2,:] = [cmll[1],cmll[3]]
    outpoly[3,:] = [cmll[1],cmll[2]]
    outpoly[4,:] = [cmll[0],cmll[2]]
    #
    return outpoly
#
#
def this_help():
    #
    helpstr = \
            '''
            %s
            _______________________________________________________________________
            insar_preprocessing.py <in_grd> 
                                -incfile [e.g. geo_T144.inc.grd]
                                -azifile [e.g. geo_T144.azi.grd]
                                -inc [35 in default, a value could be instead of incfile]
                                -azi [-12 in default, a value could be instead of azifile]
                                -un [utm zone number]
                                -ul [utm zone letter]
                                -maskext [None in default, e.g. lon1,lon2,lat1,lat2]
                                -maxblock [256 in default, should be a number as 2**n, n an integer...]
                                -minblock [4 in default, as above]
                                -std      [0.75 radian ]
                                -wavelength [0.0554 in default]
                                -km         [0 in default, if given, utm coordinates in km will be outputed]
                                -outroot [None in default]
                                -pot   [False in default]
            ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
            This was first written for use in the class. Now it is gradually used in the daily work...
            
            by Wanpeng Feng, @SYSU, 2024/10/24

            added some helpinfo @DAL, Halifax, 2025/02/24

            '''
    #
    print(helpstr % sys.argv[0])
    sys.exit(-1)
    #
def log(in_inputs,logname=None):
    #
    cnow = datetime.datetime.now().strftime('%Y-%m-%dT%H:%M:%S.%f')
    #
    if len(in_inputs)>1:
        #
        if logname is None:
           logname = in_puts[1]+'.log'
        #
        with open(logname,'a') as fid:
            #
            fid.write('%s: ' % cnow)
            #
            for i,key in enumerate(in_inputs):
                #
                if i < len(in_inputs)-1:
                   fid.write('%s ' % in_inputs[i])
                else:
                   fid.write('%s\n' % in_inputs[i])

class InSAR_preprocessing(object):
    #
    def __init__(self,km=0,wavelength=0.0556):
        #
        self.maskdata = []
        self.deramp_data = []
        self.data = []
        self.vmin = []
        self.vmax = []
        self.lons = []
        self.un = None
        self.ul = None
        self.lats = []
        self.gext = []
        self.inc = []
        self.azi = []
        self.wavelength=wavelength
        self.minblock = 128
        self.maxblock = 32
        self.quad_pts = []
        self.quad_xy = []
        self.filename = []
        self.subext   = []
        self.codes = [0,0,0]
        self.utm = []
        self.mode = []
        self.pot = False
        # 
        # input grd may have coordinates in km
        # when self.km == 1
        #
        self.km   = km
        #
    #
    def export_grd(self,outroot):
        #
        grd_write_xr(self.deramp_data, self.lons, self.lats, outroot+'_deramp.grd')
        grd_write_xr(self.data,        self.lons, self.lats, outroot+'_subset.grd')
        #
    #
    def subset(self,ext=None):
        #
        # ext, having [minlon,maxlon,minlat,maxlat]
        if ext is not None:
            #
            flag1 = self.lons>=ext[0]
            flag2 = self.lons<=ext[1]
            flag_lon = flag1 * flag2
            #
            index_lon = np.where(flag_lon)[0]
            x0 = index_lon[0]
            x1 = index_lon[-1]
            #
            flag2 = self.lats>=ext[2]
            flag3 = self.lats<=ext[3]
            flag_lat = flag2 * flag3
            #
            index_lat = np.where(flag_lat)[0]
            _y0 = index_lat[0]
            _y1 = index_lat[-1]
            #
            nlat = flag_lat.shape[0]
            #
            y1 = nlat-_y0
            y0 = nlat-_y1
            #
            self.lons = self.lons[x0:x1]
            self.lats = self.lats[_y0:_y1]
            self.maskdata = self.maskdata[y0:y1,x0:x1]
            self.data = self.data[y0:y1,x0:x1]
            #
            self.deramp_data = self.deramp_data[y0:y1,x0:x1]
            #
            #
            #

    #
    def insarPlot(self,vmin=None,vmax=3.14,cmap='jet'):
        #
        #
        if vmin is None or vmax is None:
            #
            _data = self.data
            _data = _data[~np.isnan(_data)]
            _data = _data[_data!=0]
            vmin = np.min(_data)
            vmax = np.max(_data)
            self.vmin = vmin
            self.vmax = vmax
            #
        if len(self.lons) > 0:
            # 
            # interpolation_method
            #
            subplots_n = np.sum(self.codes)
            #
            #
            plt.subplot(subplots_n,1,1)
            #
            plt.imshow(self.data,interpolation='bilinear',vmin=vmin,vmax=vmax,\
                    extent = self.gext,\
                    cmap=plt.get_cmap(cmap))
            #
            plt.colorbar()
            #
            if self.codes[1] == 1:
                plt.subplot(subplots_n,1,2)
                plt.imshow(self.maskdata,interpolation='bilinear',vmin=vmin,vmax=vmax,\
                    extent = self.gext,\
                    cmap=plt.get_cmap(cmap))
                #
                plt.colorbar()
                #
            if self.codes[2] == 1:
                #
                plt.subplot(subplots_n,1,3)
                plt.imshow(self.deramp_data,interpolation='bilinear',vmin=vmin,vmax=vmax,\
                    extent = self.gext,\
                    cmap=plt.get_cmap(cmap))
                #
                plt.colorbar()
            #
            plt.show()
            #
        else:
            print(" Err: Null data is given. Please read a grd file first!!!")

    #
    def update_pts(self,incfile=None,azifile=None):
        #
        if len(self.quad_pts)>0:
            #

            ds = xr.open_dataset(self.filename, engine='netcdf4')
            keys = ds.variables.keys()
            #
            #
            #
            if incfile is not None:
                ds_inc = xr.open_dataset(incfile, engine='netcdf4')
                inc_keys = ds_inc.variables.keys()
                #
            if incfile is not None:
                ds_azi = xr.open_dataset(azifile, engine='netcdf4')
                azi_keys = ds_azi.variables.keys()

            #
            print(" insar_preprocessing: now updating the downsampled values...")
            #
            for i in range(self.quad_pts.shape[0]):
                 #
                 #
                 if 'x' in keys:
                    da = ds.sel(x=self.quad_pts[i,0],y=self.quad_pts[i,1],method="nearest")
                 else:
                    da = ds.sel(lon=self.quad_pts[i,0],lat=self.quad_pts[i,1],method="nearest")
                 #
                 #if self.mode.upper() == 'LOS':
                 if not self.pot:
                    #
                    # print(self.wavelength / (4*np.pi))
                    if 'z' in keys:
                       self.quad_pts[i,2] = da.z.values * -1 * self.wavelength / (4*np.pi)
                    if 'Band1' in keys:
                       self.quad_pts[i,2] = da.Band1.values * -1 * self.wavelength / (4*np.pi)
                 else:
                    if 'z' in keys:
                        self.quad_pts[i,2] = da.z.values
                    else:
                        self.quad_pts[i,2] = da.Band1.values
                    #
                 #
                 if incfile is not None:
                     #
                     if 'lon' in inc_keys:
                        da = ds_inc.sel(lon=self.quad_pts[i,0],lat=self.quad_pts[i,1],method="nearest")
                     else:
                        da = ds_inc.sel(x=self.quad_pts[i,0],y=self.quad_pts[i,1],method="nearest")
                        #
                     if 'z' in inc_keys:
                        c_inc = da.z.values
                     elif 'Band1' in inc_keys:
                        c_inc = da.Band1.values
                 else:
                     c_inc = self.inc
                     #
                 if azifile is not None:
                     #
                     if 'lon' in azi_keys:
                        da = ds_azi.sel(lon=self.quad_pts[i,0],lat=self.quad_pts[i,1],method="nearest")
                     else:
                        da = ds_azi.sel(x=self.quad_pts[i,0],y=self.quad_pts[i,1],method="nearest")
                        #
                     if 'z' in azi_keys:
                        c_azi = da.z.values
                     elif 'Band1' in azi_keys:
                        c_azi = da.Band1.values
                 else:
                     c_azi = self.azi
                 #
                 nv,ev,uv = qds.losvec(inc=c_inc,azi=c_azi,looking_dir='right',mode=self.mode)
                 #print(c_inc,c_azi)
                 self.quad_pts[i,[3,4,5]] = [ev,nv,uv]
            #
    #
    #
    def export_ds(self,outname,utm=False,zone=None,letter=None):
        #
        #
        if utm:
            import pSAR
            #
            if zone is None or letter is None:
               lons,lats = self.quad_pts[:,0],self.quad_pts[:,1]
               _,_,zone,letter = pSAR.utm_conversion.from_latlon(np.mean(lats), np.mean(lons), force_zone_number=zone, force_zone_letter=letter)
               #
               outname = "%s_%d%s.dat" % (outname,zone,letter)
               print(" insar_preprocessing: rename output as %s" % outname)

            #
            md = self.quad_pts
            md[np.isnan(md[:,2]),2] = 0
            #
            x,y = pSAR.utm_conversion.lltoutm(md[:,1],md[:,0],force_zone_number=zone, force_zone_letter=letter)
            #
            with open(outname+'_utm_%d%s.dat' % (zone,letter),'w') as fid:
                for i in range(md.shape[0]):
                    #
                    if md[i,2]!=0:
                        #
                        fid.write('%f %f %f %f %f %f %f\n' % (x[i]/1000, y[i]/1000, md[i,2],md[i,3],md[i,4],md[i,5],md[i,6]))
                        #
                    #
                #
            #

        #
        #
        fid_xy = open(outname+'.xypoly.gmt','w')
        #
        #
        with open(outname,'w') as fid:
            #
            fid.write('#lon lat los(m) ev ev uv factor(nouse)\n')
            for i in range(self.quad_pts.shape[0]):
                #
                cmll = self.quad_xy[i]
                #
                opoly = quad_xy2poly(cmll)
                carea = calculate_polygon_area(opoly)
                #
                md = self.quad_pts[i,:]
                lon = md[0]
                lat = md[1]
                #
                if utm:
                    x,y,_,_ = pSAR.utm_conversion.from_latlon(lat,lon,force_zone_number=zone, force_zone_letter=letter)
                    #
                    x = x/1000
                    y = y/1000
                else:
                    x = lon
                    y = lat
                #
                if md[2]!=0:
                  #
                  # replace std using area of a single polygon, md[6]
                  # fid.write('%f %f %f %f %f %f %f\n' % (x, y, md[2],md[3],md[4],md[5],md[6]))
                  #
                  fid.write('%f %f %f %f %f %f %f\n' % (x, y, md[2],md[3],md[4],md[5],carea))
                  #
                  #
                  cmll = self.quad_xy[i]
                  #
                  fid_xy.write('>-Z%f\n' % md[2])
                  #
                  # output polygon
                  for j in range(opoly.shape[0]):
                      #
                      #
                      if utm:
                         x,y,_,_ = pSAR.utm_conversion.from_latlon(opoly[j,1],opoly[j,0],force_zone_number=zone, force_zone_letter=letter)
                         fid_xy.write('%f %f\n' % (x/1000,y/1000))
                      else:
                         fid_xy.write('%f %f\n' % (opoly[j,0],opoly[j,1]))
                      #
                  #
              #
        #
        fid_xy.close()
            #
        #
        #
    def quadtree_ds(self,maxblock=512,minblock=4,std=0.5,inc=35,azi=-169, wavelength=None, mode='los',vmin=None,vmax=None,isplot=False,pot=False):
        #
        self.mode     = mode
        self.minblock = minblock
        self.maxblock = maxblock
        self.azi = azi
        self.inc = inc
        self.pot = pot
        #
        if wavelength is not None:
           self.wavelength = wavelength
        #
        if len(self.deramp_data)==0:
            #./insar_preprocessing.py
            print(" Warning: no deramping was performed. We will use original data instead...")
            self.deramp_data = self.data
            #
        #
        #
        outxy,outdata = qds.quad_dsmdata(np.flipud(self.deramp_data),b_block = maxblock,s_block=minblock,std=std)
        mll   = qds.quad_ij2lonlat(outdata,self.lons,self.lats)
        #
        cov_factor = 1
        #
        if mode.upper()=='LOS' or mode.upper() == 'RNG':
           nv,ev,uv = qds.losvec(inc=inc,azi=azi,looking_dir='right',mode='rng')
           #
           cov_factor = -1 * self.wavelength / (4*np.pi)
        elif mode.upper() == 'AZI':
           nv,ev,uv = qds.losvec(inc=inc,azi=azi,looking_dir='right',mode='azi')
        elif mode.upper() == 'E':
            ev = 1
            nv = 0
            uv = 0
        elif mode.upper() == 'N':
            ev = 0
            nv = 1 
            uv = 0
        elif mode.upper() == 'U':
            ev = 0
            nv = 0
            uv = 1
        #
        outpts = np.zeros([len(mll),7])
        #
        for i in range(len(mll)):
            #
            outpts[i,0] = (mll[i][0] + mll[i][1])/2
            outpts[i,1] = (mll[i][2] + mll[i][3])/2
            outpts[i,2] = outdata[i][5] * cov_factor 
            outpts[i,3] = ev
            outpts[i,4] = nv
            outpts[i,5] = uv
            #
        #
        self.quad_pts = outpts
        #
        #
        if isplot:
           #
           extend = [np.min(self.lons),np.max(self.lons),\
                     np.min(self.lats),np.max(self.lats)]
           qds.quad_plot(self.deramp_data.shape,mll,extend=extend,vmin=vmin,vmax=vmax,c_color='seismic')
        #
        self.quad_xy = mll
        #
        #
        flag = 1
        #
        # return mll,outxy,outdata
            
    def deramp(self,order=2,dsc=5,dem=None,npara=3):
        #
        #
        cdata = self.maskdata[::dsc,::dsc]
        #
        if cdata.shape[0]<10:
            print(" Err: only %d points available..." % cdata.shape[0])
            #
        #
        mlons = np.zeros([self.lats.shape[0],self.lons.shape[0]])
        mlats = np.copy(mlons)
        #
        for i in range(self.lats.shape[0]):
            #
            mlons[i,:] = self.lons
        for j in range(self.lons.shape[0]):
            mlats[:,j] = self.lats
        #
        _mlons = mlons[::dsc,::dsc]
        _mlats = mlats[::dsc,::dsc]
        _mlons = _mlons[~np.isnan(cdata)]
        _mlats = _mlats[~np.isnan(cdata)]
        #
        cdata = cdata[~np.isnan(cdata)]
        _mlons = _mlons[cdata!=0]
        _mlats = _mlats[cdata!=0]
        cdata = cdata[cdata!=0]
        #
        if npara == 3:
          if dem is None:
             A = np.zeros([_mlons.shape[0],3])
          else:
             A = np.zeros([_mlons.shape[0],4])
        
          A[:,0] = _mlons
          A[:,1] = _mlats
          A[:,2] = 1
        if npara == 5:
          if dem is None:
             A = np.zeros([_mlons.shape[0],5])
          else:
             A = np.zeros([_mlons.shape[0],6])

          A[:,0] = _mlons**2
          A[:,1] = _mlats**2
          A[:,2] = _mlons
          A[:,3] = _mlats
          A[:,4] = 1

        #
        x,res,rnk,info = np.linalg.lstsq(A,cdata,rcond=None)
        #
        if npara == 3:
           simData = mlons* x[0] + mlats * x[1] + x[2]
        #
        if npara == 5:
           simData = mlons**2 * x[0] + mlats**2 * x[1] + mlons * x[2] + mlats * x[3] + x[4]
        #
        self.deramp_data = self.data - simData
        #
        self.codes[1] = 1
        self.codes[2] = 1
        #
    #
    def mask(self,subext):
        #
        self.subext = subext
        #
        ind_x_1 = self.lons>=subext[0] 
        ind_x_2 = self.lons<=subext[1]
        ind_x = ind_x_1 * ind_x_2
        #
        ind_y_1 = self.lats>=subext[2] 
        ind_y_2 = self.lats<=subext[3]
        ind_y = ind_y_1 * ind_y_2
        #
        #
        #
        index_lon = np.where(ind_x)[0]
        x0 = index_lon[0]
        x1 = index_lon[-1]
        #
        index_lat = np.where(ind_y)[0]
        #
        _y0 = index_lat[0]
        _y1 = index_lat[-1]
        #
        nlat = ind_y.shape[0]
        #
        y1 = nlat-_y0
        y0 = nlat-_y1
        #
        #
        self.maskdata[y0:y1,x0,x1] = np.nan
        #
        self.codes[1] = 1
        #
        #
    def grd_read(self,grdname,ext=None):
        #
        self.filename = grdname
        #
        xrds = xr.open_dataset(self.filename, engine='netcdf4')
        #
        keys = xrds.variables.keys()
        #
        # print(keys)
        #
        xlabel = 'lon'
        if 'x' in keys:
            xlabel = 'x'
        #
        ylabel = 'lat'
        if 'y' in keys:
            ylabel = 'y'
        #
        zlabel = 'z'
        if 'Band1' in keys:
            zlabel = 'Band1'
        #
        #
        if ext is None:
           lons = xrds[xlabel].values
           lats = xrds[ylabel].values
           data = xrds[zlabel].values
           #
        else:
           da = xrds[zlabel]
           subz = da.where((da[xlabel]>=ext[0]) & \
                      (da[xlabel]<=ext[1]) & \
                      (da[ylabel]>=ext[2]) & \
                      (da[ylabel]<=ext[3]), drop=True)
           #
           lons = subz[xlabel].values
           lats = subz[ylabel].values
           data = subz.values
        #
        data = np.flipud(data)
        #
        #
        self.maskdata = np.copy(data)
        self.lons = lons
        self.lats = lats
        self.data = data
        self.gext = [np.min(lons[~np.isnan(lons)]),np.max(lons[~np.isnan(lons)]),\
                     np.min(lats[~np.isnan(lats)]),np.max(lats[~np.isnan(lats)])]
        #
        self.codes[0] = 1
        #
    #
    #
    
#
if __name__ == '__main__':
    #
    # 
    # a simple log
    if len(sys.argv)<2:
        this_help()
        #
    #
    logname = os.path.basename(sys.argv[0])+'.log'
    #
    log(sys.argv,logname=logname)
    #
    #
    #
    km = 0 
    dsc = 10
    wavelength = 0.0556
    inc = 35
    azi = -10
    vmin = -0.25
    vmax = 0.25
    isplot = False
    #
    # an example insar data
    #
    ingrd      = 'test_strike_right_rng.grd'
    maskext    = None #[99.75, 100.3, 29.75, 30.25]
    maxblock   = 256
    minblock   = 16
    wavelength = 0.0556
    std        = 0.75
    outroot    = None
    un = None
    ul = None
    incfile = '' 
    azifile = ''
    mode='LOS'
    gohelp = False
    pot = False
    #
    #
    #
    if len(sys.argv)>1:
        ingrd = sys.argv[1]
        #
    #
    for i,key in enumerate(sys.argv):
        #
        if key == '-pot':
            pot = True
        if key == '-help':
            gohelp = True
        if key == '-incfile':
            incfile = sys.argv[i+1]
        if key == '-azifile':
            azifile = sys.argv[i+1]
        if key == '-un':
            un = int(sys.argv[i+1])
            #
        if key == '-ul':
            ul = sys.argv[i+1]
            #
        #
        if key == 'km':
            km = 1
        if key == '-plot':
            isplot = True
        if key == '-outroot':
            outroot = sys.argv[i+1]
        if key == '-inc': 
            inc = float(sys.argv[i+1])
        if key == '-azi':
            azi = float(sys.argv[i+1])
        if key == '-mode':
            mode = sys.argv[i+1]
            #
        if key == '-wavelength':
            wavelength = float(sys.argv[i+1])
            #
        if key == '-maxblock':
            maxblock = int(sys.argv[i+1])
        if key == '-minblock':
            minblock = int(sys.argv[i+1])
            #
        if key == '-std':
            std = float(sys.argv[i+1])
            #
        if key == '-maskext':
            maskext = [float(cext) for cext in sys.argv[i+1].split(',')]
    #
    if gohelp:
        this_help()
        #
    #
    if outroot is None:
        #
        outroot = '%s_quadds' % os.path.basename(ingrd).split('.grd')[0]
    #
    myinsar = InSAR_preprocessing(km=km,wavelength=wavelength)
    #
    myinsar.grd_read(ingrd)
    #
    if un is None:
        #
        meanlat = np.mean(myinsar.lats)
        meanlon = np.mean(myinsar.lons)
        #
        utminfo = get_utm_info(meanlat, meanlon)
        #
        myinsar.un = utminfo['un']
        myinsar.ul = utminfo['ul']
        #
        #
    #
    if isplot:
      myinsar.insarPlot()
    #
    #
    if maskext is not None:
      myinsar.mask(maskext)
      #
      myinsar.deramp(dsc=dsc)
      #
      if isplot:
         myinsar.insarPlot()
         #
      #
      #
      print(" *** subsetting data with a region of -R%f/%f/%f/%f" % (maskext[0],maskext[1],maskext[2],maskext[3]))
      #
      myinsar.subset(ext=maskext)
      myinsar.export_grd(outroot)
    #
    if isplot:
        myinsar.insarPlot()
    #
    #  
    myinsar.quadtree_ds(maxblock=maxblock, wavelength=wavelength, minblock=minblock, inc=inc,azi=azi, \
            mode=mode,std=std,vmin=vmin,vmax=vmax,isplot=isplot,pot=pot)
    #
    #
    if os.path.exists(incfile) and os.path.exists(azifile):
       myinsar.update_pts(incfile=incfile,azifile=azifile)
    #
    #
    if km == 1:
        myinsar.export_ds(outroot+'_km_UTM_%d%s.dat' % (myinsar.un,myinsar.ul))
    else:
        myinsar.export_ds(outroot+'_ll.dat')
        myinsar.export_ds(outroot,utm=True,zone=un,letter=ul)
    #
    #
#


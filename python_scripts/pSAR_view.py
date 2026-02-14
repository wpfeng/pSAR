#!/usr/bin/env python
#
#
####################################
import pSAR,os,sys
#import pSIMP
import numpy as np
#import pisce 
# 
# below are for geotif reader
#
import pDATA
import glob
from matplotlib import cm
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
#
def get_phs(argv):
    #
    Flag = True
    phss = []
    goflag = True
    for ni in range(1,len(argv)):
        cargv = argv[ni]
        #
        # refine by Wanpeng Feng, @SYSU, 2018-08-09
        # 
        if os.path.exists(cargv) and goflag:
           phss.append(cargv)
        else:
           goflag = False
           
           #if '-' in cargv:
           #Flag = False
    #
    return phss
####################################
if len(sys.argv) < 2:
    helpstr=\
    '''INFO    : @(#)RCM InSAR Automatic Processing System 
    
    %s <in_1> <in_2>... -minv [-5] -maxv [5] -wrap -scale 1 -r [1] -l [1]
                 -gamma [not given in default] -dims [None,None in default] -disp 
                 -dtype [NULL in default, can be float32, Int16, float64] -title [None]
                 -cmap [a colormap name, e.g. RdYlBu] -pts [None] -zerov [None in default]
                 -complex -fix3 -ccmask [None] -ccthresh [0.25 in default]
                 -utm [not given in default] -geotifdir [None in default]
                 -floor [not given in default]
                 -ext [None in default]
                 -interp [none in default]
                 -prof [none in default, profile line]
    +Input:
    
    <in_phs>   a binary phase file with a rsc header or a GAMMA-format file
    -scalar    amplify factor, None in default
    -wrap      flag for rewrapping phs
    -minv     -5 in default
    -maxv      5 in default
    -geotifdir a path to search geotif files, None in default
    -r         number of row  for subplot
    -l         number of line for subplot
    -dims      None,None in default
    -gamma     flag for checking if it is a GAMMA format (big-endien)
    -bands     1 in default
    -zerov     fill zero pixels with any given value
    -complex   optional, used to read GAMMA complex format data. return phase (atan2(imag,pwr)) only
    -fix3      optional, if given, only first two are used to generate difference between them 
    -ccmask    a cc file to mask pixels with a threshold 
    -ccthresh  0.25 in default
    -utm       keyword for showing data in utm projection with coordinate units (km).
    -floor     if given, for example, 3.14159265, the largest number after data / scalar(given by floor)
    -outfig    [None in default]
               rootname of outfile, if given, <outfig>.png and <outfig>.pdf will be produced...
    -cycle     False in default. If given, number of cycles will be return...
    -shift     False in default. If given, info['Z_SHIFT'] will be considered during reading the data...
    -grdxy     False in default, If given, ext of the data will be reset based on number of rows and lines
    -aoff      flag to remove annotations for x and y axises
    -ext       extension for plot, None in default
    -pts_label False in default. When it is given, pts_label will be set as True.
    -loop      False in default
    -simp      [None in default]
    -mod       [None in default] 
    -isce     [could be isce_int, isce_cor, isce_los...]
    #
    (c) 2015-2018, Canada Centre for Mapping and Earth Observation, Natural Resources Canada (NRCan), Canada
   
    Development History:
    -----------------------------------------------------------
    First version was done by Feng, Wanpeng, @NRCan, 2016-03-18
    -----------------------------------------------------------
    A file in GAMMA format is supported now.
    Updated by Wanpeng Feng, @NRCan, 2017-02-21. 
    Updated by Wanpeng Feng, @SYSU, 2018-08-10
    Updated by Wanpeng Feng, @SYSU, 2019-03-15
      since this version, it is allowed to search geotif folders to add geotif polygons
    Updated by Wanpeng Feng, @SYSU, 2019/07/19
      since this version, a file in geotiff format can be suppored by pSAR_view.py
    Updated by Wanpeng Feng, @SYSU, 2020/09/13
      since this version, grdxy is available
    Updated by Wanpeng FEng, @SYSU, Zhuhai, 2020/10/22
      since this version, aff is avaiable. When it is set, x and y axises will be removed....
    '''
    print(helpstr % sys.argv[0]) 
    sys.exit(1)
#
loop = False
profile = None
show_ext = None
outfig = None 
band   = [1]
l      = '1'
r      = '1'
scale  = 1
isce = None
mod = False 
interpolation_method = 'none'
minv   = -3.1415926
maxv   = 3.1415926
ccthresh = 0.25
ccmask   = 'none'
iswrap = False
dtype  = "NULL"
disp   = False
dims   = [None,None]
isgamma = False
geotifdir = None
cmap  = 'RdYlBu'
cmap = 'jet'
cycle = False
pts   = None
zerov = np.nan 
title = None
iscomplex = False
fix3      = False
utm       = False
isfloor   = False
floor_scalar = np.pi
scalar = None
shift = False
grdxy = False
aoff  = False
pts_label = False
simp      = None
####################################
# Read keywords from inputs
#
k = 0
ARGV = sys.argv
for c_keyword in ARGV:
    #
    k += 1
    #
    if c_keyword == '-isce':
       isce = sys.argv[k]
    #
    if c_keyword == '-mod':
       mod = True 
    if c_keyword == '-simp':
        simp = sys.argv[k]
        #
    if c_keyword == '-loop':
        loop = True
    if c_keyword == '-pts_label':
       pts_label = True
    if c_keyword == '-prof' or c_keyword == '-profile':
        profile = sys.argv[k]
    if c_keyword == '-aoff':
        aoff = True
    if c_keyword == '-interp':
        interpolation_method = sys.argv[k]
    if c_keyword == '-w':
        iswrap = True
    if c_keyword == '-ext':
        show_ext = ARGV[k].split(',')
        show_ext = [float(ci) for ci in show_ext]
    if c_keyword == '-grdxy':
        grdxy = True
    if c_keyword == '-shift':
       shift = True
    if c_keyword == '-cycle':
       cycle = True
       isfloor = True
    if c_keyword == '-outfig':
       outfig = sys.argv[k]
    if c_keyword == '-scalar':
       scalar = float(sys.argv[k])
    if c_keyword == '-floor':
       isfloor = True
       floor_scalar = float(sys.argv[k])
    if c_keyword == '-geotifdir':
       geotifdir = sys.argv[k]
    if c_keyword == '-utm':
       utm = True
    if c_keyword == '-ccthresh':
       ccthresh = float(sys.argv[k])
    if c_keyword == '-ccmask':
       ccmask = sys.argv[k]
    if c_keyword == '-fix3':
       fix3 = True
    if c_keyword == '-complex':
       iscomplex = True 
    if c_keyword == '-title':
       title = sys.argv[k].split(',')
    if c_keyword == '-bands':
       band = sys.argv[k].split(',')
       band = [int(cband) for cband in band]
    if c_keyword == '-zerov':
       zerov = float(sys.argv[k])
    if c_keyword == '-pts':
       pts = sys.argv[k]
    if c_keyword == '-cmap':
       cmap = sys.argv[k]
    if c_keyword == "-gamma":
       isgamma = True
    if c_keyword == "-dims":
       dims = sys.argv[k].split(',')
       dims = [int(float(dims[0])),int(float(dims[1]))]
    if c_keyword == '-disp':
       disp = True
    if c_keyword == "-dtype":
       dtype = sys.argv[k]
    if c_keyword == "-minv":
       if len(sys.argv) > k:
          if sys.argv[k]=='None':
             minv = None
          else:
             minv = float(sys.argv[k])
    #
    if c_keyword == "-maxv":
       if len(sys.argv) > k:
          maxv = float(sys.argv[k])
    #
    if c_keyword == "-wrap":
       iswrap=True
    if c_keyword == "-scale":
       if len(sys.argv) > k:
          scale = int(sys.argv[k])
    if c_keyword == "-l": 
       if len(sys.argv) > k: 
          l = sys.argv[k]
    if c_keyword == "-r":
       if len(sys.argv) > k:
          r = sys.argv[k]
####################################
# logging
pSAR.util.log(sys.argv)
#
in_file       = sys.argv[1]
in_rsc        = in_file + '.rsc'
#
if isgamma:
   dtype = ">f4"
####################################
#
counter = 0
#
internalflag = False
#
if minv is None:
   internalflag=True
#
phsfiles = get_phs(sys.argv[:])
#
print("Bands: ",band)
print("Number of inputs: %d" % len(phsfiles))
#
if len(band) < len(phsfiles):
   #outband = np.zeros(len(phsfiles))
   #
   band = [band[0] for i in range(len(phsfiles))]
#
i = 0
if fix3:
   r = 3
if profile is not None:
    prof_data = np.loadtxt(profile)
    #
    #
#
if simp is not None:
    #
    fpara,un,ul = pSIMP.import_simp(simp)
    #
    llpolys,deps = pSIMP.simp_fpara2geopoly(fpara,un,ul)
    #
for cphs in get_phs(sys.argv[:]):
    i += 1
    if os.path.isfile(cphs):
       #
       counter += 1
       plt.subplot(int(l),int(r),counter)
       #
       bname = os.path.basename(cphs)
       #
       fext = bname.split('.')[-1]
       #
       #
       if isce is not None:
           dtype = 'float32'
           #
       #
       #
       if dtype.upper() == "NULL":
          #
          if fext.upper() == 'GRD' or fext.upper() == 'NC':
             fmt = " based on GRD data itself"
          elif fext.upper() == 'TIF' or fext.upper() == 'TIFF':
             fmt = " based on GTiff data itself"
          else:
             fmt  = pSAR.roipac.roi_to_fmt(cphs,band=band[i-1])
       else:
          fmt  = dtype
       #
       # 
       print(' %s in %s' % (cphs,fmt))
       if fext.upper()=='GRD' or fext.upper()=='NC':
          cinfo,ext,data = pDATA.grd_read(cphs)
          #
          if grdxy:
              #
              ext = [0,data.shape[1],0,data.shape[0]]
          #
          fmt = data.dtype.name
          print(" fmt: %s" % fmt)
          #
       elif fext.upper() == 'TIF' or fext.upper() == 'TIFF': 
          info,ext = pDATA.geotif2info(cphs)
          data,ext,proj = pDATA.read_geotif(cphs,noband=1)
          #
       else:
          if isce is None:
             data = pSAR.roipac.roi_read(cphs,band=band[i-1],iscomplex=iscomplex,dtype=fmt,dims=dims)
             info,ext = pSAR.roipac.rsc_read(cphs+'.rsc')
             if info['Z_SHIFT'] != 'NONE' and shift:
                print(" Using a shift of %s" % info['Z_SHIFT'])
                data[data!=0] = data[data!=0] - float(info['Z_SHIFT'])
                #
          else:
             isce_mode = isce.split('_')[1]
             bnd = isce.split('_')[2]
             data,info, ext = pisce.read_isce_data(cphs,bnd=int(bnd),mode=isce_mode)
             #
       #
       #
       data[np.isnan(data)] = 0.
       data[data==np.inf]   = 0.
       if fix3:
          if i == 1:
             data1 = np.copy(data)
          if i == 2:
             data2 = np.copy(data)
       #
       crsc = cphs + '.rsc'
       if os.path.exists(crsc):
          cinfo, cext = pSAR.roipac.rsc_read(crsc)
          if grdxy:
              cext = ext
          #
          if utm:
             #
             # added by Wanpeng Feng, @SYSU, 2018-08-10
             #
             ux,uy,zoneN,zoneL = pSAR.utm_conversion.from_latlon((cext[2]+cext[3])/2.,(cext[0]+cext[1])/2.)
             #
             ux1,uy1,zoneN,zoneL = pSAR.utm_conversion.from_latlon(cext[2],cext[0],\
                     force_zone_number=zoneN,force_zone_letter=zoneL)
             ux2,uy2,zoneN,zoneL = pSAR.utm_conversion.from_latlon(cext[3],cext[1],\
                     force_zone_number=zoneN,force_zone_letter=zoneL)
             cext = [ux1/1000.,ux2/1000.,uy1/1000.,uy2/1000.]
             #
             print(' UTM projection info: ZONE %d & %s' % (zoneN,zoneL))
             
       else:
           if fext.upper()=='GRD':
               cext = ext
           elif fext.upper() == 'TIF' or fext.upper() == 'TIFF':
               cext = ext
           else:
               if isce is None:
                  cext = [0,data.shape[1]-1,0,data.shape[0]-1]
               else:
                  cext = ext
       #
       if fmt.upper() == "INT16":
          data = data.astype(float)
       #
       data = data[::scale,::scale]
       #
       if scalar is not None:
           data = data * scalar
       #
       if disp:
          wavelength = float(cinfo['WAVELENGTH'])
          data = data * -1 * wavelength * 100 / (4 * np.pi)
       if internalflag:
          #
          cdata= data[data!=0]
          if cdata.shape[0] == 0:
             print(" pSAR_view: ERROR!!! All input data are zero .")
             sys.exit(-1) 
          #
          minv = np.min(cdata)
          maxv = np.max(cdata)
       #
       if iswrap and not isfloor:
          data = pSAR.roipac.wrap(data,minv=minv,maxv=maxv,isplot=0,scale=scale)
          #
       elif isfloor:
          #
          if loop:
            #
            print(' calculating the cycle number of the phase colosure...')
            data = pSAR.roipac.loop2cycles(data)
          else:
            data = pSAR.roipac.phs2cycles(data,minv=-np.pi,maxv=np.pi,scale=scale)
       #
       #
       if zerov is not None:
           data = data.astype('float')
           data[data == 0] = np.nan
       #
       maskdata = data * 0. + 1.
       if ccmask.upper() != 'NONE':
          maskdata = pSAR.roipac.roi_read(ccmask)
          maskdata = maskdata[::scale,::scale]
          maskdata[maskdata<ccthresh] = 0.
       #
       #
       data = data * maskdata
       #
       if cext[0]>180: 
           cext[0] = cext[0]-360
       if cext[1]>180:
           cext[1] = cext[1]-360
       #
       #if show_ext is not None:
       #    cext = show_ext;
       #
       # print(cext)
       if mod:
           data = np.mod(data+np.pi,np.pi*2)-np.pi
           #
       #
       im = plt.imshow(data,extent=cext,interpolation=interpolation_method,vmin=minv,vmax=maxv,cmap=plt.get_cmap(cmap))
       #
       if geotifdir is not None:
         tifs = glob.glob(geotifdir+'/*/*.tif')
         #
         print(" ** we found %d TIFs in %s..." % (len(tifs),geotifdir))
         #
         for ctif in tifs:
           #
           try:
             exts = pDATA.geotif_extent(ctif)
             exts = np.array(exts)
             cpoly = np.zeros([5,2])
             for i in range(4):
               cpoly[i,:] = exts[i,:]
             cpoly[4,:] = cpoly[0,:]
             #
             # print(" *** now loading %s into the figure..." % os.path.basename(ctif))
             plt.plot(cpoly[:,0],cpoly[:,1],linestyle='solid',color='red',linewidth=5.0)
             plt.text(np.mean(cpoly[:,0]),np.mean(cpoly[:,1]),os.path.basename(ctif))
           except:
             print(" ERROR:  %s is corrupted!!! ***" % ctif)
         #
       if show_ext is not None:
         axid = plt.gca()
         axid.set_xlim([show_ext[0],show_ext[1]])
         axid.set_ylim([show_ext[2],show_ext[3]])

       if aoff:
           plt.axis('off')
           #
       if pts is not None:
          ptdata = np.loadtxt(pts)
          dims = len(ptdata.shape)
          if dims == 1:
              t = np.zeros([1,2])
              t[0,:] = ptdata
              ptdata = t 
              
              #
          #
          print(ptdata.shape)
          # ptdata[ptdata[:,0]>180,0] = ptdata[ptdata[:,0]>180,0] - 360
          print(ptdata)
          #
          plt.plot(ptdata[:,0],ptdata[:,1],'vr',ms=15)
          for pt_i in range(ptdata.shape[0]):
              #
              if pts_label:
                 plt.text(ptdata[pt_i,0],ptdata[pt_i,1],'No-%d' % pt_i)
              #
       if profile is not None:
           #
           plt.plot(prof_data[:,0],prof_data[:,1],'-',color='yellow',linewidth=3.5,zorder=20)
           plt.plot(prof_data[0,0],prof_data[0,1],'o',color='yellow',markersize=4.5,zorder=20)
           plt.plot(prof_data[-1,0],prof_data[-1,1],'o',color='yellow',markersize=4.5,zorder=20)
       if title is not None:
          ctitle = title[i-1]
       else:
          ctitle = os.path.basename(cphs).split('.')[0] #bname
       #
       if simp is not None:
         #
         for i in range(llpolys.shape[0]):
           cpolys = np.reshape(llpolys[i,:],[2,5])
           plt.plot(cpolys[1,:],cpolys[0,:],'-r',linewidth=3.5,label='Fault')
       #
       if 'geo_' in ctitle:
           ctitle = ctitle.split('geo_')[1]
           #
       #
       plt.title(ctitle)
       ax = plt.gca()
       divider = make_axes_locatable(ax)
       cax = divider.append_axes("right", size="5%", pad=0.05)
       cbar = plt.colorbar(im, cax=cax)
       #
       # since matplitlib (v3.1), to replace the below with the after
       # modified by Wanpeng Feng, @SYSU, Guangzhou, 2019/06/13
       #
       # cbar.set_clim(minv,maxv)
       cbar.mappable.set_clim(minv,maxv)
#
if fix3:
   #
   plt.subplot(int(l),int(r),3)
   ddata = data1-data2
   if iswrap:
      ddata = pSAR.roipac.wrap(ddata,minv=minv,maxv=maxv,isplot=0,scale=1)
       #
   if zerov is not None:
      data[data == 0] = np.nan
       #
   # minv = -0.00005
   # maxv =  0.00005
   maskdata = dddata * 0. + 1.
   if ccmask.upper() != 'NONE':
       maskdata = pSAR.roipac.roi_read(ccmask)
       maskdata = maskdata[::scale,::scale]
       maskdata[maskdata<ccthresh] = 0.
   #
   ddata = ddata * maskdata
   im = plt.imshow(ddata,extent=cext,vmin=minv,vmax=maxv,cmap=plt.get_cmap(cmap))
   plt.axis('off')
   plt.title("DIff")
   ax = plt.gca()
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(im, cax=cax)
   cbar.set_clim(minv,maxv)
   #
if outfig is None:
   #
   plt.show()
else:
   plt.tight_layout()
   fig = plt.gcf()
   fig.set_size_inches(10,6, forward=True)
   #
   plt.savefig(outfig+'.png',dpi=720)
   plt.savefig(outfig+'.pdf',dpi=720)
       

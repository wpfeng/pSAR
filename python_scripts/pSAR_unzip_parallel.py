#!/usr/bin/env python
#
#
#
import glob
import os
import sys
import multiprocessing as mpg
import pDATA
#
#
def go_unzip(i):
    #
    czip = zips[i]
    if ischeck != 0:
       flag = pDATA.zipcheck(czip)
    else:
       flag = True
    if flag:
        #
        # added by Wanpeng Feng to avoid additional unzipping
        # 2024/02/18
        #
        go_unzip = True
        #
        if rawcheck:
            #
            cdate = os.path.basename(czip).split('_')[5][0:8]
            raw_dir = glob.glob('../../raw/S*%s*.SAFE/' % cdate)
            if len(raw_dir)>0:
                go_unzip = False
                #
            #
        if go_unzip:
          goStr = 'ginsar_unzip.sh %s %s' % (objdir,czip)
          print(goStr)
          os.system(goStr)
        return ''
    else:
        return czip 
#
######################################################
#
global zips,objdir,ischeck,rawcheck
#
#######################################################
if len(sys.argv)<3:
    helpstr = \
        '''
        %s <objdir> <njob> [completeness_check, 0 in default] -rawcheck
        ++++++++++++++++++++++++++++++++++++++++++
        To unzip zip files in parallel...
        #
        <objdir> a path to save unzipped folders.
        <njob>   number of workers to work at the same time...
        <completeness check> 1 yes, 0 ignored
        Written by Wanpeng Feng,@SYSU, Zhuhai, 2021/01/07

        '''
    print(helpstr % sys.argv[0])
    sys.exit(-1)
#
#
#######################################################
#
rawcheck = False
#
objdir   = sys.argv[1]
njob     = int(sys.argv[2])
ischeck  = 0
#
for i,key in enumerate(sys.argv):
    #
    if key == '-rawcheck':
        #
        rawcheck = True
#
if len(sys.argv)>3:
    ischeck = int(sys.argv[3])
#
if not os.path.exists(objdir):
    os.makedirs(objdir)
#
if True:
    zips = glob.glob('S*.zip')
    #
    allindex = [i for i in range(len(zips))]
    #
    pool = mpg.Pool(processes=njob)
    result = pool.imap(go_unzip,allindex)
    pool.close()
    pool.join()
    #
    with open('badzip.list','w') as fid:
        for cline in result:
            if len(cline)>0:
                fid.write('%s\n' % cline)
            #
        #
    #

#!/usr/bin/env python
#
#
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import Normalize
#
#
def losvec(inc,azi,mode='rng', looking_dir='right'):
    #
    # inc and azi in degree
    #
    inc_rad = inc * np.pi / 180
    azi_rad = azi * np.pi / 180
    #
    if looking_dir=="LEFT":
       azi_rad = azi * math.pi / 180 + 3.14159265
    #
    if mode.upper() == 'RNG':
       nv =      np.sin(azi_rad) * np.sin(inc_rad)
       ev = -1 * np.cos(azi_rad) * np.sin(inc_rad)
       uv =      np.cos(inc_rad)
       #
    elif mode.upper() == 'AZI':
       nv = np.cos(azi_rad)
       ev = np.sin(azi_rad)
       uv = 0
       #
    #
    return nv,ev,uv
#
def quad_ij2lonlat(corners,lons,lats):
    #
    #
    mll = []
    #
    for i in range(len(corners)):
        #
        c4p = corners[i]
        #
        if len(c4p)==4:
           cll =[lons[c4p[0]],lons[c4p[1]],lats[c4p[2]],lats[c4p[3]],0,0]
        elif len(c4p)==5:
           cll =[lons[c4p[0]],lons[c4p[1]],lats[c4p[2]],lats[c4p[3]],c4p[4],0]
        else:
           cll =[lons[c4p[0]],lons[c4p[1]],lats[c4p[2]],lats[c4p[3]],c4p[4],c4p[5]]
        mll.append(cll)
        #
    #
    return mll
#   
#
def quad_block2polygon(outxy):
    #
    corners = []
    values = np.zeros([len(outxy),2])
    #
    for i in range(len(outxy)):
        #
        cxy = outxy[i]
        cpoly = [[outxy[i][0],outxy[i][2]],\
                 [outxy[i][0],outxy[i][3]],\
                 [outxy[i][1],outxy[i][3]],\
                 [outxy[i][1],outxy[i][2]],\
                 [outxy[i][0],outxy[i][2]]]
        #
        if len(cxy)>4:
            values[i,0] = cxy[4]
        if len(cxy)>5:
            values[i,1] = cxy[5]
            #
        #
        corners.append(np.array(cpoly))
        #
    #
    return corners, values
#
def quad_blocksValues(data,blocks):
    #
    outdata = []
    for i in range(len(blocks)):
        #
        blockxy = blocks[i]
        cdata = data[blockxy[2]:blockxy[3],blockxy[0]:blockxy[1]]
        #
        tdata = cdata[~np.isnan(cdata)]
        tdata = tdata[tdata!=0]
        #
        #
        if len(tdata>0):
            outdata.append([blockxy[0],blockxy[1],blockxy[2],blockxy[3],blockxy[4],np.mean(tdata)])
            #
        #
    #
    return outdata
#
def quad_stdINblock(data,blockxy):
    #
    cdata = data[blockxy[2]:blockxy[3],blockxy[0]:blockxy[1]]
    #
    tdata = cdata[~np.isnan(cdata)]
    tdata = tdata[tdata!=0]
    #
    if len(tdata)>0:
        #
        return np.std(tdata)
        #
    else:
        return 0
    #
def quad_dsmdata(data,b_block=512,s_block=32,std=0.4):
    #
    dims = data.shape
    #
    init_blocks = quad_data2blocks(dims,xdim_c=None,ydim_c=None,block_size=b_block,minblock=s_block*4)
    #
    undecide_blocks = init_blocks 
    #
    #
    decided_blocks = []
    tmp_undecide_blocks = []
    #
    block = b_block
    #
    #
    while len(undecide_blocks)>0 and block>=s_block:
        #
        block = int(block / 2)
        #
        tmp_undecide_blocks = []
        #
        for ii in range(len(undecide_blocks)):
            #
            blockxy = undecide_blocks[ii]
            #
            cstd = quad_stdINblock(data,blockxy)
            #
            #
            if cstd <= std:
                #
                decided_blocks.append(blockxy)
            else:
                #
                #
                refined_blocks = quad_data2blocks(dims,xdim_c=[blockxy[0],blockxy[1]],ydim_c=[blockxy[2],blockxy[3]],block_size=block)
                #
                # print(len(refined_blocks),blockxy[1]-blockxy[0],block)
                #
                for jj in range(len(refined_blocks)):
                    #
                    if block > s_block:
                      #
                      tmp_undecide_blocks.append(refined_blocks[jj])
                      #
                    else:
                      #
                      decided_blocks.append(refined_blocks[jj])
                    #
                #
            #
        #
        print( "   +++ %d blocks left under STD of %f with a searching window of %d * %d" % (len(decided_blocks),std, block,block))
        #
        undecide_blocks = tmp_undecide_blocks
        #
    #
    outdata = quad_blocksValues(data,decided_blocks)
    return decided_blocks,outdata
#
def quad_plot(dims,outxy,extend=None,c_color='jet',vmin=None,vmax=None):
    #
    corners,values = quad_block2polygon(outxy)
    #
    cmap = plt.get_cmap(c_color)
    # 
    # setting of color ranges
    #
    if vmin is None:
        vmin = np.min(values[:,1])
    if vmax is None:
        vmax = np.max(values[:,1])
        #
    #
    norm_c = Normalize(vmin=vmin, vmax=vmax)
    #
    fig, ax = plt.subplots()
    #
    for i in range(len(corners)):
        #
        pxy = corners[i]
        #
        plt.plot(pxy[:,0],pxy[:,1],color='white',linewidth=0.5)
        #
        #
        polygon = plt.Polygon(list(pxy), closed=True)
        #
        color = cmap(norm_c(values[i,1]))
        polygon.set_facecolor(color)
        #
        ax.add_patch(polygon)
    #
    ax.axis('equal')
    #
    if extend is None:
      ax.set_ylim([0,dims[0]-1])
      ax.set_xlim([0,dims[1]-1])
    else:
      ax.set_ylim([extend[2],extend[3]])
      ax.set_xlim([extend[0],extend[1]])
    #
    synvalue = np.linspace(vmin,vmax,num=100)
    #
    m = cm.ScalarMappable(cmap=c_color)
    m.set_array(synvalue)
    # cbar = plt.colorbar(m,orientation='vertical')
    #
    plt.show()
#
def quad_data2blocks(dims,xdim_c=None,ydim_c=None,block_size=512,minblock=32):
    #
    #
    # 
    # dim  = data.shape
    #
    if xdim_c is None:
        #
        xdim_c = [0,dims[1]]
        #
    if ydim_c is None:
        #
        ydim_c = [0,dims[0]]
        #
    #
    outblocks = []
    #
    nb_x = int((xdim_c[1] - xdim_c[0])/block_size)
    nb_y = int((ydim_c[1] - ydim_c[0])/block_size)
    #
    if (xdim_c[0] + nb_x*block_size) < xdim_c[1] :
        nb_x = nb_x + 1
        #
    if (ydim_c[0] + nb_y*block_size) < ydim_c[1]:
        nb_y = nb_y + 1
    #
    #
    for i in range(nb_x):
        #
        for j in range(nb_y):
            #
            scale = 2
            #
            outxy = [i*block_size + xdim_c[0], (i+1)*block_size + xdim_c[0], j*block_size + ydim_c[0], (j+1)*block_size + ydim_c[0],(block_size)]
            #
            if outxy[1] <= (dims[1]-1) and outxy[3] <= (dims[0]-1): #int(block_size/scale)<=minblock:
               #
               outblocks.append(outxy)
               #
            else:
               if int(block_size/scale)>=minblock:
                  #
                  coutxy = quad_data2blocks(dims,xdim_c=[outxy[0],outxy[1]],ydim_c=[outxy[2],outxy[3]],block_size = int(block_size/scale),minblock=minblock)
                  #
                  scale *= 2
                  #
                  #
                  for j in range(len(coutxy)):
                     #
                     if coutxy[j][1] <= (dims[1]-1) and coutxy[j][3] <= (dims[0]-1):
                        outblocks.append(coutxy[j])
                     #

            #
        #
    #
    return outblocks
#
#
if __name__ == '__main__':
    #
    #
    dsc = 10
    ingrd   = '2024055_2024067_ll.grd'
    maskext = [92.9, 93.2, 33.41, 33.7]
    #
    myinsar = InSAR_preprocessing()
    #
    myinsar.grd_read(ingrd)
    #
    myinsar.mask(maskext)
    #
    myinsar.deramp(dsc=dsc)
    #
    data = myinsar.deramp_data
    #
    # dims = [2350,3725]
    #
    block_size =512
    #
    #
    quad_dsmdata(data,b_block=512,s_block=64)
    #
    sys.exit(-1)
    #
    outxy = quad_data2blocks(dims,xdim_c=None,ydim_c=None,block_size=block_size) 
    #
    corners = quad_block2polygon(outxy)
    #
    for i in range(len(corners)):
        #
        plt.plot(corners[i][:,0],corners[i][:,1],'-b')
        #
    #
    outxy1 = quad_data2blocks(dims,xdim_c=[outxy[2][0],outxy[2][1]],ydim_c=[outxy[2][2],outxy[2][3]],block_size=block_size/2)
    corners1 = quad_block2polygon(outxy1)
    #
    for i in range(len(corners1)):
        #
        plt.plot(corners1[i][:,0],corners1[i][:,1],'-r')
    #
    ax = plt.gca()
    #
    ax.set_ylim([0,dims[0]-1])
    ax.set_xlim([0,dims[1]-1])
    #

    plt.show()

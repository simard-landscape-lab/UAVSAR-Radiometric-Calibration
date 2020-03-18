# -*- coding: utf-8 -*-
"""
Created on Thu Nov 10 15:21:17 2016

@author: mdenbina
"""
import numpy as np
import matplotlib.pyplot as plt


lutpath = '/Users/mdenbina/src/radiocal/vegetation_lut/'

hhlut = np.memmap(lutpath+'caltbl_Louisiana_GulfCoast_Wetlands_InclForest_Oct2016_HH.flt', shape=(900,900), dtype='<f4', mode='c')
hvlut = np.memmap(lutpath+'caltbl_Louisiana_GulfCoast_Wetlands_InclForest_Oct2016_HV.flt', shape=(900,900), dtype='<f4', mode='c')
vvlut = np.memmap(lutpath+'caltbl_Louisiana_GulfCoast_Wetlands_InclForest_Oct2016_VV.flt', shape=(900,900), dtype='<f4', mode='c')

plt.figure()
plt.plot(np.nanmean(hhlut,axis=0),label='HH')
plt.plot(np.nanmean(hvlut,axis=0),label='HV')
plt.plot(np.nanmean(vvlut,axis=0),label='VV')



look_low = 22
look_low_bin = int(np.floor(look_low*10))

look_high = 65
look_high_bin = int(np.floor(look_high*10))

hhlut[:,0:look_low_bin] = hhlut[:,look_low_bin,np.newaxis]
hvlut[:,0:look_low_bin] = hvlut[:,look_low_bin,np.newaxis]
vvlut[:,0:look_low_bin] = vvlut[:,look_low_bin,np.newaxis]

hhlut[:,look_high_bin+1:] = hhlut[:,look_high_bin,np.newaxis]
hvlut[:,look_high_bin+1:] = hvlut[:,look_high_bin,np.newaxis]
vvlut[:,look_high_bin+1:] = vvlut[:,look_high_bin,np.newaxis]

plt.figure()
plt.plot(np.nanmean(hhlut,axis=0),label='HH')
plt.plot(np.nanmean(hvlut,axis=0),label='HV')
plt.plot(np.nanmean(vvlut,axis=0),label='VV')


hhlut.tofile(lutpath+'caltbl_Louisiana_GulfCoast_Wetlands_InclForest_Oct2016_v2_HH.flt')
hvlut.tofile(lutpath+'caltbl_Louisiana_GulfCoast_Wetlands_InclForest_Oct2016_v2_HV.flt')
vvlut.tofile(lutpath+'caltbl_Louisiana_GulfCoast_Wetlands_InclForest_Oct2016_v2_VV.flt')
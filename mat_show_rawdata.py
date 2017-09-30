#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 20 13:22:21 2016

@author: asoulain
"""

# Import necessary libraries
import matplotlib as mpl
mpl.use('TkAgg')

import sys
from matplotlib import pyplot as plt
import pyfits
import os
import cPickle
import numpy as np
import matplotlib.gridspec as gridspec

# Save the file content in a temporary python binary file
_dBfile_matisse = '/tmp/data_matisse.dpy'
    
def readDb_matisse(name_file):
    global _dBfile_matisse
    if not os.path.exists(_dBfile_matisse):
        return None
    f = open(_dBfile_matisse)
    db = cPickle.load(f)
    f.close()
    if name_file in db.keys():
        return db[name_file]
    else:
        #print db.keys()
        return None

def writeDb_matisse(name_file, data):
    global _dBfile_matisse
    if not os.path.exists(_dBfile_matisse):
        #print ' > getData: creating', _dBfile_matisse
        db = {name_file:data}
    else:
        #print ' > getData: updating', _dBfile_matisse
        f = open(_dBfile_matisse)
        db = cPickle.load(f)
        f.close()
    db[name_file] = data
    f = open(_dBfile_matisse, 'wb')
    cPickle.dump(db, f, 2)
    f.close()
    return
    
def open_mat(name_file):
    tmp = readDb_matisse(name_file)
    if not tmp is None:
        #print ' > getData: %s already in database MATISSE'%name_file.split('/')[-1]
        return tmp
    hdu = pyfits.open(name_file)
    img_data = hdu['IMAGING_DATA']
    
    interf = img_data.data['DATA11'].mean(axis = 0)
    phot1 =  img_data.data['DATA13'].mean(axis = 0)
    phot2 =  img_data.data['DATA12'].mean(axis = 0)
    phot3 =  img_data.data['DATA9'].mean(axis = 0)
    phot4 =  img_data.data['DATA10'].mean(axis = 0)
        
    dic = {'INTERF':interf, 'PHOT1':phot1, 'PHOT2':phot2, 'PHOT3':phot3, 'PHOT4':phot4}
    writeDb_matisse(name_file, dic)
    return dic
    

#name_file = '/Users/asoulain/Documents/These/DRS_MATISSE/DATA/2016-09-28/Kappa-03/input/MATISSE_GEN_cal_kappa_N_HIGH_shut1_IMG_272_0001.fits'
#name_file = '/Users/asoulain/Documents/These/DRS_MATISSE/DATA/2016-09-21/Kappa-01/input/MATISSE_GEN_cal_kappa_N_LOW_shut1_DARK_265_0001.fits'
#name_file = '/Users/asoulain/Documents/These/DRS_MATISSE/DATA/2016-09-21/Kappa-01/test1.fits'

dis = True
try:
    name_file = sys.argv[1]
    dic = open_mat(name_file)
except IOError:
    dis= False
    print '\n #### File not found ####'
except KeyError:
    dis= False
    print "\n #### Extension 'IMAGING_DATA' not found ####"
except IndexError:
    print "\n #### Enter a MATISSE fits file name as argument ####"
    dis = False
    
if dis:
    if len(sys.argv) == 3:
        res = sys.argv[2]
    else:
        res = 'HIGH'
    
    if not res == '-l':
        plt.close('all')
        plt.figure(figsize = (9, 6))
        G = gridspec.GridSpec(2, 6)
        
        axes_1 = plt.subplot(G[:, 0])
        axes_1.imshow(dic['PHOT1'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_1.set_title('PHOT1')
        
        plt.xticks([]), plt.yticks([])
        
        axes_2 = plt.subplot(G[:,1])
        plt.xticks([]), plt.yticks([])
        axes_2.imshow(dic['PHOT2'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_2.set_title('PHOT2')
        
        axes_3 = plt.subplot(G[:,2:4])
        plt.xticks([]), plt.yticks([])
        axes_3.imshow(dic['INTERF'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_3.set_title('INTERF')
        
        axes_4 = plt.subplot(G[:,4])
        plt.xticks([]), plt.yticks([])
        axes_4.imshow(dic['PHOT3'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_4.set_title('PHOT3')
        
        axes_5 = plt.subplot(G[:,5])
        plt.xticks([]), plt.yticks([])
        axes_5.imshow(dic['PHOT4'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_5.set_title('PHOT4')
        plt.tight_layout()
        G.update(wspace=0.1)
        plt.show()
    else:
        plt.close('all')
        plt.figure(figsize = (12, 4))
        G = gridspec.GridSpec(2, 8)
        
        axes_1 = plt.subplot(G[:, 0])
        axes_1.imshow(dic['PHOT1'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_1.set_title('PHOT1')
        
        plt.xticks([]), plt.yticks([])
        
        axes_2 = plt.subplot(G[:,1])
        plt.xticks([]), plt.yticks([])
        axes_2.imshow(dic['PHOT2'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_2.set_title('PHOT2')
        
        axes_3 = plt.subplot(G[:,2:6])
        plt.xticks([]), plt.yticks([])
        axes_3.imshow(dic['INTERF'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_3.set_title('INTERF')
        
        axes_4 = plt.subplot(G[:,6])
        plt.xticks([]), plt.yticks([])
        axes_4.imshow(dic['PHOT3'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_4.set_title('PHOT3')
        
        axes_5 = plt.subplot(G[:,7])
        plt.xticks([]), plt.yticks([])
        axes_5.imshow(dic['PHOT4'], interpolation = 'nearest', cmap = 'afmhot', origin = 'down')
        axes_5.set_title('PHOT4')
        plt.tight_layout()
        G.update(wspace=0.1)
        plt.show()


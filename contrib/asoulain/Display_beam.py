# -*- coding: utf-8 -*-
"""
Created on Fri Sep 16 13:50:04 2016

@author: asoulain
"""

import numpy as np
from matplotlib import pyplot as plt
from astropy.io import fits
from glob import glob
from matplotlib.colors import PowerNorm

def Open_beam(band, beam, BCD):
    
    nbeam2 = str(beam)
    if band == 'L':
        datadir = '/Users/asoulain/Documents/These/DRS_MATISSE/DATA/2016-09-14/BeamPos-0'+nbeam2+'/input/'
    else:
        datadir = '/Users/asoulain/Documents/These/DRS_MATISSE/DATA/Beam_Nband/BeamPos-0'+nbeam2+'/'
        
    list_file = glob(datadir + '*.fits')
    list_lamp = [f for f in list_file if 'LAMP' in f]
    list_dark = [f for f in list_file if 'DARK' in f]
    
    if BCD == 'IN':
        n_file = -2
    else:
        n_file = -1
        
    hdu = fits.open(list_dark[n_file])
    hdu2 = fits.open(list_lamp[n_file])
    
    hdr = hdu2[0].header
    hdr_D = hdu[0].header
    
    if hdr['ESO INS BSL1 ST'] and hdr_D['ESO INS BSL1 ST']:# and hdr_D['ESO INS BSN1 ST']:
        nbeam = '1'
    elif hdr['ESO INS BSL2 ST'] and hdr_D['ESO INS BSL2 ST']:
        nbeam = '2'
    elif hdr['ESO INS BSL3 ST'] and hdr_D['ESO INS BSL3 ST']:
        nbeam = '3'
    elif hdr['ESO INS BSL4 ST'] and hdr_D['ESO INS BSL4 ST']:
        nbeam = '4'
    elif hdr['ESO INS BSN1 ST'] and hdr_D['ESO INS BSN1 ST']:
        nbeam = '1'
    elif hdr['ESO INS BSN2 ST'] and hdr_D['ESO INS BSN2 ST']:
        nbeam = '2'
    elif hdr['ESO INS BSN3 ST'] and hdr_D['ESO INS BSN3 ST']:
        nbeam = '3'
    elif hdr['ESO INS BSN4 ST'] and hdr_D['ESO INS BSN4 ST']:
        nbeam = '4'
    else:
        print 'prb'
        
    print '\n%s, %s, Beam %s, BCD : %s'%(hdr['OBJECT'], hdr['ESO DET NAME'], nbeam, hdr['ESO INS BCD1 NAME'])
    print '%s, %s, Beam %s, BCD : %s'%(hdr_D['OBJECT'], hdr_D['ESO DET NAME'], nbeam, hdr_D['ESO INS BCD1 NAME'])
    
    data = hdu[2].data
    data2 = hdu2[2].data
    
    dark = np.median(data['DATA5'], axis = 0)
    lamp = np.median(data2['DATA5'], axis = 0)
    
    beam = lamp - dark
    beam[beam < 0] = 0
    dic = {'lamp':lamp, 'beam':beam, 'dark':dark}
    return dic
    

def recad_beam(beam, sizex, sizey):
    c = np.where(beam == np.max(beam))
    beam_rec = beam[c[0][0]-sizex:c[0][0]+sizex, c[1][0]-sizey:c[1][0]+sizey]
    return beam_rec
    

if __name__ == '__main__':
    band = 'N'
    nbeam2 = str(1)
    BCD = 'OUT'

    plt.close('all')
    
    beam_1_in = recad_beam(Open_beam(band, 1, 'IN')['beam'], 10, 200)
    beam_1_out = recad_beam(Open_beam(band, 1, 'OUT')['beam'], 10, 200)
    beam_2_in = recad_beam(Open_beam(band, 2, 'IN')['beam'], 10, 200)
    beam_2_out = recad_beam(Open_beam(band, 2, 'OUT')['beam'], 10, 200)
    beam_3_in = recad_beam(Open_beam(band, 3, 'IN')['beam'], 10, 200)
    beam_3_out = recad_beam(Open_beam(band, 3, 'OUT')['beam'], 10, 200)
    beam_4_in = recad_beam(Open_beam(band, 4, 'IN')['beam'], 10, 200)
    beam_4_out = recad_beam(Open_beam(band, 4, 'OUT')['beam'], 10, 200)
    
    
    n_ligne = 8
    n_col   = 4
    fig = plt.figure(0,figsize = (18,18))
    ax1 = plt.subplot2grid((n_ligne,n_col), (0,0), colspan = 2)
    ax2 = plt.subplot2grid((n_ligne,n_col), (1,0), colspan = 2)
    ax3 = plt.subplot2grid((n_ligne,n_col), (2,0), colspan = 2)
    ax4 = plt.subplot2grid((n_ligne,n_col), (3,0), colspan = 2)
    ax5 = plt.subplot2grid((n_ligne,n_col), (4,0), colspan = 2)
    ax6 = plt.subplot2grid((n_ligne,n_col), (5,0), colspan = 2)
    ax7 = plt.subplot2grid((n_ligne,n_col), (6,0), colspan = 2)
    ax8 = plt.subplot2grid((n_ligne,n_col), (7,0), colspan = 2)
    
    ax9 = plt.subplot2grid((n_ligne,n_col), (0,2), colspan = 1, rowspan = 2)
    ax10 = plt.subplot2grid((n_ligne,n_col), (0,3), colspan = 1, rowspan = 2)
    
    ax11 = plt.subplot2grid((n_ligne,n_col), (2,2), colspan = 1, rowspan = 2)
    ax12 = plt.subplot2grid((n_ligne,n_col), (2,3), colspan = 1, rowspan = 2)
    
    ax13 = plt.subplot2grid((n_ligne,n_col), (4,2), colspan = 1, rowspan = 2)
    ax14 = plt.subplot2grid((n_ligne,n_col), (4,3), colspan = 1, rowspan = 2)
    
    ax15 = plt.subplot2grid((n_ligne,n_col), (6,2), colspan = 1, rowspan = 2)
    ax16 = plt.subplot2grid((n_ligne,n_col), (6,3), colspan = 1, rowspan = 2)
    
    ax1.set_title('Beam 1 - IN')
    ax2.set_title('Beam 1 - OUT')
    ax3.set_title('Beam 2 - IN')
    ax4.set_title('Beam 2 - OUT')
    ax5.set_title('Beam 3 - IN')
    ax6.set_title('Beam 3 - OUT')
    ax7.set_title('Beam 4 - IN')
    ax8.set_title('Beam 4 - OUT')
    
    ax9.set_title('No anamorphosed - IN')
    ax10.set_title('No anamorphosed - OUT')
    
    p = .5
    ax1.imshow(beam_1_in, interpolation = 'nearest', cmap = 'gist_stern', norm = PowerNorm(p), vmin = 0, vmax = 1000)
    ax2.imshow(beam_1_out, interpolation = 'nearest', cmap = 'gist_stern',norm = PowerNorm(p), vmin = 0, vmax = 1000)
    ax3.imshow(beam_2_in, interpolation = 'nearest', cmap = 'gist_stern',norm = PowerNorm(p), vmin = 0, vmax = 1000)
    #ax4.imshow(beam_2_out, interpolation = 'nearest', cmap = 'gist_stern',norm = PowerNorm(p))
    ax4.axis('off')
    ax4.set_xlim([0,400])
    ax4.set_ylim([98,118])
    ax5.imshow(beam_3_in, interpolation = 'nearest', cmap = 'gist_stern', norm = PowerNorm(p))
    ax6.imshow(beam_3_out, interpolation = 'nearest', cmap = 'gist_stern', norm = PowerNorm(p))
    ax7.imshow(beam_4_in, interpolation = 'nearest', cmap = 'gist_stern', norm = PowerNorm(p))
    ax8.imshow(beam_4_out, interpolation = 'nearest', cmap = 'gist_stern', norm = PowerNorm(p))
    
    cb9 = ax9.imshow(beam_1_in, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb9, ax = ax9)    
    cb10 = ax10.imshow(beam_1_out, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb10, ax = ax10) 
    cb11 = ax11.imshow(beam_2_in, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb11, ax = ax11) 
    ax12.axis('off')
    cb12 = ax12.imshow(beam_2_out, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    #plt.colorbar(cb12, ax = ax12)    
    ax12.set_xlim([0,400])
    ax12.set_ylim([98,118])

    cb13 = ax13.imshow(beam_3_in, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb13, ax = ax13)    
    cb14 = ax14.imshow(beam_3_out, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb14, ax = ax14)    
    cb15 = ax15.imshow(beam_4_in, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb15, ax = ax15)    
    cb16 = ax16.imshow(beam_4_out, interpolation = 'nearest', cmap = 'gist_stern', aspect = 24, norm = PowerNorm(p))
    plt.colorbar(cb16, ax = ax16)    

    plt.tight_layout()
#    plt.savefig('All_beam_'+band+'band_sqrt1.pdf')
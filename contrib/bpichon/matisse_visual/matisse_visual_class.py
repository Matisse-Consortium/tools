# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 11:41:26 2015
@author: kervella
@authot: bpichon 
"""
try:
   import pyfits
except:
   from astropy.io import fits as pyfits
import numpy as np

#==============================================================================
# Functions
#==============================================================================

def get_key_withdefault(dict, key, default):
    if key in dict.keys():
        return dict[key]
    if key in ["HIERARCH "+c for c in dict.keys()]:
        return dict[key]
    else:
        return default

class Oifits:
    def __init__(self, filename):
        # open fits file
        self.filename = filename
        matisse_oifits=pyfits.open(filename, mode='readonly')

        # Main FITS header
        self.header = matisse_oifits[0].header.copy()
        # self.exptime_sc = self.header['HIERARCH ESO DET2 SEQ1 DIT']
        self.single = False

#        if "COMBINED" in str(self.header['HIERARCH ESO FT POLA MODE']):
#            self.polarsplit = False
#        else:
#            self.polarsplit = True
        self.polarsplit = True
        
        # get the FDDL, OI_ARRAY, OI_WAVELENGTH and OI_VIS_MET extentions
        for hdu in matisse_oifits :
            # OI_TARGET
            if hdu.name == 'OI_TARGET': # test if OI_TARGET is present in header
                self.oi_target = matisse_oifits['OI_TARGET'].data
            if hdu.name == 'OI_ARRAY' :
                oi_array=hdu.data.copy()
#            if hdu.name == 'ARRAY_GEOMETRY' : # TBC, not present (2015-09)
#                oi_array=hdu.data.copy()

           
        # get the OI_VIS, OI_FLUX and T3 extensions in combined polarization mode
#essai BP
        for hdu in matisse_oifits:
#BP            header = hdu.header.copy()
            if (hdu.name == 'OI_WAVELENGTH'):
                oi_wave_sc=hdu.data.copy()
            if (hdu.name == 'OI_VIS2'):
                oi_vis2_sc=hdu.data.copy()
            if (hdu.name == 'OI_T3'):
                oi_t3_sc=hdu.data.copy()
            if (hdu.name == 'OI_VIS'):
                oi_vis_sc=hdu.data.copy()
            if (hdu.name == 'OI_FLUX'):
                oi_flux_sc=hdu.data.copy()
                self.oi_flux_sc_time = oi_flux_sc.field('TIME')[::4]/3600. # in decimal hours
                #headerflux = hdu.header.copy()
                #BP columnflux = hdu.columns.names
        
        matisse_oifits.close()
        del matisse_oifits
        # gc.collect()

        # Wavelength scales
        self.wave_sc=oi_wave_sc.field('EFF_WAVE')*1.e6 # wavelength scale in microns
        self.minwave_sc=np.min(self.wave_sc)
        self.maxwave_sc=np.max(self.wave_sc)
        self.nwave_sc = self.wave_sc.shape[0]


        # OI_ARRAY parameters
        self.array_tel_name = oi_array.field('TEL_NAME')
        self.array_sta_name = oi_array.field('STA_NAME')
        self.array_staxyz = oi_array.field('STAXYZ')
        

        self.stations = [get_key_withdefault(self.header,'HIERARCH ESO ISS CONF STATION4','S4'),\
                         get_key_withdefault(self.header,'HIERARCH ESO ISS CONF STATION3','S3'),\
                         get_key_withdefault(self.header,'HIERARCH ESO ISS CONF STATION2','S2'),\
                         get_key_withdefault(self.header,'HIERARCH ESO ISS CONF STATION1','S1')]
                         
        self.telnames = [get_key_withdefault(self.header,'HIERARCH ESO ISS CONF T4NAME','T4'),\
                         get_key_withdefault(self.header,'HIERARCH ESO ISS CONF T3NAME','T3'),\
                         get_key_withdefault(self.header,'HIERARCH ESO ISS CONF T2NAME','T2'),\
                         get_key_withdefault(self.header,'HIERARCH ESO ISS CONF T1NAME','T1')]

        # baseline names
        self.basenames = [self.stations[0]+self.stations[1],
                          self.stations[0]+self.stations[2],
                          self.stations[0]+self.stations[3],
                          self.stations[1]+self.stations[2],
                          self.stations[1]+self.stations[3],
                          self.stations[2]+self.stations[3]]

        # Time scale of the different records in the OIFITS file
        self.time = self.oi_flux_sc_time
        
        # Number of records
        self.nrecords = self.time.shape[0]

        self.oi_flux_sc = oi_flux_sc.field('FLUXDATA').reshape(self.nrecords,4,self.nwave_sc)
        self.oi_flux_err_sc = oi_flux_sc.field('FLUXERR').reshape(self.nrecords,4,self.nwave_sc)
       

        # OI_VIS and T3 table parameters
        self.oi_vis_sc_time = oi_vis_sc.field('TIME')[::6]/3600. # in decimal hours
        self.oi_vis_sc_mjd = oi_vis_sc.field('MJD')[::6]
        self.oi_vis_sc_target_id = oi_vis_sc.field('TARGET_ID')[::6]
        self.oi_vis_sc_int_time = oi_vis_sc.field('INT_TIME') # in seconds
        
        # Complex visibilities
        
        #BP self.oi_vis_sc_visdata = oi_vis_sc.field('VISDATA')
        #BP self.oi_vis_sc_viserr = oi_vis_sc.field('VISERR')
        self.oi_vis_sc_visamp = oi_vis_sc.field('VISAMP')
        self.oi_vis_sc_visamperr = oi_vis_sc.field('VISAMPERR')
        self.oi_vis_sc_visphi = oi_vis_sc.field('VISPHI')*np.pi/180. # conversion in radians (after Sept. 16)
        self.oi_vis_sc_visphierr = oi_vis_sc.field('VISPHIERR')*np.pi/180. # conversion in radians (after Sept. 16)
        self.oi_vis_sc_flag = oi_vis_sc.field('FLAG')

        # u,v coordinates different for SC and FT

        nbaseline= 6
        self.ucoord_sc = oi_vis_sc.field('UCOORD')
        self.vcoord_sc = oi_vis_sc.field('VCOORD')
        self.projbase_sc = np.zeros(nbaseline)
        self.azimbase_sc = np.zeros(nbaseline)
        self.spfreq_sc = np.zeros((nbaseline,self.nwave_sc))
        for baseline in range(0,nbaseline):
            self.projbase_sc[baseline] =  np.sqrt(np.square(self.ucoord_sc[baseline])+np.square(self.vcoord_sc[baseline]))
            self.azimbase_sc[baseline] =  np.arctan2(self.vcoord_sc[baseline], self.ucoord_sc[baseline]) * 180./np.pi
            self.spfreq_sc[baseline,:] = self.projbase_sc[baseline]/(self.wave_sc/1.E6)/(180.*3600./np.pi) # B/lambda cycles/arcsec

        # Station index common               
        
        self.sta_index = oi_vis_sc.field('STA_INDEX')
        
        # Squared visibilities
        
        self.oi_vis2_sc_vis2data = oi_vis2_sc.field('VIS2DATA')
        self.oi_vis2_sc_vis2err = oi_vis2_sc.field('VIS2ERR')
        self.oi_vis2_sc_ucoord = oi_vis2_sc.field('UCOORD')
        self.oi_vis2_sc_vcoord = oi_vis2_sc.field('VCOORD')
        self.oi_vis2_sc_sta_index = oi_vis2_sc.field('STA_INDEX')
        self.oi_vis2_sc_flag = oi_vis2_sc.field('FLAG')
        
        # Closure quantities

        self.t3_u1coord = oi_t3_sc.field('U1COORD')
        self.t3_v1coord = oi_t3_sc.field('V1COORD')
        self.t3_baseline1 = np.sqrt(np.square(self.t3_u1coord)+np.square(self.t3_v1coord))

        self.t3_u2coord = oi_t3_sc.field('U2COORD')
        self.t3_v2coord = oi_t3_sc.field('V2COORD')
        self.t3_baseline2 = np.sqrt(np.square(self.t3_u2coord)+np.square(self.t3_v2coord))

        self.t3_u3coord = self.t3_u1coord+self.t3_u2coord
        self.t3_v3coord = self.t3_v1coord+self.t3_v2coord
        self.t3_baseline3 = np.sqrt(np.square(self.t3_u3coord)+np.square(self.t3_v3coord))
        
        self.t3_baseavg = np.nanmean([self.t3_baseline1,self.t3_baseline2,self.t3_baseline3],axis=0)
        self.t3_basemax = np.nanmax([self.t3_baseline1,self.t3_baseline2,self.t3_baseline3],axis=0)
        self.t3_spfreq_sc = np.zeros((4,self.nwave_sc))
        for triplet in range(0,4):
            self.t3_spfreq_sc[triplet,:] = self.t3_basemax[triplet]/(self.wave_sc/1.E6)/(180.*3600./np.pi)
        
        self.t3phi_sc = oi_t3_sc.field('T3PHI')*np.pi/180. # conversion in radians (after Sept. 16)
        self.t3phierr_sc = oi_t3_sc.field('T3PHIERR')*np.pi/180. # conversion in radians (after Sept. 16)

        self.t3amp_sc = oi_t3_sc.field('T3AMP')
        self.t3amperr_sc = oi_t3_sc.field('T3AMPERR')

        # Station indices for the closure quantities
        self.t3_sc_staindex = oi_t3_sc.field('STA_INDEX')

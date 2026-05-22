# Tool for flux and correlated flux calibration of VLTI/MATISSE data
# Jozsef Varga, 2019-2022
# varga@strw.leidenuniv.nl
#
# Example usage:
# from fluxcal import fluxcal
# inputfile_sci = 'path/to/raw/science/fits/file.fits'
# inputfile_cal = 'path/to/raw/calibrator/fits/file.fits'
# outputfile = 'path/to/calibrated/outputfile.fits'
# cal_database_dir = 'path/to/calibrator/database/folder/'
# cal_database_paths = [cal_database_dir+'vBoekelDatabase.fits',cal_database_dir+'calib_spec_db_v10.fits',cal_database_dir+'calib_spec_db_v10_supplement.fits']
# output_fig_dir = 'path/to/figure/folder/'       
# fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_paths, mode='corrflux',output_fig_dir=output_fig_dir)
#
# Arguments:
# inputfile_sci: path to the raw science oifits file.
# inputfile_cal: path to the raw calibrator oifits file. 
# outputfile: path of the output calibrated file
# dir_caldatabases: Path to the directory containing all the calibrator spectra databases
# mode (optional): 
#   'flux': calibrates total flux (incoherently processed oifits file expected)
#           results written in the OI_FLUX table (FLUXDATA column)
#    'corrflux': calibrates correlated flux (coherently processed oifits file expected)
#                results written in the OI_VIS table (VISAMP column)
#    'both': calibrates both total and correlated fluxes

#
#TODO: calculate uncertainties
#       partly implemented
#       caveats: uncertainty in the calibrator spectrum is not taken into account
#                
#TODO: treat if the cal database cannot be opened
# treat if there is no matching source in the database: DONE
#TODO: FIX: some dec values are NULL in vBoekeldatabase 
#
########################################################################

import numpy as np
from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.io import fits
from shutil import copyfile
from matplotlib import pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.special import j0,j1,jv
from scipy.interpolate import interp1d
import math
import os
from astroquery.simbad import Simbad
from numpy.polynomial.polynomial import polyval
from astropy.convolution import Gaussian1DKernel,Box1DKernel,convolve
import scipy.stats
import glob
import shutil
import sys
from astroquery.vizier import Vizier
import subprocess

#Path to the skycal_cli executable
#skycalc_cli_cmd = '/home/amatter/.local/bin/skycalc_cli'
skycalc_cli_cmd = str(shutil.which('skycalc_cli'))

# match_radius [arcsec]
# ra, dec [degree]

def get_spectrum_cal_STARSFLUX(cal_name,out_lst,ra=np.nan,dec=np.nan,match_radius=15.0):
    try:
        print('Checking if '+cal_name+' is in MDFC and STARSFLUX')#
        c_cal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')       
        result = Vizier.query_object(cal_name, catalog='II/361')
        #result = Vizier.query_region(c_cal, catalog='II/361', radius=str(match_radius)+"s")
        cal_name_vizier=result[0]['Name'][0].replace(' ','')
        cal_name_vizier=cal_name_vizier.replace('*','')
        #extract calibrator model spectrum
        hdu=fits.open('https://home.strw.leidenuniv.nl/~gamez/catalog/fitsfilesFinal/Dec30'+cal_name_vizier+'.fits')
        #hdu=fits.open('https://home.strw.leidenuniv.nl/~gamez/catalog/fitsfiles/Aug20'+cal_name_vizier+'.fits')
        #hdu=fits.open('/home/amatter/Documents/VLTI/MATISSE/GPAO/LGS_commissioning/Run_Dec/Data_HLTau/Dec30HD26546.fits')
        wav_cal = hdu[3].data['Wavelength'] #um
        spectrum_cal = hdu[0].data
        diam_cal=hdu[0].header['ANGRMAS']*2
        diam_err_cal=diam_cal*0.1
        print('*******calibrator is in MDFC and STARSFLUX!!!****')
        newcat=True
        out_lst += [cal_name_vizier,diam_cal,diam_err_cal,
        wav_cal*1e-6,spectrum_cal,'STARSFLUX',3,
        result[0]['RAJ2000'][0],result[0]['DEJ2000'][0]]
        return True
    except:        
        print('Calibrator not in STARSFLUX')
        return False

def get_spectrum_caldb(cal_database_path,cal_name,out_lst,ra=np.nan,dec=np.nan,match_radius=15.0,band='L'):
    c_cal = SkyCoord(ra*u.deg, dec*u.deg, frame='icrs')
    caldb = fits.open(cal_database_path)
    caldb_file = os.path.basename(cal_database_path)
    if 'vBoekelDatabase' in caldb_file:
        if 'fitsold' in caldb_file:
            cal_name_lst = caldb[8].data['NAME']
            cal_ra_lst = caldb[8].data['RAEPP']
            cal_dec_lst = caldb[8].data['DECEPP']
        else: 
            cal_name_lst = caldb['SOURCES'].data['NAME']
            cal_ra_lst = caldb['SOURCES'].data['RAEPP']
            cal_dec_lst = caldb['SOURCES'].data['DECEPP']
    elif 'calib_spec_db' in caldb_file:
        cal_name_lst = caldb['SOURCES'].data['name']
        cal_ra_lst = caldb['SOURCES'].data['ra']
        cal_dec_lst = caldb['SOURCES'].data['dec']
    else:
        cal_name_lst = caldb['SOURCES'].data['NAME']
        cal_ra_lst = caldb['SOURCES'].data['RAEPP']
        cal_dec_lst = caldb['SOURCES'].data['DECEPP']

    c_lst = SkyCoord(cal_ra_lst * u.deg, cal_dec_lst * u.deg, frame='icrs')

    # search for the calibrator in the calibrator database
    sep = c_cal.separation(c_lst)
    min_sep_idx = np.nanargmin(sep)
    min_sep = sep[min_sep_idx]
    if (min_sep.value < match_radius/3600.0):
        #match
        print('Calibrator found in the database '+caldb_file+': '+cal_name_lst[min_sep_idx]+', separation: %.2f arcsec'%(3600.0*min_sep.value))
        #get calibrator diameter
        if 'vBoekelDatabase' in caldb_file:
            offset = 9
            if 'fitsold' in caldb_file:
                diam_cal =  1000.0*caldb[-2].data['DIAMETER'][min_sep_idx] #mas
                diam_err_cal = 1000.0*caldb[-2].data['DIAMETER_ERR'][min_sep_idx] #mas
            else:
                diam_cal =  1000.0*caldb['DIAMETERS'].data['DIAMETER'][min_sep_idx] #mas
                diam_err_cal = 1000.0*caldb['DIAMETERS'].data['DIAMETER_ERR'][min_sep_idx] #mas
        elif 'calib_spec_db' in caldb_file:
            offset = 2
            if 'L' in band:
                diam_cal = caldb['SOURCES'].data['UDDL_est'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_est'][min_sep_idx] #mas
            if 'N' in band:
                diam_cal = caldb['SOURCES'].data['UDDN_est'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_est'][min_sep_idx] #mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['diam_midi'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_midi'][min_sep_idx] #mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['diam_cohen'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_cohen'][min_sep_idx] #mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['UDD_meas'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['e_diam_meas'][min_sep_idx] #mas
            if math.isnan(diam_cal):
                diam_cal = caldb['SOURCES'].data['diam_gaia'][min_sep_idx] #mas
                diam_err_cal = diam_cal*0.1 #mas 
        #extract calibrator model spectrum
        wav_cal = caldb[min_sep_idx+offset].data['WAVELENGTH'] #m
        spectrum_cal = caldb[min_sep_idx+offset].data['FLUX']
        out_lst += [caldb[min_sep_idx+offset].header['NAME'],diam_cal,diam_err_cal,
            wav_cal,spectrum_cal,caldb_file,3600.0*min_sep.value,
            cal_ra_lst[min_sep_idx],cal_dec_lst[min_sep_idx]]
        caldb.close()
        return True
    else:
        print('Calibrator not found in '+caldb_file)
        print('Closest match: '+cal_name_lst[min_sep_idx]+', separation: %.2f arcsec'%(3600.0*min_sep.value))
        caldb.close()
        return False

# def plot_spec(cal_name,cal_database_path,fig_dir ='.',ra=np.nan,dec=np.nan,match_radius=15.0,wl_lim=(np.nan,np.nan),xlog=False,ylog=False):
#     if math.isnan(ra):
#         Simbad.reset_votable_fields()
#         Simbad.remove_votable_fields('coordinates')
#         Simbad.add_votable_fields('ra(d)', 'dec(d)')
#         res =  Simbad.query_object(cal_name)
#         ra = res['RA_d'][0]
#         dec = res['DEC_d'][0] 
#     out_lst =[]
#     match = get_spectrum_cal_STARSFLUX(cal_name,out_lst,ra=ra,dec=dec,match_radius=match_radius)
#     if match == False:
#         match = get_spectrum_caldb(cal_database_path,cal_name,out_lst,ra=ra,dec=dec,match_radius=match_radius,band='L')
#     if match:
#         cal_name_db,diam_cal,diam_err_cal,wav_cal,spectrum_cal,caldb_file,min_sep_arcsec,ra_cal_db,dec_cal_db = out_lst
#         fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
#         plt.plot(wav_cal*1e6, spectrum_cal, '-k',lw=2.0)  
#         plt.ylabel('Flux (Jy)')
#         plt.xlabel(r'$\lambda$ ($\mu$m)')
#         plt.title(cal_name_db+r', diam = $%.2f \pm %.2f$ mas'%(diam_cal,diam_err_cal))
#         if not math.isnan(wl_lim[0]): 
#             ax1.set_xlim(wl_lim)
#         else:
#             wl_lim=(np.min(wav_cal*1e6),np.max(wav_cal*1e6))
#         wl_idx = np.logical_and(wav_cal*1e6 > wl_lim[0], wav_cal*1e6 < wl_lim[1])
#         if xlog:
#             ax1.set_xscale('log')
#         if ylog:
#             ax1.set_yscale('log')
#         else:
#             ax1.set_ylim((0.0,1.1*np.max(spectrum_cal[wl_idx])))
#         fig.savefig(fig_dir+'/' +cal_name_db.replace(' ','_')+'_%.1d-%.1d_um.png'%(wl_lim[0],wl_lim[1]), dpi=200)
#         # plt.show()

def plot_skycalc_output(fpath_in,figpath_out,airmass,pwv):
    hdul = fits.open(fpath_in)
    wl_um = hdul[1].data['lam']*1e-3
    trans = hdul[1].data['trans']
    fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
    plt.plot(wl_um , trans,lw=1.5)  
    plt.ylabel('Transmission')
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.title('Airmass = %.3f, pwv = %.2f mm'%(airmass,pwv))
    ax1.set_ylim([0.0,1.0])
    fig.savefig(figpath_out, dpi=200)
    hdul.close()

def plot_skycalc_MATISSE(spec_skycalc_final,figpath_out,airmass,pwv,spec_MATISSE,wave_MATISSE):
    fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
    plt.plot(wave_MATISSE, spec_skycalc_final,lw=1.5,label='Skycalc convolved and resampled') 
    plt.plot(wave_MATISSE, spec_MATISSE/np.nanmax(spec_MATISSE),lw=1.5,label='Normalized MATISSE raw spectrum') 
    plt.ylabel('Transmission')
    plt.xlabel(r'$\lambda$ ($\mu$m)')
    plt.title('Airmass = %.3f, pwv = %.2f mm'%(airmass,pwv))
    ax1.set_ylim([0.0,1.0])
    plt.legend()
    #plt.show()
    fig.savefig(figpath_out, dpi=200)

    
#mode: 'flux','corrflux','both'
def fluxcal(inputfile_sci, inputfile_cal, outputfile, dir_caldatabases, 
    mode='flux',output_fig_dir='',match_radius=15.0,do_airmass_correction=False,calc_spectrum_offset=False):

    # create the output oifits file
    copyfile(inputfile_sci, outputfile)
    outhdul = fits.open(outputfile, mode='update')

    # read in the input oifits files
    try:
        inhdul_sci = fits.open(inputfile_sci)
    except FileNotFoundError as e:
        print('Target reduced data file not found.')
        outhdul.close()
        os.remove(outputfile)
        return 3    
    try:
        inhdul_cal = fits.open(inputfile_cal)
    except FileNotFoundError as e:
        print('Calibrator reduced data file not found.')
        inhdul_sci.close()
        outhdul.close()
        os.remove(outputfile)
        return 2
    
    try:
        cal_database_paths = glob.glob(dir_caldatabases+'/*.fits')
        ind_vBoekel=cal_database_paths.index(dir_caldatabases+'/vBoekelDatabase.fits')
        ind_v10=cal_database_paths.index(dir_caldatabases+'/calib_spec_db_v10.fits')
        ind_v10_supplement=cal_database_paths.index(dir_caldatabases+'/calib_spec_db_v10_supplement.fits')
        cal_database_paths = [cal_database_paths[ind_vBoekel],cal_database_paths[ind_v10],cal_database_paths[ind_v10_supplement]] 
    except FileNotFoundError as e:
        print('Calibrator spectra databases not found.')
        return 2

    airmass_sci = (inhdul_sci[0].header['HIERARCH ESO ISS AIRM START']+inhdul_sci[0].header['HIERARCH ESO ISS AIRM END'])/2.0
    pwv_sci = (inhdul_sci[0].header['HIERARCH ESO ISS AMBI IWV30D START']+inhdul_sci[0].header['HIERARCH ESO ISS AMBI IWV30D END'])/2.0
    sci_name = inhdul_sci['OI_TARGET'].data['TARGET'][0]
    header = inhdul_sci[0].header
    detector=header['HIERARCH ESO DET CHIP NAME']

    # extract calibrator information
    cal_name = inhdul_cal['OI_TARGET'].data['TARGET'][0]
    ra_cal = inhdul_cal['OI_TARGET'].data['RAEP0'][0]
    dec_cal = inhdul_cal['OI_TARGET'].data['DECEP0'][0]
    airmass_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AIRM START']+inhdul_cal[0].header['HIERARCH ESO ISS AIRM END'])/2.0
    pwv_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AMBI IWV30D START']+inhdul_cal[0].header['HIERARCH ESO ISS AMBI IWV30D END'])/2.0
    seeing_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AMBI FWHM START']+inhdul_cal[0].header['HIERARCH ESO ISS AMBI FWHM END'])/2.0
    tau0_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AMBI TAU0 START']+inhdul_cal[0].header['HIERARCH ESO ISS AMBI TAU0 END'])/2.0
    tpl_start_cal = inhdul_cal[0].header['HIERARCH ESO TPL START']
    band = inhdul_cal[0].header['HIERARCH ESO DET CHIP TYPE'] #'IR-LM' or 'IR-N'
    print('SCI: %s, airmass = %.2f'%(sci_name,airmass_sci))
    print('CAL: %s, RA = %.6f Dec = %.6f, airmass = %.2f'%(cal_name,ra_cal,dec_cal,airmass_cal))
    #print('CAL TPL START: '+tpl_start_cal)
    
    wav_cal = inhdul_cal['OI_WAVELENGTH'].data['EFF_WAVE']  # m
    wav_sci = inhdul_sci['OI_WAVELENGTH'].data['EFF_WAVE']  # m
    
    # Computation of synthetic skycalc spectra
    if do_airmass_correction:
        airmass_sci_input = airmass_sci
        airmass_cal_input = airmass_cal
        pwv_sci_input = pwv_sci
        pwv_cal_input = pwv_cal 
              
        print('Do airmass correction.')
        outputdir = os.path.dirname(outputfile)+'/skycalc/'
        if not os.path.exists(outputdir):
            os.makedirs(outputdir)
        tag_sci = os.path.splitext(os.path.basename(inputfile_sci))[0]
                
        wmin = np.min(wav_sci)*1e6-0.02 #[um]
        wmax = np.max(wav_sci)*1e6+0.02
        #margin = 0.1*(wmax-wmin)
        #wmin = wmin-margin
        #wmax = wmax+margin
        dlambda = get_dlambda(inhdul_sci)
        
        fname_in_sci = outputdir+'/skycalc_input_sci_'+tag_sci+'.txt'
        fname_skycalc_out_sci = outputdir+'/skycalc_output_sci_'+tag_sci+'.fits'
        create_skycalc_inputfile(fname_in_sci,airmass_sci_input,pwv_sci_input,wmin,wmax)
        print('Start SkyCalc (SCI).')
        os.system(skycalc_cli_cmd + ' -i ' + fname_in_sci + ' -o ' + fname_skycalc_out_sci)
        
        #figpath_out = outputdir+'/skycalc_output_sci_'+tag_sci+'.png'
        #plot_skycalc_output(fname_skycalc_out_sci,figpath_out,airmass_sci_input,pwv_sci_input)
        
        wmin = np.min(wav_sci)*1e6-0.02 #[um]
        wmax = np.max(wav_sci)*1e6+0.02
        #margin = 0.1*(wmax-wmin)
        #wmin = wmin-margin
        #wmax = wmax+margin
        dlambda = get_dlambda(inhdul_cal)

        tag_cal = os.path.splitext(os.path.basename(inputfile_cal))[0]
        fname_in_cal = outputdir+'/skycalc_input_cal_'+tag_cal+'.txt'
        fname_skycalc_out_cal = outputdir+'/skycalc_output_cal_'+tag_cal+'.fits'
        create_skycalc_inputfile(fname_in_cal,airmass_cal_input,pwv_cal_input,wmin,wmax)
        print('Start SkyCalc (CAL).')
        os.system(skycalc_cli_cmd + ' -i ' + fname_in_cal + ' -o ' + fname_skycalc_out_cal)
        #figpath_out = outputdir+'/skycalc_output_cal_'+tag_cal+'.png'
        #plot_skycalc_output(fname_skycalc_out_cal,figpath_out,airmass_cal_input,pwv_cal_input)
    
        
    # open the calibrator database which includes the spectra of calibrators
    out_lst = []
    match = get_spectrum_cal_STARSFLUX(cal_name,out_lst,ra=ra_cal,dec=dec_cal,match_radius=match_radius)

    if match == False:
        i = 0
        while (match == False and i < len(cal_database_paths)):
            match = get_spectrum_caldb(cal_database_paths[i],cal_name,out_lst,ra=ra_cal,dec=dec_cal,match_radius=match_radius,band=band)
            i = i+1
  
    if match == True:
        cal_name_db,diam_cal,diam_err_cal,wav_cal_model,spectrum_cal,caldb_file,min_sep_arcsec,ra_cal_db,dec_cal_db = out_lst
        if math.isnan(diam_cal):
            print('Calibrator diameter not found.')
            if mode != 'flux':
                inhdul_cal.close()
                inhdul_sci.close()
                outhdul.flush()  # changes are written back to fits
                outhdul.close()
                os.remove(outputfile)
                return 4
        print('Diameter = %.2f +/- %.2f mas (%s)'%(diam_cal,diam_err_cal,cal_name_db))
        
        if 'calib_spec' in caldb_file:
            wav_cal_model = np.flip(wav_cal_model)
            spectrum_cal = np.flip(spectrum_cal) #Jy
    else:
        print('Calibrator not found in any of the databases')
        inhdul_cal.close()
        inhdul_sci.close()
        outhdul.flush()  # changes are written back to fits
        outhdul.close()
        os.remove(outputfile)
        return 1

    
    if np.min(wav_cal) < 0.0:
        print('ERROR (fluxcal): Wavelength grid in oifits file invalid.')
        inhdul_cal.close()
        inhdul_sci.close()
        outhdul.flush()  # changes are written back to fits
        outhdul.close()
        os.remove(outputfile)
        return 5

    if 'L' in band:
        wav_cal = np.flip(wav_cal)
        wav_sci = np.flip(wav_sci)
        
        
    # smooth and resample calibrator spectrum to match the spectral resolution and wavelength grid of the data
    spectral_binning = get_spectral_binning(inhdul_cal)
    idx = np.logical_and(wav_cal_model > 0.99*np.nanmin(wav_cal), wav_cal_model < 1.01*np.nanmax(wav_cal))
    spectrum_cal = spectrum_cal[idx]
    wav_cal_model = wav_cal_model[idx]
    dwave_sci = np.mean(np.gradient(wav_sci))
    dwav_cal_model = np.mean(np.gradient(wav_cal_model))
    if os.path.exists(output_fig_dir):
        figpath_out = output_fig_dir+'/' + 'calibrator_'+band+'_'+cal_name.replace(' ','_')+'_spectrum_resampling.png'
        make_plot = True
    else:
        figpath_out = ''
        make_plot = False
    if dwav_cal_model < dwave_sci:  #MATISSE spectrum is sparser than synthetic spectrum of the calibrator
        spectrum_cal_resampled  = transform_spectrum_to_real_spectral_resolution(wav_cal_model*1e6,spectrum_cal,wav_sci*1e6,detector,make_plot=make_plot,figpath_out=figpath_out,ylabel='Flux (Jy)',yscale='final')
    else:
        f = interp1d(wav_cal_model, spectrum_cal,kind='cubic')
        spectrum_cal_resampled = f(wav_sci)

    # plot calibrator spectrum in the data wavelength range
    # if os.path.exists(output_fig_dir):
    #     fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
    #     if 2.0*len(wav_cal) < np.nansum(idx):
    #         plt.plot(wav_cal_model*1e6, spectrum_cal, '-',color='grey',label='Original synthetic spectrum',lw=1.0,alpha=0.66)  
    #     else:
    #         plt.plot(wav_cal_model*1e6, spectrum_cal, '-o',color='grey',label='Original synthetic spectrum',lw=1.0,alpha=0.66)
    #     plt.plot(wav_cal*1e6,spectrum_cal_resampled,'-r',label='Convolved and resampled spectrum')
        
        
    if do_airmass_correction:
        #print('Calculate airmass correction factor')
        #read in the saved atmospheric transmission curves
        hdulc = fits.open(fname_skycalc_out_cal)
        wl_um_cal = hdulc[1].data['lam']*1e-3
        trans_cal = hdulc[1].data['trans']
        hduls = fits.open(fname_skycalc_out_sci)
        wl_um_sci = hduls[1].data['lam']*1e-3
        trans_sci = hduls[1].data['trans']
        
        #transform the transmission to the wavelength grid of the data (incl. spectral resolution & spectral binning)
        figpath_out = outputdir+'/skycalc_spec_transform_sci_'+tag_sci+'.png'        
        trans_sci_final  = transform_spectrum_to_real_spectral_resolution(wl_um_sci,trans_sci,wav_sci*1e6,detector,
        make_plot=True,figpath_out=figpath_out,ylabel='Transmission')

        figpath_out = outputdir+'/skycalc_spec_transform_cal_'+tag_cal+'.png'
        trans_cal_final  = transform_spectrum_to_real_spectral_resolution(wl_um_cal,trans_cal,wav_cal*1e6,detector,
        make_plot=True,figpath_out=figpath_out,ylabel='Transmission')

        #txtpath_out = outputdir+'/skycalc__'+tag_cal+'.png'
        #wertret
        #write_spectrum_txt(trans_cal_final,wav_cal*1e6,txtpath_out)
        
        #calculate the airmass_correction_factor
        airmass_correction_factor = trans_cal_final/trans_sci_final
        
        #plot the airmass_correction_factor
        tag_sci_final = os.path.splitext(os.path.basename(outputfile))[0]
        figa, ((axc)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        plt.plot(wav_sci*1e6 ,airmass_correction_factor) 
        wl_min = np.min(wl_um_cal)
        wl_max = np.max(wl_um_cal)
        wl_range = wl_max - wl_min
        LM_gap_idx = np.logical_and(wav_sci*1e6 > 4.1 ,wav_sci*1e6 < 4.6)
        idx = np.logical_and(wav_sci*1e6 > wl_min+wl_range*0.1,wav_sci*1e6 < wl_max-wl_range*0.1)
        range_idx = np.logical_and(idx,~LM_gap_idx)
        pmin = np.nanmin(airmass_correction_factor[range_idx])  #np.nanpercentile(y_new, 10.0)
        pmax = np.nanmax(airmass_correction_factor[range_idx])  #np.nanpercentile(y_new, 95.0)
        axc.set_ylim([pmin,pmax])
        axc.set_ylabel('Airmass correction factor')
        axc.set_xlabel(r'$\lambda$ ($\mu$m)')
        #ax1.set_ylim([0.0,1.0])
        figpath_out = outputdir+'/skycalc_airmass_correction_factor_'+tag_sci_final+'.png'
        figa.savefig(figpath_out, dpi=200)
        hduls.close()        
        hdulc.close()
    else:
        airmass_correction_factor = wav_sci*0.0+1.0
    
    # calibrate total spectrum
    if mode == 'flux' or mode == 'both':
        print('Calibrate total spectrum.')
        #check if we have an 'OI_FLUX' table
        n_exp_sci = len(inhdul_sci['OI_VIS2'].data['VIS2DATA'])/6
        n_exp_cal = len(inhdul_cal['OI_VIS2'].data['VIS2DATA'])/6
        rp_list_sci=[]
        rp_list_cal=[]
        if do_airmass_correction:
            outputdir_comp = os.path.dirname(outputfile)+'/skycalc/comparison_MATISSE_skycalc/'
            if not os.path.exists(outputdir_comp):
                os.makedirs(outputdir_comp)
            
        try:
            outhdul['OI_FLUX'].header['TUNIT5'] = 'Jy'
            outhdul['OI_FLUX'].header['TUNIT6'] = 'Jy'
            nexp=len(inhdul_sci['OI_FLUX'].data['FLUXDATA'])
            for j in range(nexp):
                ind_tel = j % 4 + 1
                ind_exp = j // 4 + 1
                flux_raw_sci = inhdul_sci['OI_FLUX'].data['FLUXDATA'][j]
                flux_raw_cal = inhdul_cal['OI_FLUX'].data['FLUXDATA'][j]
                fluxerr_raw_sci = inhdul_sci['OI_FLUX'].data['FLUXERR'][j]
                fluxerr_raw_cal = inhdul_cal['OI_FLUX'].data['FLUXERR'][j]
                if 'L' in band:
                    flux_raw_sci = np.flip(flux_raw_sci)
                    flux_raw_cal = np.flip(flux_raw_cal)
                    fluxerr_raw_sci = np.flip(fluxerr_raw_sci)
                    fluxerr_raw_cal = np.flip(fluxerr_raw_cal)

                if do_airmass_correction:
                    figpath_out = outputdir_comp+'/skycalc_comparison_cal_'+tag_cal+str(ind_exp)+'_tel'+str(ind_tel)+'.png'
                    plot_skycalc_MATISSE(trans_cal_final,figpath_out,airmass_cal_input,pwv_cal_input,flux_raw_cal,wav_cal*1e6)
                    # if calc_spectrum_offset:
                    #     #calculate correlation between the raw spectrum, and the atmospheric transmission spectrum
                    #     #shift_max = int(0.025*len(trans_cal_final))
                    #     shift_max = int(7.0*spectral_binning)
                    #     rp_list_cal.append(calc_corr_offset(trans_cal_final,flux_raw_cal,shift_max))
                    #     rp_list_sci.append(calc_corr_offset(trans_sci_final,flux_raw_sci,shift_max))
                    
                flux_calibrated_sci = flux_raw_sci/flux_raw_cal*spectrum_cal_resampled*n_exp_cal/n_exp_sci*airmass_correction_factor
                if 'L' in band:
                    flux_calibrated_sci = np.flip(flux_calibrated_sci)
                fluxerr_calibrated_sci = np.abs(flux_raw_sci/flux_raw_cal*n_exp_cal/n_exp_sci)* \
                    np.sqrt((fluxerr_raw_sci*n_exp_cal/n_exp_sci/flux_raw_sci)**2 + (fluxerr_raw_cal/flux_raw_cal)**2)* \
                    spectrum_cal_resampled
                if 'L' in band:
                    fluxerr_calibrated_sci = np.flip(fluxerr_calibrated_sci)
                outhdul['OI_FLUX'].data['FLUXDATA'][j] = flux_calibrated_sci
                outhdul['OI_FLUX'].data['FLUXERR'][j] = fluxerr_calibrated_sci

            
            # if do_airmass_correction:
            #     if calc_spectrum_offset:
            #         #plot the correlation  
            #         plot_corr_offset(rp_list_sci,np.arange(-shift_max,+shift_max),outputdir+'/skycalc_correlation_flux_'+tag_sci+'.png')
            #         plot_corr_offset(rp_list_cal,np.arange(-shift_max,+shift_max),outputdir+'/skycalc_correlation_flux_'+tag_cal+'.png')

        except KeyError as e:
            print('No OI_FLUX table found.')
    

    # calibrate correlated spectrum
    if mode == 'corrflux' or mode == 'both':
        print('Calibrate correlated spectra.')
        rp_list_sci=[]
        rp_list_cal=[]
        outhdul['OI_VIS'].header['TUNIT5'] = 'Jy'
        outhdul['OI_VIS'].header['TUNIT6'] = 'Jy'        
        for j in range(len(inhdul_sci['OI_VIS'].data['VISAMP'])):
            sta_index_sci = inhdul_sci['OI_VIS'].data['STA_INDEX'][j]
            corrflux_raw_sci = inhdul_sci['OI_VIS'].data['VISAMP'][j]
            corrfluxerr_raw_sci = inhdul_sci['OI_VIS'].data['VISAMPERR'][j]
            if 'L' in band:
                corrflux_raw_sci = np.flip(corrflux_raw_sci)
                corrfluxerr_raw_sci = np.flip(corrfluxerr_raw_sci)

            # find calibrator data with matching station configuration
            sta_indices_cal = inhdul_cal['OI_VIS'].data['STA_INDEX']
            for i in range(len(sta_indices_cal)):
                if ((sta_index_sci[0] == sta_indices_cal[i][0]) and (sta_index_sci[1] == sta_indices_cal[i][1])) \
                or ((sta_index_sci[0] == sta_indices_cal[i][1]) and (sta_index_sci[1] == sta_indices_cal[i][0])):
                    idx_cal = i
                    break

            corrflux_raw_cal = inhdul_cal['OI_VIS'].data['VISAMP'][idx_cal]
            corrfluxerr_raw_cal = inhdul_cal['OI_VIS'].data['VISAMPERR'][idx_cal]
            if 'L' in band:
                corrflux_raw_cal = np.flip(corrflux_raw_cal)
                corrfluxerr_raw_cal = np.flip(corrfluxerr_raw_cal)
            uu = inhdul_cal['OI_VIS'].data['UCOORD'][idx_cal]
            vv = inhdul_cal['OI_VIS'].data['VCOORD'][idx_cal]
            B_p = np.sqrt(uu**2 + vv**2)

            diam_cal_rad = diam_cal/1000.0/3600.0*math.pi/180.0
            spatial_frequency = B_p/wav_cal
            # visibilities of the calibrator (uniform disk model)
            vis_cal = 2*j1(math.pi*diam_cal_rad*spatial_frequency) / (math.pi*diam_cal_rad*spatial_frequency)
            # plt.figure()
            # plt.plot(wav_cal, vis_cal, '-b')
            # plt.show()
            if do_airmass_correction:
                if calc_spectrum_offset:
                    #calculate correlation between the raw spectrum, and the atmospheric transmission spectrum
                    #shift_max = int(0.025*len(trans_cal_final))
                    shift_max = int(7.0*spectral_binning)
                    if len(trans_cal_final) == len(corrflux_raw_cal):
                        rp_list_cal.append(calc_corr_offset(trans_cal_final,corrflux_raw_cal,shift_max))
                        rp_list_sci.append(calc_corr_offset(trans_sci_final,corrflux_raw_sci,shift_max))

            corrflux_calibrated_sci = corrflux_raw_sci/corrflux_raw_cal*vis_cal*spectrum_cal_resampled*airmass_correction_factor
            if 'L' in band:
                corrflux_calibrated_sci = np.flip(corrflux_calibrated_sci)
            corrfluxerr_calibrated_sci = np.abs(corrflux_raw_sci/corrflux_raw_cal)* \
                np.sqrt((corrfluxerr_raw_sci/corrflux_raw_sci)**2 + (corrfluxerr_raw_cal/corrflux_raw_cal)**2)* \
                vis_cal*spectrum_cal_resampled #uncertainty in calibrator diameter is not taken into account
            if 'L' in band:
                corrfluxerr_calibrated_sci = np.flip(corrfluxerr_calibrated_sci)
            #if os.path.exists(output_fig_dir):
            #    plt.plot(wav_cal*1e6,vis_cal*spectrum_cal_resampled,'--',label='Resampl. corr., B_p = %.2f m'%B_p)
            # plt.figure()
            # plt.plot(wav_cal, vis_cal, '-b')
            # plt.show()
            outhdul['OI_VIS'].data['VISAMP'][j] = corrflux_calibrated_sci
            outhdul['OI_VIS'].data['VISAMPERR'][j] = corrfluxerr_calibrated_sci

        if do_airmass_correction:
            #plot the correlation  
            if calc_spectrum_offset:
                plot_corr_offset(rp_list_sci,np.arange(-shift_max,+shift_max),outputdir+'/skycalc_correlation_corrflux_'+tag_sci+'.png')
                plot_corr_offset(rp_list_cal,np.arange(-shift_max,+shift_max),outputdir+'/skycalc_correlation_corrflux_'+tag_cal+'.png')
            
    outhdul[0].header['HIERARCH ESO PRO CAL NAME'] = cal_name
    outhdul[0].header['HIERARCH ESO PRO CAL RA'] = (ra_cal, '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL DEC'] = (dec_cal, '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL AIRM'] = airmass_cal
    outhdul[0].header['HIERARCH ESO PRO CAL IWV'] = pwv_cal
    outhdul[0].header['HIERARCH ESO PRO CAL FWHM'] = (seeing_cal, '[arcsec]')
    outhdul[0].header['HIERARCH ESO PRO CAL TAU0'] = (tau0_cal, 'Coherence time [s]')
    outhdul[0].header['HIERARCH ESO PRO CAL TPL START'] = tpl_start_cal
    outhdul[0].header['HIERARCH ESO PRO CAL DB NAME'] = (cal_name_db, 'Name of calibrator in cal database.')
    outhdul[0].header['HIERARCH ESO PRO CAL DB DBNAME'] = (caldb_file, 'Name of cal database')
    outhdul[0].header['HIERARCH ESO PRO CAL DB RA'] = (ra_cal_db, '[deg]') 
    outhdul[0].header['HIERARCH ESO PRO CAL DB DEC'] = (dec_cal_db , '[deg]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB DIAM'] = (diam_cal, 'Calibrator diameter [mas]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB ERRDIAM'] = (diam_err_cal, 'Error in calibrator diameter [mas]')
    outhdul[0].header['HIERARCH ESO PRO CAL DB SEP'] = (min_sep_arcsec, 'Separation (input coord - calDB coord) [arcsec]')
    outhdul[0].header['HIERARCH ESO PRO CATG'] = 'TARGET_FLUXCAL_INT'
    
    # if os.path.exists(output_fig_dir):    
    #     plt.ylabel('Flux (Jy)')
    #     plt.xlabel(r'$\lambda$ ($\mu$m)')
    #     plt.title(cal_name)
    #     ax1.set_xlim([np.min(wav_cal)*1e6,np.max(wav_cal)*1e6])
    #     ax1.set_ylim([0.7*np.min(spectrum_cal_resampled),1.1*np.max(spectrum_cal_resampled)])
    #     ax1.legend(loc='best', fontsize=8, fancybox=True, framealpha=0.5)
    #     fig.savefig(output_fig_dir+'/' + 'calibrator_'+band+'_'+cal_name.replace(' ','_')+'_spectrum.png', dpi=200)
    #     # plt.show()

    inhdul_cal.close()
    inhdul_sci.close()
    outhdul.flush()  # changes are written back to fits
    outhdul.close()
    return 0

def update_corrflux_from_vis2(inputfile, outputfile, flux_Jy):
    # create the output oifits file
    print(os.path.basename(inputfile))
    copyfile(inputfile, outputfile)
    inhdul = fits.open(inputfile)
    outhdul = fits.open(outputfile, mode='update')
    for j in range(len(outhdul['OI_VIS'].data['VISAMP'])):
        vis = np.sqrt(inhdul['OI_VIS2'].data['VIS2DATA'][j])
        viserr = np.abs(0.5*inhdul['OI_VIS2'].data['VIS2ERR'][j]/np.sqrt(inhdul['OI_VIS2'].data['VIS2DATA'][j]))
        corrflux = vis*flux_Jy
        corrfluxerr = viserr*flux_Jy
        outhdul['OI_VIS'].data['VISAMP'][j] = corrflux
        outhdul['OI_VIS'].data['VISAMPERR'][j] = corrfluxerr
    outhdul.flush()  # changes are written back to fits
    outhdul.close()    

def update_vis2_from_corrflux(inputfile, outputfile, total_flux):
    # create the output oifits file
    print(os.path.basename(inputfile))
    copyfile(inputfile, outputfile)
    inhdul = fits.open(inputfile)
    outhdul = fits.open(outputfile, mode='update')
    for j in range(len(outhdul['OI_VIS2'].data['VIS2DATA'])):
        corrflux = inhdul['OI_VIS'].data['VISAMP'][j] 
        corrfluxerr = inhdul['OI_VIS'].data['VISAMPERR'][j] 
        vis = corrflux/total_flux
        viserr = corrfluxerr/total_flux
        outhdul['OI_VIS2'].data['VIS2DATA'][j] = vis**2
        outhdul['OI_VIS2'].data['VIS2ERR'][j] = 2.0*vis*viserr
    outhdul.flush()  # changes are written back to fits
    outhdul.close()  

# #wmin, wmax [nm]
# #wgrid_mode: 'fixed_wavelength_step' or 'fixed_spectral_resolution'
# #lsf_type: 'none' or 'Gaussian' or 'Boxcar'
# def create_skycalc_inputfile(fname,airmass,pwv,wmin,wmax,
#         wgrid_mode='fixed_wavelength_step',wdelta=0.0001,wres=300000.0,lsf_type='none'):
#     pwv_allowed_values = [0.05,0.1,0.25,0.5,1.0,1.5,2.5,3.5,5.0,7.5,10.0,20.0,30.0]
#     #print(pwv)
#     ofile = open(fname,'w') 
#     mystr='airmass         :  '+'%f\n'%airmass
#     mystr+= 'pwv_mode        :  pwv \n'
#     mystr+='season          :  0 \n'
#     mystr+='time            :  0 \n'
#     idx = find_nearest_idx(pwv_allowed_values, pwv)
#     #print(idx)
#     mystr+= 'pwv             :  ' +'%f\n'%pwv_allowed_values[idx]
#     mystr+= 'msolflux        :  130.0\n'
#     mystr+= 'incl_moon       :  Y\n'
#     mystr+= 'moon_sun_sep    :  90.0\n'
#     mystr+= 'moon_target_sep :  45.0\n'
#     mystr+= 'moon_alt        :  45.0\n'
#     mystr+= 'moon_earth_dist :  1.0\n'
#     mystr+= 'incl_starlight  :  Y\n'
#     mystr+= 'incl_zodiacal   :  Y\n'
#     mystr+= 'ecl_lon         :  135.0\n'
#     mystr+= 'ecl_lat         :  90.0\n'
#     mystr+= 'incl_loweratm   :  Y\n'
#     mystr+= 'incl_upperatm   :  Y\n'
#     mystr+= 'incl_airglow    :  Y\n'
#     mystr+= 'incl_therm      :  N\n'
#     mystr+= 'therm_t1        :  0.0\n'
#     mystr+= 'therm_e1        :  0.0\n'
#     mystr+= 'therm_t2        :  0.0\n'
#     mystr+= 'therm_e2        :  0.0\n'
#     mystr+= 'therm_t3        :  0.0\n'
#     mystr+= 'therm_e3        :  0.0\n'
#     mystr+= 'vacair          :  vac\n'
#     mystr+= 'wmin            :  '+'%f\n'%wmin
#     mystr+= 'wmax            :  '+'%f\n'%wmax
#     mystr+= 'wgrid_mode      :  '+'%s\n'%wgrid_mode
#     mystr+= 'wdelta          :  '+'%f\n'%wdelta
#     mystr+= 'wres            :  '+'%f\n'%wres
#     mystr+= 'lsf_type        :  '+'%s\n'%lsf_type
# #    mystr+= 'lsf_gauss_fwhm  :  '+'%f\n'%lsf_gauss_fwhm
#     mystr+= 'lsf_boxcar_fwhm :  5.0\n'
#     mystr+= 'observatory     :  paranal'
#     ofile.write(mystr)
#     ofile.close()
    
def create_skycalc_inputfile(
    fname, airmass, pwv, wmin, wmax,
    wgrid_mode='fixed_wavelength_step',
    wdelta=0.0001,          # µm (high-res)
    wres=300000.0,
    lsf_type='none'
):
    pwv_allowed = np.array([0.05,0.1,0.25,0.5,1.0,1.5,2.5,3.5,5.0,7.5,10.0,20.0,30.0])
    pwv = pwv_allowed[np.argmin(np.abs(pwv_allowed - pwv))]

    with open(fname, 'w') as f:
        f.write(f"""
airmass         : {airmass}
pwv_mode        : pwv
pwv             : {pwv}
season          : 0
time            : 0
msolflux        : 130.0
incl_loweratm   : Y
incl_upperatm   : Y
incl_airglow    : Y
incl_therm      : N
vacair          : vac
wmin            : {wmin*1e3}
wmax            : {wmax*1e3}
wgrid_mode      : {wgrid_mode}
wdelta          : {wdelta*1e3}
wres            : {wres}
lsf_type        : {lsf_type}
observatory     : paranal
""")
    

#return dlambda [nm]
def get_dlambda(inhdul):
    header = inhdul[0].header
    dl = np.nan
    if 'AQUARIUS' in header['HIERARCH ESO DET CHIP NAME']:
        #band = 'N'
        dispname = header['HIERARCH ESO INS DIN NAME']
        if 'LOW' in dispname:
            dl = 30.0  # in nm (factor 10 oversampling compared to true MATISSE spectral resolution)
        if 'HIGH' in dispname:
            dl = 3.0  # in nm (factor 10 oversampling compared to true MATISSE spectral resolution)
    if 'HAWAII' in header['HIERARCH ESO DET CHIP NAME']:
        #band = 'LM'
        dispname = header['HIERARCH ESO INS DIL NAME']
        if 'LOW' in dispname:
            dl = 8.0 # in nm (factor 10 oversampling compared to true MATISSE spectral resolution)
        if 'MED' in dispname:
            dl = 0.6 # in nm (factor 10 oversampling compared to true MATISSE spectral resolution)
        if 'HIGH' in dispname:
            if '+' in dispname:
                dl = 0.1 # in nm (factor 10 oversampling compared to true MATISSE spectral resolution)
            else: 
                dl = 0.3 # in nm (factor 10 oversampling compared to true MATISSE spectral resolution)
    return dl

def get_spectral_binning(inhdul):
    header = inhdul[0].header
    spectral_binning = np.nan
    for i in range(1,20):
        hdrkey = 'HIERARCH ESO PRO REC1 PARAM'+'%d'%i+' NAME' 
        if hdrkey in header:
            if 'spectralBinning' in header[hdrkey]:
                hdrkey2 = 'HIERARCH ESO PRO REC1 PARAM'+'%d'%i+' VALUE' 
                spectral_binning = float(header[hdrkey2])
                #print(spectral_binning)
    return spectral_binning

def get_dl_coeffs(inhdul):
    header = inhdul[0].header
    dl = np.nan
    if 'AQUARIUS' in header['HIERARCH ESO DET CHIP NAME']:
        #band = 'N'
        dispname = header['HIERARCH ESO INS DIN NAME']
        if 'LOW' in dispname:
            dl_coeffs = [ -43.925423443201, 0.163383818251532, -0.000113452755620802]
        if 'HIGH' in dispname:
            dl_coeffs = [7.43089207111461, 0.00585108729841899, 2.56662505487663E-07]
    if 'HAWAII' in header['HIERARCH ESO DET CHIP NAME']:
        #band = 'LM'
        dispname = header['HIERARCH ESO INS DIL NAME']
        if 'LOW' in dispname:
            dl_coeffs = [15.5921852886677, -0.227526028058492, 0.00191185004041472, -7.48990363774737E-06, 1.05762333960091E-08]
        if 'MED' in dispname:
            dl_coeffs = [ 5.355347142078, -1.429557244816e-03,  1.328209861214e-9]
        if 'HIGH' in dispname:
            if '+' in dispname:
                filtername=header['HIERARCH ESO INS FIL NAME']
                if 'M' in filtername:
                    dl_coeffs = [5.15608712708281, -0.000297668023525613, 4.00969202840045E-09, -1.72964617418913E-14, 9.31730249019226E-20]
                else:
                    dl_coeffs = [4.13933317850807, -0.000251977632, 1.13811555E-08]
            else: 
                dl_coeffs = [4.28216626046355, -0.000791610248938923, 3.02152795807675E-08]
    return dl_coeffs


def transform_spectrum_to_real_spectral_resolution(
        wl_orig,
        spec_orig,
        wl_final,
        detector,
        nsigma=5,
        make_plot=False,
        figpath_out='spectrum_convolved_resampled.png',
        ylabel='',
        yscale='none'):

    spec_final = np.zeros_like(wl_final)

    dlam = np.abs(np.gradient(wl_final))
    for i, lam0 in enumerate(wl_final):
        #lam0 = lam_out[i] + 0.5*dlam[i]

        # local MATISSE resolution element
        if 'HAWAII' in detector:
            fwhm = 4.92 * dlam[i]
        else:
            fwhm = 7.87 * dlam[i]
        sigma = fwhm / 2.355

        mask = np.abs(wl_orig - lam0) < nsigma * sigma

        l = wl_orig[mask]
        t = spec_orig[mask]

        if len(l) < 5:
            spec_final[i] = np.nan
            continue

        w = np.exp(-0.5 * ((l - lam0)/sigma)**2)
        w /= np.sum(w)

        spec_final[i] = np.sum(t * w)
        
    if make_plot:
        fig, ax = plt.subplots()
        ax.plot(wl_orig, spec_orig, label='Original synthetic spectrum')
        ax.plot(wl_final, spec_final, label='Convolved and resampled spectrum')
        ax.set_ylabel(ylabel)
        ax.set_xlabel(r'$\lambda$ ($\mu$m)')
        ax.legend()
        #plt.show()
        fig.savefig(figpath_out, dpi=200)
        plt.close(fig)
        
    return spec_final


    
# def transform_spectrum_to_real_spectral_resolution(wl_orig,spec_orig,dl_coeffs,kernel_width_px,wl_final,spectral_binning,
#     make_plot=False,figpath_out='spectrum_plots.png',ylabel=''):
#     #make an uneven wavelength grid
#     min_wl = np.min(wl_orig)
#     max_wl = np.max(wl_orig)
#     wl_new=[]
#     wl = min_wl
#     wl_new.append(wl)
#     while wl < max_wl:
#         wl = wl + polyval(wl,dl_coeffs)/kernel_width_px
#         wl_new.append(wl)
#     wl_new = np.array(wl_new)
#     #print(wl_orig,wl_new)
#     #interpolate the original spectrum to the new grid
#     f_spec_new = interp1d(wl_orig, spec_orig,kind='cubic',fill_value='extrapolate')
#     spec_new = f_spec_new(wl_new)
#     #convolve with Gaussian kernel 
#     kernel = Gaussian1DKernel(stddev=kernel_width_px/(2.0*np.sqrt(2.0*np.log(2.0))))
#     #spec_convolved = convolve(spec_new,kernel)
#     spec_new[0] = np.nanmedian(spec_new[0:int(kernel.dimension/2.0)])
#     spec_new[-1] = np.nanmedian(spec_new[-1:-int(kernel.dimension/2.0)])
#     spec_convolved = convolve(spec_new,kernel,boundary='extend')
#     #interpolate the convolved spectrum to the input wavelength grid
#     f_spec_new = interp1d(wl_new, spec_convolved,kind='cubic',fill_value='extrapolate')
#     spec_interp = f_spec_new(wl_final)
#     #apply spectral binning: convolve with a top-hat kernel of size spectral_binning
#     if spectral_binning > 1:
#         kernel = Box1DKernel(spectral_binning)
#         spec_final = convolve(spec_interp,kernel,boundary='extend')
#     else:
#         spec_final = spec_interp
        
#     if make_plot:
#         fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(20, 10))
#         plt.plot(wl_orig, spec_orig,label='1) original sp.')  
#         plt.plot(wl_new, spec_new,label='2) regridded sp.') 
#         plt.plot(wl_new, spec_convolved,label=r'3) convolved sp.'+r' ($R \approx %.0f$)'%np.nanmean(wl_final/polyval(wl_final,dl_coeffs))) 
#         plt.plot(wl_final, spec_interp,label='4) convolved interp sp.')
#         plt.plot(wl_final, spec_final,label='5) binned sp.'+' (%d px)'%spectral_binning)  
#         plt.ylabel(ylabel)
#         plt.xlabelr('$\lambda$ ($\mu$m)')
#         #ax1.set_ylim([0.0,1.0])
#         plt.legend()
#         fig.savefig(figpath_out, dpi=200)
      
#     return spec_final

def calc_corr_offset(spectrum1,spectrum2,shift_max):
    spectrum2_nonan=spectrum2[np.isfinite(spectrum2)]
    spectrum1_nonan=spectrum1[np.isfinite(spectrum2)]
    Ntr  = len(spectrum1_nonan)
    Ntr2 = len(spectrum2_nonan)
    rp = []
    for k in range(-shift_max,+shift_max):
        if k < 0:
            rp.append(scipy.stats.pearsonr(spectrum1_nonan[0:(Ntr+k)],spectrum2_nonan[-k:])[0]) #Pearson's r
        else:
            rp.append(scipy.stats.pearsonr(spectrum1_nonan[k:],spectrum2_nonan[0:(Ntr-k)])[0]) #Pearson's r
    return rp


def plot_corr_offset(rp_list,x,figpath_out):
    if len(rp_list) > 0:     
        fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        for rp in rp_list:
            plt.step(x,np.array(rp),where='mid') 
        plt.ylabel('Correlation')
        plt.xlabel('Wavelength offset (px)')
        tlocs, tlabels = plt.xticks()
        tspacing = tlocs[1]-tlocs[0]
        if tspacing < 20.0:
            minor_locator = MultipleLocator(1)
        elif tspacing < 40.0:
            minor_locator = MultipleLocator(2)
        elif tspacing < 80.0:
            minor_locator = MultipleLocator(5)
        else:
            minor_locator = MultipleLocator(10)
        ax1.xaxis.set_minor_locator(minor_locator)
        plt.grid(axis='x',color='0.95',which='minor')
        plt.grid(axis='x',color='0.8',which='major')
        fig.savefig(figpath_out, dpi=200)
                    
def find_nearest_idx(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

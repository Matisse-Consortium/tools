# Tool for flux and correlated flux calibration of VLTI/MATISSE data
# Jozsef Varga, 2019
# varga@strw.leidenuniv.nl
#
# Example usage:
# from fluxcal import fluxcal
# inputfile_sci = 'path/to/raw/science/fits/file.fits'
# inputfile_cal = 'path/to/raw/calibrator/fits/file.fits'
# outputfile = 'path/to/calibrated/outputfile.fits'
# cal_database_dir = 'path/to/calibrator/database/folder/'
# cal_database_paths = [cal_database_dir+'vBoekelDatabase.fits',cal_database_dir+'calib_spec_db.fits',cal_database_dir+'calib_spec_db_new.fits']
# output_fig_dir = 'path/to/figure/folder/'       
# fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_paths, mode='corrflux',output_fig_dir=output_fig_dir)
#
# Arguments:
# inputfile_sci: path to the raw science oifits file.
# inputfile_cal: path to the raw calibrator oifits file. 
# outputfile: path of the output calibrated file
# cal_database_paths: list of paths to the calibrator databases, e.g., [caldb1_path,caldb2_path] 
# mode (optional): 
#   'flux': calibrates total flux (incoherently processed oifits file expected)
#           results written in the OI_FLUX table (FLUXDATA column)
#    'corrflux': calibrates correlated flux (coherently processed oifits file expected)
#                results written in the OI_VIS table (VISAMP column)
#    'both': calibrates both total and correlated fluxes
# output_fig_dir (optional): if it is a valid path, the script will make a plot of the calibrator model spectrum there, 
#                            deafult: '' (= no figure made)
# 
#
# load vanBoekeldatabase: DONE
# save calibrator spectrum in a plot: DONE
# airmass correction: DONE
#TODO: calculate uncertainties
#       partly implemented
#       caveats: uncertainty in the calibrator spectrum is not taken into account
#                uncertainty in calibrator diameter is not taken into account
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
from scipy.special import j0,j1
from scipy.interpolate import interp1d
import math
import os

#mode: 'flux','corrflux','both'
def fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_paths, mode='flux',output_fig_dir=''):

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

    airmass_sci = (inhdul_sci[0].header['HIERARCH ESO ISS AIRM START']+inhdul_sci[0].header['HIERARCH ESO ISS AIRM END'])/2.0
    sci_name = inhdul_sci['OI_TARGET'].data['TARGET'][0]

    # extract calibrator information
    cal_name = inhdul_cal['OI_TARGET'].data['TARGET'][0]
    ra_cal = inhdul_cal['OI_TARGET'].data['RAEP0'][0]
    dec_cal = inhdul_cal['OI_TARGET'].data['DECEP0'][0]
    airmass_cal = (inhdul_cal[0].header['HIERARCH ESO ISS AIRM START']+inhdul_cal[0].header['HIERARCH ESO ISS AIRM END'])/2.0
    band = inhdul_cal[0].header['HIERARCH ESO DET CHIP TYPE'] #'IR-LM' or 'IR-N'
    c_cal = SkyCoord(ra_cal*u.deg, dec_cal*u.deg, frame='icrs')
    print('SCI: %s, airmass = %.2f CAL: %s, RA = %.6f Dec = %.6f, airmass = %.2f'%(sci_name,airmass_sci,cal_name,ra_cal,dec_cal,airmass_cal))

    # open the calibrator database which includes the spectra of calibrators
    match = False
    for cal_database_path in cal_database_paths:
        # print(cal_database_path)
        caldb = fits.open(cal_database_path)
        caldb_file = os.path.basename(cal_database_path)
        if 'fitsold' in caldb_file:
            cal_name_lst = caldb[8].data['NAME']
            cal_ra_lst = caldb[8].data['RAEPP']
            cal_dec_lst = caldb[8].data['DECEPP']
        else:            
            cal_name_lst = caldb['SOURCES'].data['NAME']
            cal_ra_lst = caldb['SOURCES'].data['RAEPP']
            cal_dec_lst = caldb['SOURCES'].data['DECEPP']
        c_lst = SkyCoord(cal_ra_lst * u.deg, cal_dec_lst * u.deg, frame='icrs')
        # print(c_lst)
        # search for the calibrator in the calibrator database
        sep = c_cal.separation(c_lst)
        min_sep_idx = np.nanargmin(sep)
        min_sep = sep[min_sep_idx]
        if (min_sep < 5.0*u.deg/3600.0):
            match = True
            break
        else:
            print('Calibrator not found in '+caldb_file)
            print('Closest match: '+cal_name_lst[min_sep_idx]+', separation: %.2f arcsec'%(3600.0*min_sep.value))
            
    if match == True:
        #match
        print('Calibrator found in the database '+caldb_file+': '+cal_name_lst[min_sep_idx]+', separation: %.2f arcsec'%(3600.0*min_sep.value))
        #get calibrator diameter
        if 'vBoekel' in caldb_file:
            offset = 9
            if 'fitsold' in caldb_file:
                diam_cal =  1000.0*caldb[-2].data['DIAMETER'][min_sep_idx] #mas
                diam_err_cal = 1000.0*caldb[-2].data['DIAMETER_ERR'][min_sep_idx] #mas
            else:
                diam_cal =  1000.0*caldb['DIAMETERS'].data['DIAMETER'][min_sep_idx] #mas
                diam_err_cal = 1000.0*caldb['DIAMETERS'].data['DIAMETER_ERR'][min_sep_idx] #mas
        if 'calib_spec' in caldb_file:
            offset = 2
            if 'L' in band:
                diam_cal = caldb['SOURCES'].data['UDDL'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['E_LDD'][min_sep_idx] #mas
            if 'N' in band:
                diam_cal = caldb['SOURCES'].data['UDDN'][min_sep_idx] #mas
                diam_err_cal = caldb['SOURCES'].data['E_LDD'][min_sep_idx] #mas
        print('Diameter = %.1f +/- %.1f mas (%s)'%(diam_cal,diam_err_cal,caldb[min_sep_idx+offset].header['NAME']))
        wav_cal = caldb[min_sep_idx+offset].data['WAVELENGTH'] #m
        if 'calib_spec' in caldb_file:
            wav_cal = np.flip(wav_cal)
        wl = np.concatenate((np.array([wav_cal[0] - (wav_cal[1] - wav_cal[0])]), wav_cal))
        wh = np.concatenate((wav_cal, np.array([wav_cal[-1] + (wav_cal[-1] - wav_cal[-2])])))
        wm = (wh + wl) / 2
        wav_bin_lower_cal = wm[:-1]
        wav_bin_upper_cal = wm[1:]
        d_wav_cal = wav_bin_upper_cal - wav_bin_lower_cal

        #extract calibrator model spectrum
        spectrum_cal = caldb[min_sep_idx+offset].data['FLUX']
        if 'calib_spec' in caldb_file:
            spectrum_cal = np.flip(spectrum_cal) #Jy
    else:
        print('Calibrator not found in any of the databases')
        caldb.close()
        inhdul_cal.close()
        inhdul_sci.close()
        outhdul.flush()  # changes are written back to fits
        outhdul.close()
        os.remove(outputfile)
        return 1

    wav = inhdul_cal['OI_WAVELENGTH'].data['EFF_WAVE']  # m
    if 'L' in band:
        wav = np.flip(wav)

    idx=np.logical_and(wav_cal > np.nanmin(wav) , wav_cal < np.nanmax(wav))
    # print(len(wav_cal),np.nansum(idx),len(wav))
    if 2.0*len(wav) < np.nansum(idx):
        #if the sampling of the MATISSE spectrum is much sparser than the sampling of the calibrator spectrum:
        wl = np.concatenate((np.array([wav[0] - (wav[1] - wav[0])]), wav))
        wh = np.concatenate((wav,np.array([wav[-1] + (wav[-1] - wav[-2])])))
        wm = (wh+wl)/2
        wav_bin_lower = wm[:-1]
        wav_bin_upper = wm[1:]
        d_wav = wav_bin_upper - wav_bin_lower
        # print(wav*1e6)
        # print(wav_cal*1e6)

        #resample model spectrum to the wavelengths of the observation
        spectrum_cal_resampled = wav*0.0
        flux_calibrated_sci = wav*0.0
        corrflux_calibrated_sci = wav*0.0
        vis_cal = wav*0.0
        # print(wav)
        # idx=np.logical_and(wav_cal<13.0e-6,wav_cal>8.0e-6)
        # print(wav_cal[idx]*1e6)
        # print(wav_bin_upper_cal,wav_bin_lower_cal)
        for i in range(len(wav)):
            # print(wav_bin_lower[i],wav_bin_upper[i])
            wu = (wav_bin_upper_cal - wav_bin_lower[i])
            wl = (wav_bin_upper[i] - wav_bin_lower_cal)
            wi = np.where(np.logical_and(wu > 0.0,wl > 0.0))
            # print(wav_cal[wi])
            # print(wi)
            wi = wi[0]
            # print(wi)

            #sum up the spectral values within the wavelength bin with weighting
            sum = 0.0
            sum = sum + (wav_bin_upper_cal[wi[0]] - wav_bin_lower[i])*spectrum_cal[wi[0]]
            sum = sum + (wav_bin_upper[i] - wav_bin_lower_cal[wi[-1]]) * spectrum_cal[wi[-1]]
            for j in range(1,len(wi)-1):
                sum = sum + spectrum_cal[wi[j]]*d_wav_cal[wi[j]]
            spectrum_cal_resampled[i] = sum/d_wav[i]
            # print(i,len(wi),sum,d_wav[i],spectrum_cal_resampled[i])
    else:
        #if the sampling of the MATISSE spectrum comparable to or denser than the sampling of the calibrator spectrum
        #do an interpolation
        f = interp1d(wav_cal, spectrum_cal,kind='cubic')
        spectrum_cal_resampled = f(wav)

    # print(wav_cal[idx])
    # print(wav)
    # print(spectrum_cal[idx])
    # print(spectrum_cal_resampled)
    # plot calibrator spectrum in the data wavelength range
    if os.path.exists(output_fig_dir):
        fig, ((ax1)) = plt.subplots(1, 1, sharey=False, sharex=False, figsize=(5, 5))
        if 2.0*len(wav) < np.nansum(idx):
            plt.plot(wav_cal*1e6, spectrum_cal, '-',color='grey',label='Original sp.',lw=1.0,alpha=0.66)  
        else:
            plt.plot(wav_cal*1e6, spectrum_cal, '-o',color='grey',label='Original sp.',lw=1.0,alpha=0.66)
        plt.plot(wav*1e6,spectrum_cal_resampled,'-r',label='Resampled total sp.')

    # calibrate total spectrum
    if mode == 'flux' or mode == 'both':
        #check if we have an 'OI_FLUX' table
        try:
            for j in range(len(inhdul_sci['OI_FLUX'].data['FLUXDATA'])):
                flux_raw_sci = inhdul_sci['OI_FLUX'].data['FLUXDATA'][j]*np.exp(airmass_sci)
                flux_raw_cal = inhdul_cal['OI_FLUX'].data['FLUXDATA'][j]*np.exp(airmass_cal)
                fluxerr_raw_sci = inhdul_sci['OI_FLUX'].data['FLUXERR'][j]
                fluxerr_raw_cal = inhdul_cal['OI_FLUX'].data['FLUXERR'][j]
                if 'L' in band:
                    flux_raw_sci = np.flip(flux_raw_sci)
                    flux_raw_cal = np.flip(flux_raw_cal)
                    fluxerr_raw_sci = np.flip(fluxerr_raw_sci)
                    fluxerr_raw_cal = np.flip(fluxerr_raw_cal)

                flux_calibrated_sci = flux_raw_sci/flux_raw_cal*spectrum_cal_resampled
                if 'L' in band:
                    flux_calibrated_sci = np.flip(flux_calibrated_sci)
                fluxerr_calibrated_sci = np.abs(flux_raw_sci/flux_raw_cal)* \
                    np.sqrt((fluxerr_raw_sci/flux_raw_sci)**2 + (fluxerr_raw_cal/flux_raw_cal)**2)* \
                    spectrum_cal_resampled
                if 'L' in band:
                    fluxerr_calibrated_sci = np.flip(fluxerr_calibrated_sci)
                outhdul['OI_FLUX'].data['FLUXDATA'][j] = flux_calibrated_sci
                outhdul['OI_FLUX'].data['FLUXERR'][j] = fluxerr_calibrated_sci
            outhdul['OI_FLUX'].header['TUNIT5'] = 'Jy'
            outhdul['OI_FLUX'].header['TUNIT6'] = 'Jy'
        except KeyError as e:
            print('No OI_FLUX table found.')

    # calibrate correlated spectrum
    if mode == 'corrflux' or mode == 'both':
        for j in range(len(inhdul_sci['OI_VIS'].data['VISAMP'])):
            sta_index_sci = inhdul_sci['OI_VIS'].data['STA_INDEX'][j]
            corrflux_raw_sci = inhdul_sci['OI_VIS'].data['VISAMP'][j]*np.exp(airmass_sci)
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

            corrflux_raw_cal = inhdul_cal['OI_VIS'].data['VISAMP'][idx_cal]*np.exp(airmass_cal)
            corrfluxerr_raw_cal = inhdul_cal['OI_VIS'].data['VISAMPERR'][idx_cal]
            if 'L' in band:
                corrflux_raw_cal = np.flip(corrflux_raw_cal)
                corrfluxerr_raw_cal = np.flip(corrfluxerr_raw_cal)
            uu = inhdul_cal['OI_VIS'].data['UCOORD'][idx_cal]
            vv = inhdul_cal['OI_VIS'].data['VCOORD'][idx_cal]
            B_p = np.sqrt(uu**2 + vv**2)

            diam_cal_rad = diam_cal/1000.0/3600.0*math.pi/180.0
            spatial_frequency = B_p/wav
            # visibilities of the calibrator (uniform disk model)
            vis_cal = 2*j1(math.pi*diam_cal_rad*spatial_frequency) / (math.pi*diam_cal_rad*spatial_frequency)
            # plt.figure()
            # plt.plot(wav, vis_cal, '-b')
            # plt.show()
            corrflux_calibrated_sci = corrflux_raw_sci/corrflux_raw_cal*vis_cal*spectrum_cal_resampled
            if 'L' in band:
                corrflux_calibrated_sci = np.flip(corrflux_calibrated_sci)
            corrfluxerr_calibrated_sci = np.abs(corrflux_raw_sci/corrflux_raw_cal)* \
                np.sqrt((corrfluxerr_raw_sci/corrflux_raw_sci)**2 + (corrfluxerr_raw_cal/corrflux_raw_cal)**2)* \
                vis_cal*spectrum_cal_resampled #uncertainty in calibrator diameter is not taken into account
            if 'L' in band:
                corrfluxerr_calibrated_sci = np.flip(corrfluxerr_calibrated_sci)
            if os.path.exists(output_fig_dir):
                plt.plot(wav*1e6,vis_cal*spectrum_cal_resampled,'--',label='Resampl. corr., B_p = %.2f m'%B_p)
            # plt.figure()
            # plt.plot(wav, vis_cal, '-b')
            # plt.show()
            outhdul['OI_VIS'].data['VISAMP'][j] = corrflux_calibrated_sci
            outhdul['OI_VIS'].data['VISAMPERR'][j] = corrfluxerr_calibrated_sci
        outhdul['OI_VIS'].header['TUNIT5'] = 'Jy'
        outhdul['OI_VIS'].header['TUNIT6'] = 'Jy'

    if os.path.exists(output_fig_dir):    
        plt.ylabel('Flux (Jy)')
        plt.xlabel('$\lambda$ ($\mu$m)')
        plt.title(cal_name)
        ax1.set_xlim([np.min(wav)*1e6,np.max(wav)*1e6])
        ax1.set_ylim([0.7*np.min(spectrum_cal_resampled),1.1*np.max(spectrum_cal_resampled)])
        ax1.legend(loc='best', fontsize=8, fancybox=True, framealpha=0.5)
        fig.savefig(output_fig_dir+'/' + 'calibrator_'+cal_name.replace(' ','_')+'_spectrum.png', dpi=200)
        # plt.show()

    caldb.close()
    inhdul_cal.close()
    inhdul_sci.close()
    outhdul.flush()  # changes are written back to fits
    outhdul.close()
    return 0

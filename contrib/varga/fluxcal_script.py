from fluxcal import fluxcal
from shutil import copyfile
import numpy as np
import glob 
import os
from show_allred import show_allred

basedir = '/allegro1/matisse/varga/'
# basedir = '/data2/archive/data/MATISSE'
reduced_dir = basedir + '/matisse_red5/'
# DATADIR = '/data2/archive/data/MATISSE/rawdata/'
# RESDIR_MAIN = '/data2/archive/data/MATISSE/matisse_red/'
targetdir = basedir +'/targets/'
targetlist_path = targetdir+'mat_target_list.txt'

# cal_database_dir = basedir+ '/CalibMap/'
# cal_database_paths = [cal_database_dir+'vBoekelDatabase.fits',cal_database_dir+'calib_spec_db.fits',cal_database_dir+'calib_spec_db_new.fits']
# inputfile_sci='/allegro1/matisse/varga/matisse_red5/coherent/2019-05-14/2019-05-15T08_36_04/Iter1/mat_raw_estimates.2019-05-15T08_36_04.HAWAII-2RG.rb/TARGET_RAW_INT_0001.fits'
# inputfile_cal='/allegro1/matisse/varga/matisse_red5/coherent/2019-05-14/2019-05-15T09_59_06/Iter1/mat_raw_estimates.2019-05-15T09_59_06.HAWAII-2RG.rb/CALIB_RAW_INT_0001.fits'
# outputfile = inputfile_sci.replace('RAW','FLUXCAL')
# output_fig_dir = '/allegro1/matisse/varga/matisse_red5/coherent/2019-05-14/2019-05-15T08_36_04/Iter1/mat_raw_estimates.2019-05-15T08_36_04.HAWAII-2RG.rb/'
# fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_paths, mode='corrflux',output_fig_dir=output_fig_dir)
# EXTERMINATE

do_calibration = True
do_make_plots = True
do_copy = True
do_copy_plots = False

#read target list
data = np.genfromtxt(targetlist_path,comments="#",delimiter=';',encoding=None,
    dtype=[('f0', '<U32'), ('f1', '<U128'), ('f2', '<f8'), ('f3', '<f8'), ('f4', '<f8'), ('f5', '<f8'),
    ('f6', '<U32'), ('f7', '<f8'), ('f8', '<f8'), ('f9', '<f8'), ('f10', '<f8'), ('f11', '<f8'),
    ('f12', '<f8'), ('f13', '<f8'), ('f14', '<f8'), ('f15', '<U128'), ('f16', '<U32')])
target_names = data['f0']
target_tags = data['f1']
#print(target_names.size)
try:
    len(target_names)
except TypeError as e:
    target_names = [target_names]
    target_tags = [target_tags]

# flst_sci = ['2019-03-23T084119_HD163296_IR-LM_IN_IN.fits',
# '2019-03-23T084119_HD163296_IR-LM_IN_OUT.fits',
# '2019-03-23T084119_HD163296_IR-LM_OUT_IN.fits',
# '2019-03-23T084119_HD163296_IR-LM_OUT_OUT.fits']

# flst_cal = ['2019-03-23T085341_delSgr_IR-LM_IN_IN.fits',
# '2019-03-23T085341_delSgr_IR-LM_IN_OUT.fits',
# '2019-03-23T085341_delSgr_IR-LM_OUT_IN.fits',
# '2019-03-23T085341_delSgr_IR-LM_OUT_OUT.fits']

#obs_lst: list of dics: [obsname_sci,night_sci,tplstart_sci,obsname_Lcal,night_Lcal,tplstart_Lcal,obsname_Ncal,night_Ncal,tplstart_Ncal]
obs_lst=[

    {'night_sci':'2019-09-18','tpl_start_sci':'2019-09-19T00:39:18',
    'obs_targ_name_sci':'HD_166191','obs_name_sci':'SCI_HD166191',
    'night_Lcal':'2019-09-18','tpl_start_Lcal':'2019-09-19T00:11:00',
    'obs_targ_name_Lcal':'HD_170016','obs_name_Lcal':'CAL_HD170016',
    'night_Ncal':'2019-09-18','tpl_start_Ncal':'2019-09-19T00:11:00',
    'obs_targ_name_Ncal':'HD_170016','obs_name_Ncal':'CAL_HD170016'},

    # {'night_sci':'2019-05-14','tpl_start_sci':'2019-05-15T00:10:56',
    # 'obs_targ_name_sci':'HD97048','obs_name_sci':'SCI_HD97048_b6ms',
    # 'night_Lcal':'','tpl_start_Lcal':'2019-05-15T01:36:00',
    # 'obs_targ_name_Lcal':'kap_Cha','obs_name_Lcal':'CAL_HD97048_kap_Cha_LN_b6ms',
    # 'night_Ncal':'','tpl_start_Ncal':'2019-05-15T01:36:00',
    # 'obs_targ_name_Ncal':'kap_Cha','obs_name_Ncal':'CAL_HD97048_kap_Cha_LN_b6ms'},

    # {'night_sci':'2019-05-13','tpl_start_sci':'2019-05-14T03:06:51',
    # 'obs_targ_name_sci':'RULup','obs_name_sci':'SCI_RULup',
    # 'night_Lcal':'2019-05-13','tpl_start_Lcal':'2019-05-14T03:39:23',
    # 'obs_targ_name_Lcal':'HD138492','obs_name_Lcal':'CAL_RULup_HD138492_LN',
    # 'night_Ncal':'2019-05-13','tpl_start_Ncal':'2019-05-14T03:39:23',
    # 'obs_targ_name_Ncal':'HD138492','obs_name_Ncal':'CAL_RULup_HD138492_LN'},

    # {'night_sci':'','tpl_start_sci':'',
    # 'obs_targ_name_sci':'','obs_name_sci':'',
    # 'night_Lcal':'','tpl_start_Lcal':'',
    # 'obs_targ_name_Lcal':'','obs_name_Lcal':'',
    # 'night_Ncal':'','tpl_start_Ncal':'',
    # 'obs_targ_name_Ncal':'','obs_name_Ncal':''},

    # {'night_sci':'2019-03-22','tpl_start_sci':'2019-03-23T08:41:19',
    # 'obs_targ_name_sci':'HD163296','obs_name_sci':'SCI_HD163296',
    # 'night_Lcal':'2019-03-22','tpl_start_Lcal':'2019-03-23T08:53:41',
    # 'obs_targ_name_Lcal':'delSgr','obs_name_Lcal':'CAL_HD163296_delSgr',
    # 'night_Ncal':'2019-03-22','tpl_start_Ncal':'2019-03-23T08:53:41',
    # 'obs_targ_name_Ncal':'delSgr','obs_name_Ncal':'CAL_HD163296_delSgr'},

    # {'night_sci':'2019-05-05','tpl_start_sci':'2019-05-06T08:19:51',
    # 'obs_targ_name_sci':'HD163296','obs_name_sci':'SCI_HD163296',
    # 'night_Lcal':'2019-05-05','tpl_start_Lcal':'2019-05-06T08:07:33',
    # 'obs_targ_name_Lcal':'HD_156637','obs_name_Lcal':'CAL_HD163296_HD156637_L',
    # 'night_Ncal':'2019-05-05','tpl_start_Ncal':'2019-05-06T08:38:14',
    # 'obs_targ_name_Ncal':'HD168454','obs_name_Ncal':'CAL_HD163296_HD168454_N'},
    ]


bands = ['L','N'] #,'N'] #'N', 'L'
for obs in obs_lst:
    night = obs['night_sci']
    obs_targ_name_sci = obs['obs_targ_name_sci']
    tpl_start_sci = obs['tpl_start_sci']
    tpl_start_Lcal = obs['tpl_start_Lcal']
    tpl_start_Ncal = obs['tpl_start_Ncal']
    for band in bands:
        print('Calibrating '+obs_targ_name_sci+', '+band+' band, TPL_START: '+tpl_start_sci)
        if band == 'L':
            tpl_sci = tpl_start_sci.replace(':','_')
            tpl_cal = tpl_start_Lcal.replace(':','_')
            tpldir_sci = '/'+tpl_sci+'/Iter1/mat_raw_estimates.'+tpl_sci+'.HAWAII-2RG.rb/'
            tpldir_cal = '/'+tpl_cal+'/Iter1/mat_raw_estimates.'+tpl_cal+'.HAWAII-2RG.rb/'
        if band == 'N':
            tpl_sci = tpl_start_sci.replace(':','_')
            tpl_cal = tpl_start_Ncal.replace(':','_')
            tpldir_sci = '/'+tpl_sci+'/Iter1/mat_raw_estimates.'+tpl_sci+'.AQUARIUS.rb/'
            tpldir_cal = '/'+tpl_cal+'/Iter1/mat_raw_estimates.'+tpl_cal+'.AQUARIUS.rb/'
        cohs= ['coherent','incoherent'] #['coherent','incoherent']
        cal_database_dir = basedir+ '/CalibMap/' #calib_spec_db.fits'
        if band == 'L':
            cal_database_paths = [cal_database_dir+'vBoekelDatabase.fits',cal_database_dir+'calib_spec_db.fits',cal_database_dir+'calib_spec_db_new.fits']
        if band == 'N':
            cal_database_paths = [cal_database_dir+'vBoekelDatabase.fits',cal_database_dir+'vBoekelDatabase.fitsold']

        for coh in cohs:
            #print(reduced_dir+coh+'/'+night+tpldir_sci)
            flst_sci = glob.glob(reduced_dir+coh+'/'+night+tpldir_sci+'/*_RAW_INT_*.fits')
            flst_sci = sorted(flst_sci)
            # flst_cal = glob.glob(reduced_dir+coh+'/'+night+tpldir_cal+'/*_RAW_INT_*.fits')
            # flst_cal = sorted(flst_cal)
            #print(flst_sci)
            if not flst_sci:
                print('WARNING: No reduced oifits files found in sci directory.')
            else:
                for i in range(len(flst_sci)):
                    inputfile_sci = flst_sci[i]
                    inputfile_cal = reduced_dir+coh+'/'+night+tpldir_cal+'/'+os.path.basename(flst_sci[i]).replace('TARGET','CALIB')

                    outputfile = reduced_dir+coh+'/'+night+tpldir_sci+'/TARGET_FLUXCAL_INT_'+os.path.basename(flst_sci[i])[15:19]+'.fits'
                    output_fig_dir = reduced_dir+coh+'/'+night+tpldir_sci
                    if coh == 'coherent':
                        mode = 'corrflux'
                    if coh == 'incoherent':
                        mode = 'flux'
                    if do_calibration:
                        fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_paths, mode=mode,output_fig_dir=output_fig_dir)

                inputdir = reduced_dir+coh+'/'+night+tpldir_sci
                if do_make_plots:
                    show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='TARGET_FLUXCAL_INT',nbProc=6)

                if do_copy == True:
                    for k in range(len(target_names)):
                        tags = target_tags[k].split('.')
                        for tag in tags:
                            if tag == obs_targ_name_sci:
                                #copy files in the destination directory  
                                dest_dir = targetdir + '/' + target_names[k]
                                # print(dest_dir)
                                band = ''
                                if 'HAWAII' in tpldir_sci:
                                    band = 'L'
                                if 'AQUARIUS' in tpldir_sci:
                                    band = 'N'
                                dest_filetag = target_names[k] + '_' + tpl_start_sci.replace(':','_') + '_' + band + '_' + coh
                                if not os.path.exists(dest_dir):
                                    os.mkdir(dest_dir)
                                #search for calibrated files:
                                cal_int_files= glob.glob(reduced_dir+coh+'/'+night+tpldir_sci+'/TARGET_FLUXCAL_INT_*.fits')
                                # print(tpldir_sci+'/TARGET_FLUXCAL_INT_*.fits')
                                # print(cal_int_files)
                                for file in cal_int_files:
                                    src_path = file
                                    dst_path = dest_dir+'/'+dest_filetag+'_'+os.path.basename(file)
                                    copyfile(src_path,dst_path)
                                #copy plots:
                                if do_copy_plots:
                                    cal_int_plots = glob.glob(reduced_dir+coh+'/'+night+tpldir_sci+'/plots/TARGET_FLUXCAL_INT_*_corrflux.png')
                                    for file in cal_int_plots:
                                        src_path = file
                                        dst_path = dest_dir+'/'+dest_filetag+'_'+os.path.basename(file)
                                        copyfile(src_path,dst_path)
                                    cal_int_plots = glob.glob(reduced_dir+coh+'/'+night+tpldir_sci+'/plots/TARGET_FLUXCAL_INT_*_corrflux_vs_pbl*')
                                    for file in cal_int_plots:
                                        src_path = file
                                        dst_path = dest_dir+'/'+dest_filetag+'_'+os.path.basename(file)
                                        copyfile(src_path,dst_path)
                                    cal_int_plots = glob.glob(reduced_dir+coh+'/'+night+tpldir_sci+'/plots/TARGET_FLUXCAL_INT_*_flux.png')
                                    for file in cal_int_plots:
                                        src_path = file
                                        dst_path = dest_dir+'/'+dest_filetag+'_'+os.path.basename(file)
                                        copyfile(src_path,dst_path)    
                                    uv_plots = glob.glob(reduced_dir+coh+'/'+night+tpldir_sci + '/plots/TARGET_FLUXCAL_INT_*_uv.png')
                                    for file in uv_plots:
                                        src_path = file
                                        dst_path = dest_dir + '/' + dest_filetag + '_uv.png'  
                                        copyfile(src_path, dst_path)

# cohs= ['coherent'] #['coherent','incoherent']
# for coh in cohs:
#     for i in range(0):
#         numstr = '%04d' % (i+1)
#         inputfile_sci = '/Users/jvarga/Data/MATISSE/matisse_red/'+coh+'/2019-03-22/mat_raw_estimates.2019-03-23T08_41_19.HAWAII-2RG.rb/TARGET_CAL_INT_'+numstr+'.fits'
#         inputfile_cal = '/Users/jvarga/Data/MATISSE/matisse_red/'+coh+'/2019-03-22/mat_raw_estimates.2019-03-23T08_53_41.HAWAII-2RG.rb/CALIB_RAW_INT_'+numstr+'.fits'
#         cal_database_file = '/Users/jvarga/Data/calib_spec_db.fits'
#         outputfile = '/Users/jvarga/Data/MATISSE/matisse_red/'+coh+'/2019-03-22/mat_raw_estimates.2019-03-23T08_41_19.HAWAII-2RG.rb/TARGET_FLUXCAL_INT_'+numstr+'.fits'
#
#         if coh == 'coherent':
#             mode = 'corrflux'
#         if coh == 'incoherent':
#             mode = 'flux'
#         fluxcal(inputfile_sci, inputfile_cal, outputfile, cal_database_file, mode=mode)
#
#     inputdir = '/Users/jvarga/Data/MATISSE/matisse_red/'+coh+'/2019-03-22/mat_raw_estimates.2019-03-23T08_41_19.HAWAII-2RG.rb/'
#     show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='TARGET_FLUXCAL', verbose=False,sel_wls=[3.2,3.6,4.0],bandwidths=[0.2,0.2,0.2])
#     #inputdir = '/Users/jvarga/Data/MATISSE/matisse_red/'+coh+'/2019-03-22/mat_raw_estimates.2019-03-23T08_53_41.HAWAII-2RG.rb/'
#     #show_allred(inputdir, outputdir=inputdir + '/plots/',fn_pattern='', verbose=False,sel_wls=[3.2,3.6,4.0],bandwidths=[0.2,0.2,0.2])
print('READY')

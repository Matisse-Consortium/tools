#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-
"""
Created on Mon May 25 13:51:41 2015
@author: bpichon
"""
import datetime
import numpy
import dislin
import sys
from matisse_visual_class import get_key_withdefault
#==============================================================================
def produce_oifits_report(oifits,filename,band):
    nwave = oifits.nwave_sc
#
# these parameters may be modified if necessary
#
# BEGIN BEGIN BEGIN of BLOCK
#
    if ( nwave < 500 ) :
        stride_err = 20
    elif ( nwave < 1000 ) :
        stride_err = 50
    else : 
        stride_err = 100
    fact_shift_err = max( stride_err / 20 , 1 )
    #
#    scale_flux = 1.0e+05
    scale_flux = 1.0e+04
    #
    uv_plane_max = 105.0
    uv_plane_step = 10.0
    #
    if ( band[0] == "LM" ) :
       tickx_wavelength = 0.10
       xmax_spatial_freq = 180.0
       tickx_spatial_freq = 10.0
    elif (band[0] == "N" ) :
       tickx_wavelength = 0.25
       xmax_spatial_freq = 65.0
       tickx_spatial_freq = 5.0
    else : 
       print( " IMPOSSIBLE : LM or N " )
       return
    #
    limit_min_vis = -0.1
    limit_max_vis = +1.1
    limit_err_vis = 0.2
    limit_min_phi = -180.0
    limit_max_phi = +180.0
    limit_err_phi = 720.0
    fact_err = 800.0
#
# END END END of BLOCK
# 
    device = "pdf"
#    device = "xwin"       # plus petit
#    device = "cons"       # plus grand
#    device = "virt"
    #
    if device == "pdf" :
        dislin.setfil(filename+".pdf")
        dislin.filmod("Delete")
        dislin.errdev("Append")
    #
    dislin.metafl(device)
    dislin.setpag("DA4P")
    dislin.disini()
#    dislin.pagera()
    if device == "pdf" :
        dislin.errmod("warnings","on")
        dislin.errmod("check","off")
        dislin.errmod("protocol","off")
#        dislin.errmod("protocol","file")

#    dislin.pagera()
    dislin.hwfont()
    dislin.texmod("ON")
    dislin.namdis(20,"XY")    # default 30
    dislin.labdig(-2,'X')
###
### DISLIN : 100 pts by cm 
###
    nxoffset_1 = 80
    nxoffset_2 = 1900
    nyoffset = 2920
    #
    #==============================================================================
    # Others parameters
    #==============================================================================
    #
    base_list=["12", "13", "14", "23", "24", "34"]
    nb_base = len(base_list)
    #
    triplet_list = [str(oifits.t3_sc_staindex[0,:]),str(oifits.t3_sc_staindex[1,:]),\
                    str(oifits.t3_sc_staindex[2,:]),str(oifits.t3_sc_staindex[3,:])]
    nb_triplet = len(triplet_list)
    #
    telescope_list = [ "Tel_1" , "Tel_2" , "Tel_3" , "Tel_4" ]
    nb_telescope = len(telescope_list)
    #
    minwave = numpy.nanmin(oifits.wave_sc) # minimum wavelength   # pourquoi nanmin et pas min 
    maxwave = numpy.nanmax(oifits.wave_sc) # maximum wavelength
    #
    if ( band[0] == "LM" ) :
        print ("Wavelength limits for the data (microns) from %.4f to %.4f "%(minwave,maxwave))
        waveSC_0 = band[1][0]
        minwaveSC_0 = waveSC_0[0]
        maxwaveSC_0 = waveSC_0[1]
        print ("Wavelength limits for band L in plots (microns) from %.4f to %.4f "%(minwaveSC_0,maxwaveSC_0))
        waveSC_1 = band[1][1]
        minwaveSC_1 = waveSC_1[0]
        maxwaveSC_1 = waveSC_1[1]
        print ("Wavelength limits for band M in plots (microns) from %.4f to %.4f "%(minwaveSC_1,maxwaveSC_1))
        minwaveSC = minwaveSC_0
        maxwaveSC = maxwaveSC_1
        waveok_0 = ( oifits.wave_sc>=minwaveSC_0 )
        waveok_1 = ( oifits.wave_sc<=maxwaveSC_0 )
        waveok_2 = ( oifits.wave_sc>=minwaveSC_1 )
        waveok_3 = ( oifits.wave_sc<=maxwaveSC_1 )
#        waveok = ( waveok_0 & waveok_3 )
        waveok = ( waveok_0 & waveok_1 ) | ( waveok_2 & waveok_3 ) 
#        print( " len of waveok = ", waveok.size )
#        print( len(oifits.wave_sc) , len(waveok) )
    elif ( band[0] == "N" ) :
        print ("Wavelength limits for the data (microns) from %.4f to %.4f "%(minwave,maxwave))
        waveSC = band[1][0]
        minwaveSC = waveSC[0]
        maxwaveSC = waveSC[1]
        print ("Wavelength limits for band N in plots (microns) from %.4f to %.4f "%(minwaveSC,maxwaveSC))
        waveok_0 = ( oifits.wave_sc>=minwaveSC )
        waveok_1 = ( oifits.wave_sc<=maxwaveSC )
        waveok = ( waveok_0 & waveok_1 )
#        print( " len of waveok = ", waveok.size )
#        print( len(oifits.wave_sc) , len(waveok) )
    else:
        print( " IMPOSSIBLE : LM or N " )
        return
    #
    #==============================================================================
    # Selection of wavelength range for plots
    #==============================================================================
    #
    oifits.wave_sc = oifits.wave_sc[waveok]
    oifits.nwave_sc = len(oifits.wave_sc)
#    print ( " len of arrays : ", oifits.nwave_sc, len(waveok) )
    oifits.spfreq_sc = oifits.spfreq_sc[:,waveok]
    oifits.t3_spfreq_sc = oifits.t3_spfreq_sc[:,waveok]
    oifits.t3phi_spfreq_sc = oifits.t3phi_spfreq_sc[:,waveok]
    #
    oifits.oi_flux_sc = oifits.oi_flux_sc[:,:,waveok]
    oifits.oi_vis2_sc_vis2data = oifits.oi_vis2_sc_vis2data[:,waveok]
    oifits.oi_vis2_sc_vis2err = oifits.oi_vis2_sc_vis2err[:,waveok]
    oifits.t3phi_sc = oifits.t3phi_sc[:,waveok]
    oifits.t3phierr_sc = oifits.t3phierr_sc[:,waveok]
    oifits.oi_vis_sc_visamp=oifits.oi_vis_sc_visamp[:,waveok]
    oifits.oi_vis_sc_visamperr=oifits.oi_vis_sc_visamperr[:,waveok]
    oifits.oi_vis_sc_visphi=oifits.oi_vis_sc_visphi[:,waveok]
    oifits.oi_vis_sc_visphierr=oifits.oi_vis_sc_visphierr[:,waveok]
    #
    nb_expo_0 = int((oifits.oi_vis2_sc_vis2data.shape)[0]) / 6
    nb_expo_1 = int((oifits.t3phi_sc.shape)[0]) / 4
    nb_expo_2 = int((oifits.oi_vis_sc_visamp.shape)[0]) / 6
    nb_expo_3 = int((oifits.oi_vis_sc_visphi.shape)[0]) / 6
    nb_expo_all = [ nb_expo_0, nb_expo_1, nb_expo_2, nb_expo_3 ]
    nb_expo_min = min( nb_expo_all )
    nb_expo_max = max( nb_expo_all )
    if ( nb_expo_min != nb_expo_max ) : 
        print( " BIG PROBLEM : ", nb_expo_all )
        return
    nb_expo = nb_expo_min
    print( "Number of exposures : nb_expo = %i "%nb_expo )
    #
    if ( nb_expo > 1 ) :
        b_color = 30
        a_color = (255-2*b_color) / (nb_expo-1)
    else :
        dislin.color("RED")
        b_color = dislin.getclr()
        dislin.color("FORE")
        a_color = 255             # any value
    #
    #==============================================================================
    # summary and start report
    #==============================================================================
    #
    if hasattr(oifits,'stations'):
        config = oifits.telnames[0]+'='+oifits.stations[0]+'  '+ \
                 oifits.telnames[1]+'='+oifits.stations[1]+'  '+ \
                 oifits.telnames[2]+'='+oifits.stations[2]+'  '+ \
                 oifits.telnames[3]+'='+oifits.stations[3]
    else:
        config = ''
    #
    lamp = band[2]
    if ( lamp ) :
        objname  = "SOURCE (Artificial source)" 
        objalpha = "     "          #BP FAUX BIEN SUR
        objdelta = "     "          #BP FAUX BIEN SUR
#        objmagL  = "-7.0"           #BP FAUX BIEN SUR
#        objmagM  = "-8.0"           #BP FAUX BIEN SUR
#        objmagN  = "-9.0"           #BP FAUX BIEN SUR
        objmagL  = "    "           #BP FAUX BIEN SUR
        objmagM  = "    "           #BP FAUX BIEN SUR
        objmagN  = "    "           #BP FAUX BIEN SUR
        ambi = ' not defined '
    else:
        object = str( get_key_withdefault(oifits.header,'OBJECT','Pichon star') )
        objname =  str( get_key_withdefault(oifits.header,'HIERARCH ESO OBS TARG NAME','Pichon star') )
        if ( object == "STD" ) : 
            objname += " [STD]"
        objalpha = str( get_key_withdefault(oifits.header,'RA','25.98') )
        objdelta = str( get_key_withdefault(oifits.header,'DEC','99.98') )
#        objmag = str(get_key_withdefault(product.header,'HIERARCH ESO FT ROBJ MAG',''))
#        objmagL  = "-7.1"           #BP FAUX BIEN SUR
#        objmagM  = "-8.1"           #BP FAUX BIEN SUR
#        objmagN  = "-9.1"           #BP FAUX BIEN SUR
        objmagL  = "TBD"           #BP FAUX BIEN SUR
        objmagM  = "TBD"           #BP FAUX BIEN SUR
        objmagN  = "TBD"           #BP FAUX BIEN SUR
        seeing_0 =  str( get_key_withdefault(oifits.header,'HIERARCH ESO ISS AMBI FWHM START','99.9') )
        seeing_1 =  str( get_key_withdefault(oifits.header,'HIERARCH ESO ISS AMBI FWHM END','99.9') )
        windspeed = str( get_key_withdefault(oifits.header,'HIERARCH ESO ISS AMBI WINDSP','99.99') )
        tau0V_0 =  str( get_key_withdefault(oifits.header,'HIERARCH ESO ISS AMBI TAU0 START','9.999') )
        tau0V_1 =  str( get_key_withdefault(oifits.header,'HIERARCH ESO ISS AMBI TAU0 END','9.999') )
        ambi = seeing_0 + " --> " + seeing_1 + "  ;  " + windspeed + "  ;  " + tau0V_0 + " --> " + tau0V_1
    #
    xoffset_1 = 0.02
    xoffset_2 = 0.45 
    xoffset_3 = xoffset_1 + xoffset_2
	#
    dislin.axspos(100,800)
    dislin.axslen(1900,650)
    dislin.titlin("MATISSE OIFITS Quality Control Report",4)

    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
#    dislin.grid(1,10)
    dislin.title()
#    dislin.curve( (xoffset_2,xoffset_2) , (-0.05,1.0) , 2 )
    dislin.height(24)
    
    yoffset = 0.95
    dislin.rlmess("Filename",xoffset_1,yoffset)
    dislin.rlmess(filename+".fits",xoffset_3,yoffset)
    
    yoffset = 0.85
    dislin.rlmess("Observing date",xoffset_1,yoffset)
    dislin.rlmess(oifits.header['DATE-OBS'],xoffset_3,yoffset)

    yoffset = 0.75
    dislin.rlmess("Processing/report date",xoffset_1,yoffset)
    dislin.rlmess(oifits.header['DATE']+'   '+datetime.datetime.now().strftime("%Y-%m-%dT%H:%M:%S"),xoffset_3,yoffset)
    
    yoffset = 0.65
    dislin.rlmess("Product category, Chip name",xoffset_1,yoffset)
    dislin.rlmess(oifits.header['HIERARCH ESO PRO CATG']+', '+oifits.header['HIERARCH ESO DET CHIP NAME'],xoffset_3,yoffset)

    yoffset = 0.55
    if band[0] == "LM":
        dislin.rlmess("DIL, PIL, POL, FIL, SFL, BCD1, BCD2",xoffset_1,yoffset)
        dislin.rlmess(oifits.header['HIERARCH ESO INS DIL ID']+', '+oifits.header['HIERARCH ESO INS PIL ID']+', '+
                    oifits.header['HIERARCH ESO INS POL ID']+', '+oifits.header['HIERARCH ESO INS FIL ID']+', '+
                    oifits.header['HIERARCH ESO INS SFL ID']+', '+oifits.header['HIERARCH ESO INS BCD1 ID']+', '+
                    oifits.header['HIERARCH ESO INS BCD2 ID'],xoffset_3,yoffset)
    elif band[0] == "N":
        dislin.rlmess("DIN, PIN, PON, FIN, SFN, BCD1, BCD2",xoffset_1,yoffset)
        dislin.rlmess(oifits.header['HIERARCH ESO INS DIN ID']+', '+oifits.header['HIERARCH ESO INS PIN ID']+', '+
                    oifits.header['HIERARCH ESO INS PON ID']+', '+oifits.header['HIERARCH ESO INS FIN ID']+', '+
                    oifits.header['HIERARCH ESO INS SFN ID']+', '+oifits.header['HIERARCH ESO INS BCD1 ID']+', '+
                    oifits.header['HIERARCH ESO INS BCD2 ID'],xoffset_3,yoffset)
    else:
        print ("IMPOSSIBLE !!!")
    #
    yoffset = 0.45
    ndit = oifits.header['HIERARCH ESO DET NDIT']
    dit = oifits.header['HIERARCH ESO DET SEQ1 DIT']
    time_tot = ndit * dit
    dislin.rlmess("NDIT x DIT ;  time_tot ;  nb_expo ;  nwave",xoffset_1,yoffset)
    dislin.rlmess( str(ndit) + " x " + str(dit) + " s  ;   " + str(time_tot) + " s  ;   " + str(nb_expo) + "  ;   " + str(nwave),xoffset_3,yoffset)
    #
    yoffset = 0.35
    dislin.rlmess("Object name",xoffset_1,yoffset)
    dislin.rlmess(objname,xoffset_3,yoffset)
    #
    yoffset = 0.25
    if band[0] == "LM":
        dislin.rlmess("Object RA, Dec, L, M",xoffset_1,yoffset)
        dislin.rlmess(objalpha +'  '+ objdelta +'  L = '+ objmagL +'  M = '+ objmagM,xoffset_3,yoffset)
    elif band[0] == "N":
        dislin.rlmess("Object RA, Dec, N",xoffset_1,yoffset)
        dislin.rlmess(objalpha +'  '+ objdelta +'  N = '+ objmagN,xoffset_3,yoffset)
    else:
        print ("IMPOSSIBLE !!!")
    #
    yoffset = 0.15
    dislin.rlmess("Telescope stations",xoffset_1,yoffset)
    dislin.rlmess(config,xoffset_3,yoffset)
    #
    yoffset = 0.05
    dislin.rlmess("Seeing (arcsec) ; Wind (m/s) ; T0 in V (s) ",xoffset_1,yoffset)
    dislin.rlmess(ambi,xoffset_3,yoffset)
    #
    dislin.sendbf()
    dislin.endgrf()
    #
    if (  not lamp ) : 
##    #==============================================================================
##    # UV-plane
##    #==============================================================================
        dislin.setgrf("name","name","ticks","ticks")    # restore default
        dislin.axends("none","XY")                      # restore defaults
        dislin.ticks(1,"XY")
        dislin.axspos(100,1200)
        dislin.axslen(1900,350)
    #    dislin.nograf()
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
        #
        dislin.height(50)
        yoffset = 0.85
        xoffset_1 = 0.05
        dislin.rlmess("expo ==> color",xoffset_1,yoffset)
        yoffset_1 = 0.3
        yoffset_2 = 0.5
        xoffset = 0.0
        for expo in range(nb_expo) :
            xoffset = xoffset + 0.05
            n_color = a_color * expo + b_color
            dislin.setclr(n_color)
            dislin.rlmess(str(expo),xoffset,yoffset_1)
            nsymb = 21
            dislin.rlsymb(nsymb,xoffset,yoffset_2)
        #
        dislin.color("FORE")
        dislin.sendbf()
        dislin.endgrf()
        #
        dislin.reset("nograf")
        dislin.setgrf("name","name","ticks","ticks")    # restore default
        dislin.axends("none","XY")                      # restore defaults
        dislin.ticks(1,"XY")
        dislin.height(30)
        #
        dislin.axspos(300,2800)
        dislin.axslen(1500,1500)
        #
        dislin.setgrf("name","name","labels","labels")
        uv_max = uv_plane_max
        uv_min = - uv_max
        uv_ori = uv_max
        if ( int(uv_ori)%10 != 0 ) :
            uv_ori = 10.0 * float( int(uv_ori)/10 )
        uv_step = uv_plane_step
        dislin.height(24)
        dislin.graf( uv_min,uv_max,-uv_ori,uv_step , uv_min,uv_max,-uv_ori,uv_step )
        dislin.dot()
        dislin.grid(1,1)
        dislin.solid()
        #
#        dislin.height(35)
#        nx_uv = 100
#        ny_uv = 2425
#        dislin.messag( "UV Plane (in m) ==> ",nx_uv,ny_uv )
        #
        u_coords = oifits.ucoord_sc
        v_coords = oifits.vcoord_sc
        dislin.height(30)
        nb_coords =  len(u_coords) 
        if ( nb_base * nb_expo != nb_coords ) :
            print( " impossible : ", nb_base, nb_expo, nb_coords )
            return
        #
        for base in range(nb_base) :
            for expo in range(nb_expo) :
                n = base + expo * nb_base
                n_color = a_color * expo + b_color
                dislin.setclr(n_color)
                nsymb = 8
                dislin.rlsymb( nsymb ,+u_coords[n] ,+v_coords[n] )
                dislin.rlsymb( nsymb ,-u_coords[n] ,-v_coords[n] )
#        for n in range(nb_coords):
#            dislin.rlsymb( 21 ,+u_coords[n] ,+v_coords[n] )
#            dislin.rlsymb( 21 ,-u_coords[n] ,-v_coords[n] )
    #
    dislin.color("FORE")
    dislin.height(30)
    page = 1
    print ( " END OF PAGE ", page )
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    dislin.sendbf()
    dislin.endgrf()
    dislin.newpag()
#
#    #==============================================================================
#    # Tables of numerical values									PAGES 2-x et 2
#    #==============================================================================
#
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    values_vis2    = numpy.empty( (nb_expo,nb_base,7) )
    values_visamp  = numpy.empty( (nb_expo,nb_base,7) )
    values_visphi  = numpy.empty( (nb_expo,nb_base,7) )
    values_t3phi = numpy.empty( (nb_expo,nb_triplet,2) )
    #
    # LOG
    #
#    log = open("mean_log.txt","a")
#    string = " fact_err = %.2f"%(fact_err)
#    log.write(string+"\n")
    #
    les_posy = [0.54 , 0.44 , 0.34 , 0.24 , 0.14 , 0.04 ]
    les_xoffset_base = [0.02 , 0.16 , 0.42 , 0.54 , 0.66 , 0.78 , 0.90 ]
    #
    for expo in range(nb_expo) :
        #
        dislin.height(36)
        dislin.axspos(150,600)
        dislin.axslen(1900,350)
        string = "Exposure number %i"%(expo)
        dislin.titlin("",1)
        dislin.titlin("",2)
        dislin.titlin(string,3)
        dislin.titlin("",4)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.nograf()
        dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
        #dislin.height(24)
        dislin.title()
        dislin.reset("nograf")
        dislin.sendbf()
        dislin.endgrf()
        #
        n_color = a_color * expo + b_color
        dislin.setclr(n_color)
        dislin.shdpat(16)
        dislin.rectan(0135,70,700,65)
        dislin.rectan(1400,70,600,65)
        dislin.color("FORE")
        #
        # vis^2
        #
        dislin.axspos(100,800)
        dislin.axslen(1900,650)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
        dislin.height(30)
        xoffset_1 = les_xoffset_base[0]
        yoffset = 0.95
        dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
        yoffset = 0.88
        dislin.rlmess("Col 2 : Average squared visibility per baseline (vis^2 $\pm$ err) ==> page 3",xoffset_1,yoffset)
        yoffset = 0.81
        dislin.rlmess("Cols 3 --> 7 : Fraction of points Ok , points with value<limit_min , value>limit_max ",xoffset_1,yoffset)
        yoffset = 0.74
        dislin.rlmess("                            points with error(err)>limit_err , error(tol)>limit_tol ",xoffset_1,yoffset)
        #
        dislin.height(34)
        yoffset = 0.64
        dislin.rlmess("Baseline",les_xoffset_base[0],yoffset)
        dislin.rlmess("vis^2   ",les_xoffset_base[1],yoffset)
        dislin.rlmess("frac_ok ",les_xoffset_base[2],yoffset)
        dislin.rlmess("frac_min",les_xoffset_base[3],yoffset)
        dislin.rlmess("frac_max",les_xoffset_base[4],yoffset)
        dislin.rlmess("frac_err",les_xoffset_base[5],yoffset)
        dislin.rlmess("frac_tol",les_xoffset_base[6],yoffset)
        #
        for base in range(nb_base) :
            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
        #
        dislin.height(36)
        dislin.fixspc(0.50)
        for base in range(nb_base) :
#            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
            index = base + expo * nb_base
            yy = [ float(elem) for elem in oifits.oi_vis2_sc_vis2data[index,:] ]
            err = [ float(elem) for elem in oifits.oi_vis2_sc_vis2err[index,:] ]
            (mean,incer,n_limits) = my_mean( yy , err , limit_err_vis , limit_min_vis , limit_max_vis , fact_err )
            n_tot = len(yy)         #                      sum(n_limits[0:6])
            n_ok = n_limits[0]
            n_min = n_limits[1]
            n_max = n_limits[2]
            n_err = n_limits[3]
            n_tol = n_limits[4]
            frac_ok = float(n_ok)/float(n_tot)
            frac_min = float(n_min)/float(n_tot)
            frac_max = float(n_max)/float(n_tot)
            frac_err = float(n_err)/float(n_tot)
            frac_tol = float(n_tol)/float(n_tot)
            string = "%.3f"%(mean)+" $\pm$ "+"%.3f"%(incer)
            dislin.rlmess(string,les_xoffset_base[1],les_posy[base])
            string = "%.3f"%(frac_ok)
            dislin.rlmess(string,les_xoffset_base[2],les_posy[base])
            string = "%.3f"%(frac_min)
            dislin.rlmess(string,les_xoffset_base[3],les_posy[base])
            string = "%.3f"%(frac_max)
            dislin.rlmess(string,les_xoffset_base[4],les_posy[base])
            string = "%.3f"%(frac_err)
            dislin.rlmess(string,les_xoffset_base[5],les_posy[base])
            string = "%.3f"%(frac_tol)
            dislin.rlmess(string,les_xoffset_base[6],les_posy[base])
            values_vis2[expo,base,0] = mean
            values_vis2[expo,base,1] = incer
            values_vis2[expo,base,2] = frac_ok
            values_vis2[expo,base,3] = frac_min
            values_vis2[expo,base,4] = frac_max
            values_vis2[expo,base,5] = frac_err
            values_vis2[expo,base,6] = frac_tol
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # vis
        #
        dislin.axspos(100,1500)
        dislin.axslen(1900,650)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
		#
        xoffset_1 = 0.02
        dislin.height(30)
        yoffset = 0.95
        dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
        yoffset = 0.88
        dislin.rlmess("Col 2 : Average visibility amplitude per baseline (vis $\pm$ err) ==> page 4",xoffset_1,yoffset)
        yoffset = 0.81
        dislin.rlmess("Cols 3 --> 7 : Fraction of points Ok , points with value<limit_min , value>limit_max ",xoffset_1,yoffset)
        yoffset = 0.74
        dislin.rlmess("                            points with error(err)>limit_err , error(tol)>limit_tol ",xoffset_1,yoffset)
        #
        dislin.height(34)
        yoffset = 0.64
        dislin.rlmess("Baseline",les_xoffset_base[0],yoffset)
        dislin.rlmess("vis     ",les_xoffset_base[1],yoffset)
        dislin.rlmess("frac_ok ",les_xoffset_base[2],yoffset)
        dislin.rlmess("frac_min",les_xoffset_base[3],yoffset)
        dislin.rlmess("frac_max",les_xoffset_base[4],yoffset)
        dislin.rlmess("frac_err",les_xoffset_base[5],yoffset)
        dislin.rlmess("frac_tol",les_xoffset_base[6],yoffset)
        #
        for base in range(nb_base) :
            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
        #
        dislin.height(36)
        dislin.fixspc(0.50)
        for base in range(nb_base) :
#            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
            index = base + expo * nb_base
            yy = [ float(elem) for elem in oifits.oi_vis_sc_visamp[index,:] ]
            err = [ float(elem) for elem in oifits.oi_vis_sc_visamperr[index,:] ]
            (mean,incer,n_limits) = my_mean( yy , err , limit_err_vis , limit_min_vis , limit_max_vis , fact_err )
            n_tot = len(yy)         #                     = sum(n_limits[0:6])
            n_ok = n_limits[0]
            n_min = n_limits[1]
            n_max = n_limits[2]
            n_err = n_limits[3]
            n_tol = n_limits[4]
            frac_ok = float(n_ok)/float(n_tot)
            frac_min = float(n_min)/float(n_tot)
            frac_max = float(n_max)/float(n_tot)
            frac_err = float(n_err)/float(n_tot)
            frac_tol = float(n_tol)/float(n_tot)
            string = "%.3f"%(mean)+" $\pm$ "+"%.3f"%(incer)
            dislin.rlmess(string,les_xoffset_base[1],les_posy[base])
            string = "%.3f"%(frac_ok)
            dislin.rlmess(string,les_xoffset_base[2],les_posy[base])
            string = "%.3f"%(frac_min)
            dislin.rlmess(string,les_xoffset_base[3],les_posy[base])
            string = "%.3f"%(frac_max)
            dislin.rlmess(string,les_xoffset_base[4],les_posy[base])
            string = "%.3f"%(frac_err)
            dislin.rlmess(string,les_xoffset_base[5],les_posy[base])
            string = "%.3f"%(frac_tol)
            dislin.rlmess(string,les_xoffset_base[6],les_posy[base])
            values_visamp[expo,base,0] = mean
            values_visamp[expo,base,1] = incer
            values_visamp[expo,base,2] = frac_ok
            values_visamp[expo,base,3] = frac_min
            values_visamp[expo,base,4] = frac_max
            values_visamp[expo,base,5] = frac_err
            values_visamp[expo,base,6] = frac_tol
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # vis_phi
        #
        dislin.axspos(100,2200)
        dislin.axslen(1900,650)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
    
        xoffset_1 = 0.02
        dislin.height(30)
        yoffset = 0.95
        dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
        yoffset = 0.88
        dislin.rlmess("Col 2 : Average differential phase per baseline (visphi $\pm$ err), in degrees ==> page 6",xoffset_1,yoffset)
        yoffset = 0.81
        dislin.rlmess("Cols 3 --> 7 : Fraction of points Ok , points with value<limit_min , value>limit_max ",xoffset_1,yoffset)
        yoffset = 0.74
        dislin.rlmess("                            points with error(err)>limit_err , error(tol)>limit_tol ",xoffset_1,yoffset)
        #
        dislin.height(34)
        yoffset = 0.64
        dislin.rlmess("Baseline",les_xoffset_base[0],yoffset)
        dislin.rlmess("vis_phi ",les_xoffset_base[1],yoffset)
        dislin.rlmess("frac_ok ",les_xoffset_base[2],yoffset)
        dislin.rlmess("frac_min",les_xoffset_base[3],yoffset)
        dislin.rlmess("frac_max",les_xoffset_base[4],yoffset)
        dislin.rlmess("frac_err",les_xoffset_base[5],yoffset)
        dislin.rlmess("frac_tol",les_xoffset_base[6],yoffset)
        #
        for base in range(nb_base) :
            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
        #
        dislin.height(36)
        dislin.fixspc(0.50)
        for base in range(nb_base) :
#            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
            index = base + expo * nb_base
            yy = [ float(elem) for elem in oifits.oi_vis_sc_visphi[index,:] ]
            err = [ float(elem) for elem in oifits.oi_vis_sc_visphierr[index,:] ]
            (mean,incer,n_limits) = my_mean( yy , err , limit_err_phi , limit_min_phi , limit_max_phi , fact_err )
            n_tot = len(yy)         #                     = sum(n_limits[0:6])
            n_ok = n_limits[0]
            n_min = n_limits[1]
            n_max = n_limits[2]
            n_err = n_limits[3]
            n_tol = n_limits[4]
            frac_ok = float(n_ok)/float(n_tot)
            frac_min = float(n_min)/float(n_tot)
            frac_max = float(n_max)/float(n_tot)
            frac_err = float(n_err)/float(n_tot)
            frac_tol = float(n_tol)/float(n_tot)
            string = "%+8.3f"%(mean)+" $\pm$ "+"%7.3f"%(incer)
            dislin.rlmess(string,les_xoffset_base[1]-0.05,les_posy[base])
            string = "%.3f"%(frac_ok)
            dislin.rlmess(string,les_xoffset_base[2],les_posy[base])
            string = "%.3f"%(frac_min)
            dislin.rlmess(string,les_xoffset_base[3],les_posy[base])
            string = "%.3f"%(frac_max)
            dislin.rlmess(string,les_xoffset_base[4],les_posy[base])
            string = "%.3f"%(frac_err)
            dislin.rlmess(string,les_xoffset_base[5],les_posy[base])
            string = "%.3f"%(frac_tol)
            dislin.rlmess(string,les_xoffset_base[6],les_posy[base])
            values_visphi[expo,base,0] = mean
            values_visphi[expo,base,1] = incer
            values_visphi[expo,base,2] = frac_ok
            values_visphi[expo,base,3] = frac_min
            values_visphi[expo,base,4] = frac_max
            values_visphi[expo,base,5] = frac_err
            values_visphi[expo,base,6] = frac_tol
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # average closure phase
        #
        dislin.height(24)
        #
        dislin.reset("nograf")
        dislin.setgrf("name","name","ticks","ticks")    # restore default
        dislin.axends("none","XY")                      # restore defaults
        dislin.ticks(1,"XY")
        #
        dislin.axspos(100,2500)
        dislin.axslen(1900,250)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
        #
        dislin.height(30)
        yoffset = 0.90
        dislin.rlmess("Average closure phase per triplet (t3phi $\pm$ err), in degrees ==> page 5",xoffset_1,yoffset)
    
        dislin.height(34)
        xoffset = 0.02
        yoffset_1 = 0.55
        yoffset_2 = 0.20
        dislin.rlmess("Triplet",xoffset,yoffset_1)
        dislin.rlmess("Phi(deg)",xoffset,yoffset_2)
        les_posx = [0.16 , 0.37 , 0.58 , 0.80 ]
        dislin.fixspc(0.5)
        for triplet in range(nb_triplet) :
            dislin.rlmess(triplet_list[triplet],les_posx[triplet],yoffset_1)
            index = triplet + expo * nb_triplet
            yy = [ float(elem) for elem in oifits.t3phi_sc[index,:] ]
            err = [ float(elem) for elem in oifits.t3phierr_sc[index,:] ]
            (mean,incer,_) = my_mean( yy , err , limit_err_phi , limit_min_phi , limit_max_phi , fact_err )
            string = "%+.3f"%(mean)+" $\pm$ "+"%.3f"%(incer)
            dislin.rlmess(string,les_posx[triplet],yoffset_2)
            values_t3phi[expo,triplet,0] = mean
            values_t3phi[expo,triplet,1] = incer
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        if ( nb_expo == 1 ) :
            #
            # average photometric flux
            #
            dislin.reset("nograf")
            dislin.setgrf("name","name","ticks","ticks")    # restore default
            dislin.axends("none","XY")                      # restore defaults
            dislin.ticks(1,"XY")
            #
            dislin.axspos(100,2800)
            dislin.axslen(1900,250)
            dislin.setgrf("none","none","none","none")
            dislin.axends("noends","XY")
            dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
            #
            dislin.height(36)
            xoffset = 0.02
            yoffset = 0.92
            dislin.rlmess("Average photometric flux ("+"%.1e"%(scale_flux)+" photo-e-/s/sp.channel $\pm$ std) ==> page 7",xoffset,yoffset)
            #
            dislin.height(36)
            xoffset = 0.02
            yoffset_1 = 0.55
            yoffset_2 = 0.20
            y_scale = 1.0 / scale_flux
            dislin.rlmess("Telescope",xoffset,yoffset_1)
            dislin.rlmess("Flux",xoffset,yoffset_2)
            #
            les_posx = [0.20 , 0.40 , 0.60 , 0.80 ]
            dislin.fixspc(0.50)
            for tel in range(nb_telescope) :
                dislin.rlmess(telescope_list[tel],les_posx[tel],yoffset_1)
                index = tel + expo * nb_telescope
#               yy = [ float(elem) for elem in oifits.oi_flux_sc[:,index,:] ]
#               err = [ float(elem) for elem in oifits.oi_flux_err_sc[:,index,:] ]
#               (mean,incer,_) = my_mean( yy , err , limit_err_phi , limit_min_phi , limit_max_phi , fact_err )
#               string = "%+.3f"%(mean)+" $\pm$ "+"%.3f"%(incer)
#               dislin.rlmess(string,les_posx[triplet],yoffset_2)
#               z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,index,:])
#               z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,index,:])
                z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,tel,:])
                z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,tel,:])
                dislin.rlmess("%.3f"%(z1)+" $\pm$ "+"%.3f"%(z2),les_posx[tel],yoffset_2)
            #
            dislin.reset("fixspc")
            dislin.sendbf()
            dislin.endgrf()
        #
        dislin.color("FORE")
        dislin.height(30)
        page = 2
        if ( nb_expo == 1 ) : 
            print ( " END OF PAGE ", page )
            dislin.messag(filename,nxoffset_1,nyoffset)
            dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
        else : 
            print ( " END OF PAGE ", page , expo)
            dislin.messag(filename,nxoffset_1,nyoffset)
            dislin.messag("Page : "+str(page)+"-"+str(expo),nxoffset_2-40,nyoffset)
        #
        dislin.sendbf()
        dislin.endgrf()
        dislin.newpag()
    #
    # SUMMARY OF ALL EXPOSURES
    #
    if ( nb_expo > 1 ):
        dislin.height(36)
        dislin.axspos(150,600)
        dislin.axslen(1900,350)
        string = "Summary of all exposures"
        dislin.titlin("",1)   # to erase previous call
        dislin.titlin("",2)
        dislin.titlin(string,3)   # to erase previous call
        dislin.titlin("",4)   # to erase previous call
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.nograf()
        dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
        dislin.title()
        dislin.reset("nograf")
        dislin.sendbf()
        dislin.endgrf()
		#
        dislin.axspos(100,800)
        dislin.axslen(1900,650)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
		#
        xoffset_1 = 0.02
        dislin.height(30)
        yoffset = 0.95
        dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
        yoffset = 0.88
        dislin.rlmess("Col 2 : Average squared visibility per baseline (vis^2 $\pm$ err) ==> page 3",xoffset_1,yoffset)
        yoffset = 0.81
        dislin.rlmess("Cols 3 --> 7 : Fraction of points Ok , points with value<limit_min , value>limit_max ",xoffset_1,yoffset)
        yoffset = 0.74
        dislin.rlmess("                            points with error(err)>limit_err , error(tol)>limit_tol ",xoffset_1,yoffset)
        #
        dislin.height(34)
        yoffset = 0.64
        dislin.rlmess("Baseline",les_xoffset_base[0],yoffset)
        dislin.rlmess("vis^2"   ,les_xoffset_base[1],yoffset)
        dislin.rlmess("frac_ok ",les_xoffset_base[2],yoffset)
        dislin.rlmess("frac_min",les_xoffset_base[3],yoffset)
        dislin.rlmess("frac_max",les_xoffset_base[4],yoffset)
        dislin.rlmess("frac_err",les_xoffset_base[5],yoffset)
        dislin.rlmess("frac_tol",les_xoffset_base[6],yoffset)
        #
        dislin.height(36)
        dislin.fixspc(0.50)
        for base in range(nb_base) :
            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
            mean     = numpy.mean( values_vis2[:,base,0] )
            std      = numpy.std( values_vis2[:,base,0] )
            incer    = numpy.mean( values_vis2[:,base,1] )
            frac_ok  = numpy.mean( values_vis2[:,base,2] )
            frac_min = numpy.mean( values_vis2[:,base,3] )
            frac_max = numpy.mean( values_vis2[:,base,4] )
            frac_err = numpy.mean( values_vis2[:,base,5] )
            frac_tol = numpy.mean( values_vis2[:,base,6] )
            string = "%.3f"%(mean)+" $\pm$ "+"%.3f"%(std)+" $\pm$ "+"%.3f"%(incer)
            dislin.rlmess(string,les_xoffset_base[1]-0.06,les_posy[base])
            string = "%.3f"%(frac_ok)
            dislin.rlmess(string,les_xoffset_base[2],les_posy[base])
            string = "%.3f"%(frac_min)
            dislin.rlmess(string,les_xoffset_base[3],les_posy[base])
            string = "%.3f"%(frac_max)
            dislin.rlmess(string,les_xoffset_base[4],les_posy[base])
            string = "%.3f"%(frac_err)
            dislin.rlmess(string,les_xoffset_base[5],les_posy[base])
            string = "%.3f"%(frac_tol)
            dislin.rlmess(string,les_xoffset_base[6],les_posy[base])
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # vis
        #
        dislin.axspos(100,1500)
        dislin.axslen(1900,650)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
		#
        xoffset_1 = 0.02
        dislin.height(30)
        yoffset = 0.95
        dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
        yoffset = 0.88
        dislin.rlmess("Col 2 : Average visibility amplitude per baseline (vis $\pm$ err) ==> page 4",xoffset_1,yoffset)
        yoffset = 0.81
        dislin.rlmess("Cols 3 --> 7 : Fraction of points Ok , points with value<limit_min , value>limit_max ",xoffset_1,yoffset)
        yoffset = 0.74
        dislin.rlmess("                            points with error(err)>limit_err , error(tol)>limit_tol ",xoffset_1,yoffset)
        #
        dislin.height(34)
        yoffset = 0.64
        dislin.rlmess("Baseline",les_xoffset_base[0],yoffset)
        dislin.rlmess("vis     ",les_xoffset_base[1],yoffset)
        dislin.rlmess("frac_ok ",les_xoffset_base[2],yoffset)
        dislin.rlmess("frac_min",les_xoffset_base[3],yoffset)
        dislin.rlmess("frac_max",les_xoffset_base[4],yoffset)
        dislin.rlmess("frac_err",les_xoffset_base[5],yoffset)
        dislin.rlmess("frac_tol",les_xoffset_base[6],yoffset)
        #
        for base in range(nb_base) :
            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
        #
        dislin.height(36)
        dislin.fixspc(0.50)
        for base in range(nb_base) :
#            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
            mean     = numpy.mean( values_visamp[:,base,0] )
            std      = numpy.std( values_visamp[:,base,0] )
            incer    = numpy.mean( values_visamp[:,base,1] )
            frac_ok  = numpy.mean( values_visamp[:,base,2] )
            frac_min = numpy.mean( values_visamp[:,base,3] )
            frac_max = numpy.mean( values_visamp[:,base,4] )
            frac_err = numpy.mean( values_visamp[:,base,5] )
            frac_tol = numpy.mean( values_visamp[:,base,6] )
            string = "%.3f"%(mean)+" $\pm$ "+"%.3f"%(std)+" $\pm$ "+"%.3f"%(incer)
            dislin.rlmess(string,les_xoffset_base[1]-0.06,les_posy[base])
            string = "%.3f"%(frac_ok)
            dislin.rlmess(string,les_xoffset_base[2],les_posy[base])
            string = "%.3f"%(frac_min)
            dislin.rlmess(string,les_xoffset_base[3],les_posy[base])
            string = "%.3f"%(frac_max)
            dislin.rlmess(string,les_xoffset_base[4],les_posy[base])
            string = "%.3f"%(frac_err)
            dislin.rlmess(string,les_xoffset_base[5],les_posy[base])
            string = "%.3f"%(frac_tol)
            dislin.rlmess(string,les_xoffset_base[6],les_posy[base])
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # visphi
        #
        dislin.axspos(100,2200)
        dislin.axslen(1900,650)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
		#
        xoffset_1 = 0.02
        dislin.height(30)
        yoffset = 0.95
        dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
        yoffset = 0.88
        dislin.rlmess("Col 2 : Average differential phase per baseline (visphi $\pm$ err), in degrees ==> page 6",xoffset_1,yoffset)
        yoffset = 0.81
        dislin.rlmess("Cols 3 --> 7 : Fraction of points Ok , points with value<limit_min , value>limit_max ",xoffset_1,yoffset)
        yoffset = 0.74
        dislin.rlmess("                            points with error(err)>limit_err , error(tol)>limit_tol ",xoffset_1,yoffset)
        #
        dislin.height(34)
        yoffset = 0.64
        dislin.rlmess("Baseline",les_xoffset_base[0],yoffset)
        dislin.rlmess("vis_phi ",les_xoffset_base[1],yoffset)
        dislin.rlmess("frac_ok ",les_xoffset_base[2],yoffset)
        dislin.rlmess("frac_min",les_xoffset_base[3],yoffset)
        dislin.rlmess("frac_max",les_xoffset_base[4],yoffset)
        dislin.rlmess("frac_err",les_xoffset_base[5],yoffset)
        dislin.rlmess("frac_tol",les_xoffset_base[6],yoffset)
        #
        for base in range(nb_base) :
            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
        #
        dislin.height(36)
        dislin.fixspc(0.50)
        for base in range(nb_base) :
#            dislin.rlmess(base_list[base],les_xoffset_base[0],les_posy[base])
            mean     = numpy.mean( values_visphi[:,base,0] )
            std      = numpy.std( values_visphi[:,base,0] )
            incer    = numpy.mean( values_visphi[:,base,1] )
            frac_ok  = numpy.mean( values_visphi[:,base,2] )
            frac_min = numpy.mean( values_visphi[:,base,3] )
            frac_max = numpy.mean( values_visphi[:,base,4] )
            frac_err = numpy.mean( values_visphi[:,base,5] )
            frac_tol = numpy.mean( values_visphi[:,base,6] )
            string = "%+7.3f"%(mean)+" $\pm$ "+"%.3f"%(std)+" $\pm$ "+"%.3f"%(incer)
            dislin.rlmess(string,les_xoffset_base[1]-0.10,les_posy[base])
            string = "%.3f"%(frac_ok)
            dislin.rlmess(string,les_xoffset_base[2],les_posy[base])
            string = "%.3f"%(frac_min)
            dislin.rlmess(string,les_xoffset_base[3],les_posy[base])
            string = "%.3f"%(frac_max)
            dislin.rlmess(string,les_xoffset_base[4],les_posy[base])
            string = "%.3f"%(frac_err)
            dislin.rlmess(string,les_xoffset_base[5],les_posy[base])
            string = "%.3f"%(frac_tol)
            dislin.rlmess(string,les_xoffset_base[6],les_posy[base])
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # average closure phase
        #
        dislin.height(24)
        #
        dislin.reset("nograf")
        dislin.setgrf("name","name","ticks","ticks")    # restore default
        dislin.axends("none","XY")                      # restore defaults
        dislin.ticks(1,"XY")
        #
        dislin.axspos(100,2500)
        dislin.axslen(1900,250)
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
        #
        dislin.height(30)
        yoffset = 0.90
        dislin.rlmess("Average closure phase per triplet (t3phi $\pm$ err), in degrees ==> page 5",xoffset_1,yoffset)
        #
        dislin.height(34)
        xoffset = 0.02
        yoffset_1 = 0.65
        yoffset_2 = 0.30
        dislin.rlmess("Triplet",xoffset,yoffset_1)
        dislin.rlmess("Phi(deg)",xoffset,yoffset_2)
        les_posx = [0.15 , 0.33 , 0.51 , 0.69 ]
        #
        for triplet in range(nb_triplet) :
            dislin.rlmess(triplet_list[triplet],les_posx[triplet],yoffset_1)
        #
        dislin.fixspc(0.5)
        yoffset_delta_0 = 0.1
        for triplet in range(nb_triplet) :
#            dislin.rlmess(triplet_list[triplet],les_posx[triplet],yoffset_1)
            mean     = numpy.mean( values_t3phi[:,triplet,0] )
            std      = numpy.std( values_t3phi[:,triplet,0] )
            incer    = numpy.mean( values_t3phi[:,triplet,1] )
            string = "%+.3f"%(mean)+" $\pm$ "+"%.3f"%(std)+" $\pm$ "+"%.3f"%(incer)
            yoffset_delta = yoffset_delta_0 * float(1-2*(triplet%2))
            dislin.rlmess(string,les_posx[triplet],yoffset_2+yoffset_delta)
            values_t3phi[expo,triplet,0] = mean
            values_t3phi[expo,triplet,1] = incer
        #
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
        #
        # average photometric flux
        #
        dislin.reset("nograf")
        dislin.setgrf("name","name","ticks","ticks")    # restore default
        dislin.axends("none","XY")                      # restore defaults
        dislin.ticks(1,"XY")
        #
        dislin.axspos(100,2800)
        dislin.axslen(1900,250)
        #
        dislin.setgrf("none","none","none","none")
        dislin.axends("noends","XY")
        dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
        #
        dislin.height(36)
        xoffset = 0.02
        yoffset = 0.92
        dislin.rlmess("Average photometric flux ("+"%.1e"%(scale_flux)+" photo-e-/s/sp.channel $\pm$ std) ==> page 7",xoffset,yoffset)
        #
        dislin.height(36)
        xoffset = 0.02
        yoffset_1 = 0.55
        yoffset_2 = 0.20
        y_scale = 1.0 / scale_flux
        dislin.rlmess("Telescope",xoffset,yoffset_1)
        dislin.rlmess("Flux",xoffset,yoffset_2)
        #
        les_posx = [0.20 , 0.40 , 0.60 , 0.80 ]
        dislin.height(36)
        dislin.fixspc(0.5)
        for tel in range(nb_telescope) :
            dislin.rlmess(telescope_list[tel],les_posx[tel],yoffset_1)
#            index = tel + expo * nb_telescope
#            yy = [ float(elem) for elem in oifits.oi_flux_sc[:,index,:] ]
#            err = [ float(elem) for elem in oifits.oi_flux_err_sc[:,index,:] ]
#            (mean,incer,_) = my_mean( yy , err , limit_err_phi , limit_min_phi , limit_max_phi , fact_err )
#            string = "%+.3f"%(mean)+" $\pm$ "+"%.3f"%(incer)
#            dislin.rlmess(string,les_posx[triplet],yoffset_2)
#            z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,index,:])
#            z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,index,:])
            z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,tel,:])
            z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,tel,:])
            dislin.rlmess("%.3f"%(z1)+" $\pm$ "+"%.3f"%(z2),les_posx[tel],yoffset_2)
		#
        dislin.reset("fixspc")
        dislin.sendbf()
        dislin.endgrf()
		#
        dislin.color("FORE")
        dislin.height(30)
        page = 2
        print ( " END OF PAGE ", page )
        dislin.messag(filename,nxoffset_1,nyoffset)
        dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
        #
        dislin.sendbf()
        dislin.endgrf()
        dislin.newpag()

#    log.close()
    #
#    dislin.sendbf()
#    dislin.endgrf()
#    #
    dislin.reset("fixspc")
    dislin.height(24)
    #
    dislin.reset("nograf")
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
##    #==============================================================================
##    # Time averaged spectrum of VIS2                                      PAGE 3
##    #==============================================================================
    #
    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
    dislin.titlin("",1)
    dislin.titlin("Squared visibility vs wavelength (in microns) ==> VIS2DATA",2)
    dislin.titlin("",3)
    dislin.titlin("",4)
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
    #dislin.height(24)
    dislin.title()
    dislin.reset("nograf")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.height(16)
    posx =  120
    posy =  580
    lenx = 1850
    leny =  330
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    yminval = -0.1
    ymaxval = +1.2
    tickx = tickx_wavelength
    ticky = 0.1
    posy0 = posy
    posy1 = 440    # 450
    #dislin.grid(1,1)
    #
    for base in range(nb_base) :
        posy = posy0 + posy1 * base
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if base == 0 :
            dislin.setgrf("name","name","labels","labels")    # restore default
        else :
            dislin.setgrf("name","name","ticks","labels")    # restore default
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
        #
        for expo in range(nb_expo) :
            index = base + expo * nb_base
            xx = [ float(elem) for elem in oifits.wave_sc ]
            yy = [ float(elem) for elem in oifits.oi_vis2_sc_vis2data[index,:] ]
            err = [ float(elem) for elem in oifits.oi_vis2_sc_vis2err[index,:] ]
            n_color = a_color * expo + b_color
            dislin.setclr(n_color)
            my_curve( xx , yy , len(xx) , band[1] )
            dislin.color("FORE")
            shift = fact_shift_err * expo
            xx = xx[shift::stride_err]
            yy = yy[shift::stride_err]
            err = err[shift::stride_err]
            dislin.hsymbl(1)
            dislin.marker(3)
            dislin.errbar( xx , yy , err , err , len(xx) )
#            xshift = float(expo) * 0.005 * (xmax-xmin)
#            xx = xx[::stride_err]
#            xx_shift = [ elem + xshift for elem in xx ]
#            yy = yy[::stride_err]
#            err = err[::stride_err]
#            dislin.errbar( xx_shift , yy , err , err , len(xx) )
            dislin.hsymbl(35)
            dislin.color("FORE")
        #
        xpos = xmin + 0.05 * ( xmax - xmin )
        ypos = ymaxval - 0.15 * ( ymaxval - yminval )
        dislin.height(36)
        dislin.rlmess(base_list[base],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page += 1
    print ( " END OF PAGE ", page )
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    dislin.newpag()
    #
##    #==============================================================================
##    # Time averaged spectrum of VISAMP                                    PAGE 4
##    #==============================================================================
    #
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.height(16)
    dislin.ticks(1,"XY")
    #
    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
    dislin.titlin("",1)
    dislin.titlin("Time averaged visibility amp. vs wavelength (in microns) ==> VISAMP",2)
    dislin.titlin("",3)
    dislin.titlin("",4)
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )
    dislin.title()
    dislin.reset("nograf")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.height(16)
    posx =  120
    posy =  580
    lenx = 1850
    leny =  330
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    yminval = -0.1
    ymaxval = +1.6
    tickx = tickx_wavelength
    ticky = 0.2
    posy0 = posy
    posy1 = 440     # 450
##    dislin.grid(1,1)
    #
    for base in range(nb_base) :
        posy = posy0 + posy1 * base
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if base == 0 :
            dislin.setgrf("name","name","labels","labels")    # restore default
        else :
            dislin.setgrf("name","name","ticks","labels")    # restore default
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,0.0,ticky)
        #
        for expo in range(nb_expo) : 
            index = base + expo * nb_base
            xx = [ float(elem) for elem in oifits.wave_sc ]
            yy = [ float(elem) for elem in oifits.oi_vis_sc_visamp[index,:] ]
            err = [ float(elem) for elem in oifits.oi_vis_sc_visamperr[index,:] ]
            n_color = a_color * expo + b_color
            dislin.setclr(n_color)
            my_curve( xx , yy , len(xx) , band[1] )
            dislin.color("FORE")
            shift = fact_shift_err * expo
            xx = xx[shift::stride_err]
            yy = yy[shift::stride_err]
            err = err[shift::stride_err]
            dislin.hsymbl(1)
            dislin.marker(3)
            dislin.errbar( xx , yy , err , err , len(xx) )
            dislin.hsymbl(35)
            dislin.color("FORE")
        #
        xpos = xmin + 0.05 * ( xmax - xmin )
        ypos = ymaxval - 0.15 * ( ymaxval - yminval )
        dislin.height(36)
        dislin.rlmess(base_list[base],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page += 1
    print ( " END OF PAGE ", page )
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    dislin.newpag()
    #
##    #==============================================================================
##    # T3PHI spectrum                                                  PAGE 5
##    #==============================================================================
    #
    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
    dislin.titlin("",1)   # to erase previous call
    dislin.titlin("Closure phase (in degres) vs wavelength (in microns) ==> T3PHI",2)
    dislin.titlin("",3)   # to erase previous call
    dislin.titlin("",4)   # to erase previous call
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )
    dislin.title()
    dislin.reset("nograf")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.height(16)
    posx = 120
    posy = 800
    lenx = 1850
    leny =  550
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    ymaxval = numpy.nanmax(numpy.abs(oifits.t3phi_sc))
    yminval = -ymaxval
    tickx = tickx_wavelength
    if ymaxval <= 5.0 :
        ticky = 0.5
    elif ymaxval <= 10.0 :
        ticky = 1.0
    elif ymaxval <= 20.0 :
        ticky = 2.0
    elif ymaxval <= 50.0 :
        ticky = 5.0
    elif ymaxval <= 100.0 :
        ticky = 10.0
    else :
        ticky = 20.0
    ymaxval = ticky * numpy.ceil(ymaxval/ticky)
    yminval = -ymaxval
    posy0 = posy
    posy1 = 650
    #
    for triplet in range(nb_triplet) :
        posy = posy0 + posy1 * triplet
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if triplet == 0 :
            dislin.setgrf("name","name","labels","labels")
        else :
            dislin.setgrf("name","name","ticks","labels")
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
        for expo in range(nb_expo) : 
            index = triplet + expo * nb_triplet
            xx = [ float(elem) for elem in oifits.wave_sc ]
            yy = [ float(elem) for elem in oifits.t3phi_sc[index,:] ]
            err = [ float(elem) for elem in oifits.t3phierr_sc[index,:] ]
            n_color = a_color * expo + b_color
            dislin.setclr(n_color)
            my_curve( xx , yy , len(xx) , band[1] )
            dislin.color("FORE")
            shift = fact_shift_err * expo
            xx = xx[shift::stride_err]
            yy = yy[shift::stride_err]
            err = err[shift::stride_err]
            dislin.hsymbl(1)
            dislin.marker(3)
            dislin.errbar( xx , yy , err , err , len(xx) )
            dislin.hsymbl(35)
            dislin.color("FORE")
        #
        xpos = xmin + 0.1 * ( xmax - xmin )
        ypos = ymaxval - 0.1 * ( ymaxval - yminval )
        dislin.height(36)
        dislin.rlmess(triplet_list[triplet],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page += 1
    print ( " END OF PAGE ", page )
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    dislin.newpag()
    #
##    #==============================================================================
##    # VISPHI spectrum                                                     PAGE 6
##    #==============================================================================
    #
    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
    dislin.titlin("",1)
    dislin.titlin("Differential closure phase (in degres) vs wavelength (in microns) ==> VISPHI",2)
    dislin.titlin("",3)
    dislin.titlin("",4)
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )
    dislin.title()
    dislin.reset("nograf")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.height(16)
    posx0 = 125
    posx1 = 1000
    posy0 = 1050
    posy1 =  875
    lenx = 850
    leny = 775
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    ymaxval = numpy.nanmax(numpy.abs(oifits.oi_vis_sc_visphi))
    yminval = -ymaxval
    tickx = 2.0 * tickx_wavelength
    if ymaxval <= 30.0 :
        ticky = 5.0
    elif ymaxval <= 100.0 :
        ticky = 10.0
    else :
        ticky = 20.0
    ymaxval = ticky * numpy.ceil(ymaxval/ticky)
    yminval = - ymaxval
    #
    for base in range(nb_base) : 
        posx = posx0 + posx1 * (base/3)
        posy = posy0 + posy1 * (base%3)
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if base == 0 :
            dislin.setgrf("name","name","labels","ticks")
        elif base in (1,2) :
            dislin.setgrf("name","name","ticks","ticks")
        elif base == 3 :
            dislin.setgrf("name","name","labels","labels")
        elif base in (4,5) :
            dislin.setgrf("name","name","ticks","labels")
        else : 
            print("PAS PREVU")
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
        #
        for expo in range(nb_expo) :
            index = base + expo * nb_base
            xx = [ float(elem) for elem in oifits.wave_sc ]
            yy = [ float(elem) for elem in oifits.oi_vis_sc_visphi[index,:] ]
            err = [ float(elem) for elem in oifits.oi_vis_sc_visphierr[triplet,:] ]
            n_color = a_color * expo + b_color
            dislin.setclr(n_color)
            my_curve( xx , yy , len(xx) , band[1] )
            dislin.color("FORE")
            shift = fact_shift_err * expo
            xx = xx[shift::stride_err]
            yy = yy[shift::stride_err]
            err = err[shift::stride_err]
            dislin.hsymbl(1)
            dislin.marker(3)
            dislin.errbar( xx , yy , err , err , len(xx) )
            dislin.hsymbl(35)
            dislin.color("FORE")
        #
        xpos = xmin + 0.05 * ( xmax - xmin )
        ypos = ymaxval - 0.10 * ( ymaxval - yminval )
        dislin.height(36)
        dislin.rlmess(base_list[base],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page += 1
    print ( " END OF PAGE ", page )
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    dislin.newpag()
    #
##    #==============================================================================
##    # Flux signals vs wavelength + spatial frequencies                       PAGE 7
##    #==============================================================================
    #
    pos_x_1 =  200
    pos_x_2 = 1750
    #
    pos_y_1 = 1000
    pos_y_2 = 1900
    pos_y_3 = 2800
    #
    len_x_1 = 1400
    len_x_2 = 250
    #
    len_y_1 = 800
    len_y_2 = 650
    #
    dislin.height(16)
    posx = pos_x_1
    posy = pos_y_1
    lenx = len_x_1
    leny = len_y_1
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    avgflux = numpy.zeros((4,oifits.nwave_sc),dtype='d')
    for tel in range(nb_telescope):
        avgflux[tel,:] = numpy.mean(oifits.oi_flux_sc[:,tel,:],axis=0) # average flux over time
    y_scale = 1.0 / scale_flux
    ymaxval = y_scale * numpy.nanmax(avgflux)    # NANMAX
    yminval = 0.0
    tickx = tickx_wavelength
    if ymaxval <= 1.0 :
        ticky = 0.05
    elif ymaxval <= 5.0 :
        ticky = 0.2
    elif ymaxval <= 10.0 :
        ticky = 0.5
    else :
        ticky = 1.0
    ymaxval = ticky * numpy.ceil(ymaxval/ticky)
    #
    dislin.labdig(-2,'Y')
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.setgrf("name","name","labels","labels")
    dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
    #
    dislin.titlin("",1)
    dislin.titlin("",2)
    dislin.titlin("Average spectrum (in "+"%.1e"%(scale_flux)+" photo-e/DIT) vs wavelength (in microns)",3)
    dislin.titlin("",4)
    dislin.height(36)
    dislin.title()
    dislin.height(16)
    #
    colors = [ "RED" , "ORANGE" , "BLUE" , "GREEN" ]
    #
    xx = [ float(elem) for elem in oifits.wave_sc ]
    #
    for tel in range(nb_telescope) :
        yy = [ y_scale * float(elem) for elem in avgflux[tel,:] ]
        dislin.color( colors[tel] )
        my_curve( xx , yy , len(xx) , band[1] )
    #
    dislin.color("FORE")
    dislin.sendbf()
    dislin.endgrf()
    dislin.labdig(1,'Y')
    #
    posx = pos_x_2
    posy = pos_y_1
    lenx = len_x_2
    leny = len_y_1
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.setgrf("none","none","none","none")
    dislin.axends("none","XY")                      # restore defaults
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
    #
    dislin.height(30)
    xoffset_1 = 0.2
    xoffset_2 = 0.7
    yoffset = 0.95
    dislin.rlmess("telescope",xoffset_1,yoffset)
    yoffset = yoffset - 0.10
    for tel in range(nb_telescope) :
        dislin.color("FORE")
        dislin.rlmess("tel_"+str(tel+1),xoffset_1,yoffset)
        dislin.color(colors[tel])
        nsymb = 21
        dislin.rlsymb(nsymb,xoffset_2,yoffset-0.10)
        yoffset = yoffset - 0.25
    #
    dislin.color("FORE")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.height(30)
    page += 1
    print ( " END OF PAGE ", page )
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    lamp = band[2]
    if lamp:
        dislin.disfin()
        if device == "pdf" :
           reportname = filename+".pdf"
           print("Create the PDF report : "+reportname)
        else : 
            print(" DONE : "+filename) 
        print ("==========================================================================")
        return
    #
    # Squarred visibility vs spatial frequencies
    # 
    dislin.height(16)
    posx = pos_x_1
    posy = pos_y_2
    lenx = len_x_1
    leny = len_y_2
    #
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.titlin("",1)
    dislin.titlin("",2)
    dislin.titlin("",3)
    dislin.titlin(" Squarred visibility vs spatial frequencies ",4)
    dislin.setgrf("name","name","labels","labels")
    dislin.color("FORE")
    dislin.axends("none","XY")                      # restore defaults
    #
    xmin = 0.0
    xmax = xmax_spatial_freq
    tickx = tickx_spatial_freq
    ymin = -0.1
    ymax = +1.2
    ticky = 0.1
    dislin.graf( xmin,xmax,xmin,tickx, ymin,ymax,ymin,ticky)
    dislin.height(36)
    dislin.title()
    dislin.height(16)
    #
    for base in range(nb_base) :
        spfreq = oifits.spfreq_sc[base]
        vis2data = oifits.oi_vis2_sc_vis2data[base,:]
#        vis2err = oifits.oi_vis2_sc_vis2err[base,:]
        xx = [ float(elem) for elem in spfreq ]
        yy = [ float(elem) for elem in vis2data ]
#        zz = [ float(elem) for elem in vis2err ]
        n_color = 39*base+30
        dislin.setclr(n_color)
        dislin.curve( xx , yy , len(xx) )
    #
    dislin.color("FORE")
    dislin.height(16)
    dislin.sendbf()
    dislin.endgrf()
    #
    posx = pos_x_2
    posy = pos_y_2
    lenx = len_x_2
    leny = len_y_2
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
    #
    dislin.height(25)
    xoffset_1 = 0.1
    xoffset_2 = 0.7
    yoffset = 0.95
    dislin.rlmess("base color",xoffset_1,yoffset)
    dislin.height(30)
    xoffset_1 = 0.2
    yoffset = yoffset - 0.15
    for base in range(nb_base) :
        dislin.color("FORE")
        dislin.rlmess(base_list[base],xoffset_1,yoffset)
        n_color = 39*base+30
        dislin.setclr(n_color)
        nsymb = 21
        dislin.rlsymb(nsymb,xoffset_2,yoffset)
        yoffset = yoffset - 0.15
    #
    dislin.color("FORE")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.color("FORE")
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    # Phase closure vs maximal spatial frequencies
    #
    dislin.height(16)
    posx =  pos_x_1
    posy = pos_y_3
    lenx = len_x_1
    leny = len_y_2
    #
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.titlin("",1)
    dislin.titlin("",2)
    dislin.titlin("",3)
    dislin.titlin(" Phase closure vs maximal spatial frequencies ",4)
    dislin.setgrf("name","name","labels","labels")
    #
    ymax = +180.0
    ymin = -ymax
    ticky = 20.0
    dislin.graf( xmin,xmax,xmin,tickx, ymin,ymax,ymin,ticky)
    dislin.height(36)
    dislin.title()
    dislin.height(16)
    #
    for triplet in range(nb_triplet) :
        t3_spfreq = oifits.t3_spfreq_sc[triplet]
        t3phi = oifits.t3phi_sc[triplet]
#        t3phierr = oifits.t3phierr_sc[triplet]
        xx = [ float(elem) for elem in t3_spfreq ]
        yy = [ float(elem) for elem in t3phi ]
#        zz = [ float(elem) for elem in t3phierr ]
        n_color = 65*triplet+30
        dislin.setclr(n_color)
        dislin.curve( xx , yy , len(xx) )
    #
    dislin.color("FORE")
    dislin.sendbf()
    dislin.endgrf()
    #
    posx = pos_x_2
    posy = pos_y_3
    lenx = len_x_2
    leny = len_y_2
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
    #
    dislin.height(24)
    xoffset_1 = 0.1
    xoffset_2 = 0.7
    yoffset = 0.95
    dislin.rlmess("triplet color",xoffset_1,yoffset)
    dislin.height(30)
    xoffset_1 = 0.15
    yoffset = yoffset - 0.10
    for triplet in range(nb_triplet) :
        dislin.color("FORE")
        dislin.rlmess(triplet_list[triplet],xoffset_1,yoffset)
        n_color = 65*triplet+30
        dislin.setclr(n_color)
        nsymb = 21
        dislin.rlsymb(nsymb,xoffset_2,yoffset-0.10)
        yoffset = yoffset - 0.25
    #
    dislin.color("FORE")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.disfin()
    #
    if ( device == "pdf" ) :
       reportname = filename+".pdf"
       print("Create the PDF report : "+reportname)
    else : 
       print(" DONE : "+filename) 
    print ("==========================================================================")
#
#==============================================================================
#==============================================================================
#
def my_curve ( xx , yy , dim , waves_all ):
	#
    nb_waves = len( sum(waves_all,() ) )
    nb_bands = nb_waves / 2
    #
    for m in range (nb_bands) :
        waves = waves_all[m]
        minwave = waves[0]
        maxwave = waves[1]
        xx_new = []
        yy_new = []
        for n in range(dim) :
            x = xx[n]
            if ( ( x >= minwave ) and ( x <= maxwave ) ) :
                xx_new.append(x)
                yy_new.append(yy[n])
        dim_new = len(xx_new)
        dislin.curve(xx_new,yy_new,dim_new)
#
#==============================================================================
#==============================================================================
#
def my_mean ( yy , ee , limit_err = +sys.float_info.max , 
             limit_min = -sys.float_info.max , limit_max = +sys.float_info.max , 
             fact_e = 4.0 , fact_tol = 1000.0 , stride = 10 ):
    #
    # fact_tol provisoirement pour pallier a l'absence de certaines erreurs prises encore a zero 
    #                de meme pour stride 
    # faire les securites
    #
    ee_max = max(ee)
    if ( ( ee_max <= 0.0 ) & ( fact_tol <= 1.0 ) ) :
        print( " PROBLEM for the value of fact_tol" )
    #
    dim = len(yy)
    yy_prov = yy[::stride]
    m_prov_0 = numpy.mean(yy_prov)
    yy_prov = yy[stride/2::stride]
    m_prov_1 = numpy.mean(yy_prov)
    m = 0.5 * ( m_prov_0 + m_prov_1 )
    e = abs( m_prov_0 - m_prov_1 )
    e_tol = fact_tol * e
    #
    # recherche du premier "bon" point    
    #
    m_prov = m
    n_start = 0
    for n in range(dim) :
        y = yy[n]
        e = ee[n]
        tol = max(fact_e*e,e_tol)
        if ( ( y < limit_min ) | ( y > limit_max ) ) :
            continue
        if ( abs(y-m_prov) <= tol ) : break
        break
    #
    n_start = n
    #
    # calcul de la moyenne robuste
    #
    n_ok = 0
    n_min = 0
    n_max = 0
    n_err = 0
    n_tol = 0
    mean = 0.0
    err = 0.0
    for n in range(n_start,dim) :
        y = yy[n]
        e = ee[n]
        tol = max(fact_e*e,e_tol)
        if ( y < limit_min )  :
            n_min += 1
            continue
        if ( y > limit_max ) :
            n_max += 1
            continue
        if ( e > limit_err ) :
            n_err += 1
            continue
        if ( abs(y-mean) <= tol ) :
            n_ok += 1
            mean += (y-mean) / float(n_ok)
            err  += (e-err)  / float(n_ok)
        else :
            n_tol += 1
    #
    n_limits = (n_ok,n_min,n_max,n_err,n_tol,n_start)
    return(mean,err,n_limits)
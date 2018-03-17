#! /usr/bin/env python
# -*- coding: iso-8859-15 -*-
"""
Created on Mon May 25 13:51:41 2015
@author: bpichon
"""
import datetime
import numpy
import dislin

#==============================================================================
# Preparation of the report using reportlab
#==============================================================================
def produce_oifits_report_common(oifits,filename,band):

    base_list=["12", "13", "14", "23", "24", "34"]
    tel_list = numpy.array([[0,1],[0,2],[0,3],[1,2],[1,3],[2,3]])
    nbase = len(base_list)            # avant code en dir ==> Kervella
    ntel = len(tel_list)              # avant code en dir ==> Kervella
#    blank_10 = "          "
#    blank_20 = "                    "
#    blank_30 = "                              "
#    blank_40 = "                                        "

    device = "pdf"
#    device = "xwin"
#    device = "cons"
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
        dislin.errmod("protocol","file")

#    dislin.pagera()
    dislin.hwfont()
    dislin.texmod("ON")
    dislin.namdis(20,"XY")    # default 30
    dislin.labdig(-2,'X')

    nxoffset_1 = 80
    nxoffset_2 = 1900
    nyoffset = 2920

    #==============================================================================
    # Global parameters
    #==============================================================================
    
#BP    sobjname =  get_key_withdefault(oifits.header,'HIERARCH ESO INS SOBJ NAME','')
#BP    objname  =  get_key_withdefault(oifits.header,'HIERARCH ESO FT ROBJ NAME','')
    
    triplet_list = [str(oifits.t3_sc_staindex[0,:]),str(oifits.t3_sc_staindex[1,:]),\
                    str(oifits.t3_sc_staindex[2,:]),str(oifits.t3_sc_staindex[3,:])]
    telescope_list = [ "Tel_1" , "Tel_2" , "Tel_3" , "Tel_4" ]

    minwave = numpy.nanmin(oifits.wave_sc) # minimum wavelength for FT plots
    maxwave = numpy.nanmax(oifits.wave_sc) # maximum wavelength for FT plots
    minwaveSC = band[1]                 # different from FT if specified by user
    maxwaveSC = band[2]

    print ("Wavelength limits for the data (microns) from %.4f to %.4f "%(minwave,maxwave))        # BP
    print ("Wavelength limits for    plots (microns) from %.4f to %.4f "%(minwaveSC,maxwaveSC))    # BP

    #==============================================================================
    # Selection of wavelength range for plots
    #==============================================================================

    waveok = numpy.squeeze(numpy.where((oifits.wave_sc>=minwaveSC) * (oifits.wave_sc<=maxwaveSC)))

    oifits.wave_sc = oifits.wave_sc[waveok]
    oifits.nwave_sc = len(waveok)
    oifits.spfreq_sc = oifits.spfreq_sc[:,waveok]
    oifits.t3_spfreq_sc = oifits.t3_spfreq_sc[:,waveok]
    
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
#    #==============================================================================
#    # Resampling factor to avoid too large PDF files
#    #==============================================================================
#    

#   resamp = 3 if oifits.nwave_sc > 500 else 1

    scale_flux = 1.0e+05
    stride_err = 10
#
#    #==============================================================================
#    # Start Story
#    #==============================================================================

    if hasattr(oifits,'stations'):
        config =   'GV1='+oifits.telnames[0]+'='+oifits.stations[0]+ \
                 ', GV2='+oifits.telnames[1]+'='+oifits.stations[1]+ \
                 ', GV3='+oifits.telnames[2]+'='+oifits.stations[2]+ \
                 ', GV4='+oifits.telnames[3]+'='+oifits.stations[3]
    else:
        config = ''
#BP
#BP qq lignes
#BP
    ambi = ' not yet defined '
    
    objname  = "Artifical Source"    #BP FAUX BIEN SUR
    objalpha = "25.99"          #BP FAUX BIEN SUR
    objdelta = "99.99"          #BP FAUX BIEN SUR
    objmagL  = "-7.0"           #BP FAUX BIEN SUR
    objmagM  = "-8.0"           #BP FAUX BIEN SUR
    objmagN  = "-9.0"           #BP FAUX BIEN SUR

    xoffset_1 = 0.02
    xoffset_2 = 0.45 
    xoffset_3 = xoffset_1 + xoffset_2
#
# 100 pts par cm 
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
    dislin.height(26)
    
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
    
    yoffset = 0.45
    dislin.rlmess("NDIT x DIT",xoffset_1,yoffset)
    dislin.rlmess( str(oifits.header['HIERARCH ESO DET NDIT'])+' x '+str(oifits.header['HIERARCH ESO DET SEQ1 DIT'])+' s',xoffset_3,yoffset)

    yoffset = 0.35
    dislin.rlmess("Object name",xoffset_1,yoffset)
    dislin.rlmess(objname,xoffset_3,yoffset)

    yoffset = 0.25
    if band[0] == "LM":
        dislin.rlmess("Object RA, Dec, L, M",xoffset_1,yoffset)
        dislin.rlmess(objalpha +'  '+ objdelta +'  L = '+ objmagL +'  M = '+ objmagM,xoffset_3,yoffset)
    elif band[0] == "N":
        dislin.rlmess("Object RA, Dec, N",xoffset_1,yoffset)
        dislin.rlmess(objalpha +'  '+ objdelta +'  N = '+ objmagN,xoffset_3,yoffset)
    else:
        print ("IMPOSSIBLE !!!")

    yoffset = 0.15
    dislin.rlmess("Telescope stations",xoffset_1,yoffset)
    dislin.rlmess(config,xoffset_3,yoffset)

    yoffset = 0.05
    dislin.rlmess("Seeing Wind T0(V) T0(K)",xoffset_1,yoffset)
    dislin.rlmess(ambi,xoffset_3,yoffset)

    dislin.sendbf()
    dislin.endgrf()

#    #==============================================================================
#    # Tables of numerical values
#    #==============================================================================
#
    #
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")

    dislin.axspos(100,1500)
    dislin.axslen(1900,650)
#    dislin.nograf()

    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )

    xoffset_1 = 0.02
#    xoffset_2 = 0.40 
#    xoffset_3 = xoffset_1 + xoffset_2

    dislin.height(30)
    yoffset = 0.95
    dislin.rlmess("Col 1 : Baseline",xoffset_1,yoffset)
    yoffset = 0.88
    dislin.rlmess("Col 2 : Average squared visibility per baseline (vis^2 $\pm$ std) ==> page 2",xoffset_1,yoffset)
    yoffset = 0.81
    dislin.rlmess("Col 3 : Average differential visibility per baseline (vis $\pm$ std) ==> page 3",xoffset_1,yoffset)
    yoffset = 0.74
    dislin.rlmess("Col 4 : Average diffential phase per baseline (visphi $\pm$ std), in degrees ==> page 5",xoffset_1,yoffset)
#
    ixwave=[]
    ixwaveedge=[]
    ixwaveM=[]
    flagM=False
    if (oifits.header['HIERARCH ESO DET CHIP NAME'] == "AQUARIUS"):
        for i in range(len(oifits.wave_sc)):
            if (oifits.wave_sc[i]>=9.5 and oifits.wave_sc[i]<=10.):
                ixwave.append(i)
            if (oifits.wave_sc[i]>=8.25 and oifits.wave_sc[i]<=8.75):
                ixwaveedge.append(i)
    if (oifits.header['HIERARCH ESO DET CHIP NAME'] == "HAWAII-2RG"):
        for i in range(len(oifits.wave_sc)):
            if (oifits.wave_sc[i]>=3.3 and oifits.wave_sc[i]<=3.5):
                ixwave.append(i)
            if ((oifits.wave_sc[i]>=3.0 and oifits.wave_sc[i]<=3.1) or (oifits.wave_sc[i]>=4.0 and oifits.wave_sc[i]<=4.1)):
                ixwaveedge.append(i)
            if (oifits.wave_sc[i]>=4.6 and oifits.wave_sc[i]<=4.7):
                ixwaveM.append(i)
                flagM=True
    #
    dislin.height(36)
    xoffset_1 = 0.02
    xoffset_2 = 0.12
    xoffset_3 = 0.65
    xoffset_4 = 0.80
    yoffset = 0.64
    dislin.rlmess("Base",xoffset_1,yoffset)
    if (flagM == False):
        dislin.rlmess("      vis^2(central/edge)",xoffset_2,yoffset)
    else:
        dislin.rlmess("      vis^2(central/edge/M)",xoffset_2,yoffset)
    dislin.rlmess("      vis_diff",xoffset_3,yoffset)
    dislin.rlmess("      phi_diff",xoffset_4,yoffset)

#    les_posy = [0.60 , 0.50 , 0.40 , 0.30 , 0.20 , 0.10 ]
    les_posy = [0.54 , 0.44 , 0.34 , 0.24 , 0.14 , 0.04 ]
    #

    vis2log=numpy.zeros((6,3),dtype=numpy.float)
    visdifflog=numpy.zeros(6,dtype=numpy.float)
    phidifflog=numpy.zeros(6,dtype=numpy.float)
    philog=numpy.zeros((4,3),dtype=numpy.float)
    
    dislin.fixspc(0.50)
    for base in range(6) :
        dislin.rlmess(base_list[base],xoffset_1,les_posy[base])
#        z1 = numpy.mean(oifits.oi_vis2_sc_vis2data[base,:])
        z1 = numpy.average(oifits.oi_vis2_sc_vis2data[base,ixwave],weights=oifits.oi_vis2_sc_vis2err[base,ixwave]**(-2))
        z2 = numpy.median(oifits.oi_vis2_sc_vis2err[base,ixwave])
        z3 = numpy.average(oifits.oi_vis2_sc_vis2data[base,ixwaveedge],weights=oifits.oi_vis2_sc_vis2err[base,ixwaveedge]**(-2))
        z4 = numpy.median(oifits.oi_vis2_sc_vis2err[base,ixwaveedge])
        vis2log[base,0]=z1
        vis2log[base,1]=z3
        if (flagM == False):
            dislin.rlmess("%.3f"%(z1)+" $\pm$"+"%.2f"%(z2)+" / %.3f"%(z3)+" $\pm$"+"%.2f"%(z4),xoffset_2,les_posy[base])
        else:
            z5 = numpy.average(oifits.oi_vis2_sc_vis2data[base,ixwaveM],weights=oifits.oi_vis2_sc_vis2err[base,ixwaveM]**(-2))
            z6 = numpy.median(oifits.oi_vis2_sc_vis2err[base,ixwaveM])
            vis2log[base,2]=z5
            dislin.rlmess("%.3f"%(z1)+" $\pm$"+"%.2f"%(z2)+" / %.3f"%(z3)+" $\pm$"+"%.2f"%(z4)+" / %.3f"%(z5)+" $\pm$"+"%.2f"%(z6),xoffset_2,les_posy[base])
            
#    ===> Code modifie car les erreurs sont toutes nulles dans les donnees  <=== 
        z1 = numpy.average(oifits.oi_vis_sc_visamp[base,ixwaveedge])
        z2 = numpy.median(oifits.oi_vis_sc_visamperr[base,ixwaveedge])
        visdifflog[base]=z1
        dislin.rlmess("%+.3f"%(z1)+" $\pm$"+"%.2f"%(z2),xoffset_3,les_posy[base])
        if (numpy.max(oifits.oi_vis_sc_visphierr[base,ixwave]) == 0):
            z1 = numpy.average(oifits.oi_vis_sc_visphi[base,ixwave])
        else:
            z1 = numpy.average(oifits.oi_vis_sc_visphi[base,ixwaveedge],weights=oifits.oi_vis_sc_visphierr[base,ixwaveedge] ) * (180.0/numpy.pi)
        z2 = numpy.median(oifits.oi_vis_sc_visphierr[base,ixwaveedge]) * (180.0/numpy.pi)
        phidifflog[base]=z1
        dislin.rlmess("%+8.3f"%(z1)+" $\pm$"+"%.2f"%(z2),xoffset_4,les_posy[base])
    #
    dislin.reset("fixspc")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.height(24)
    #
    dislin.reset("nograf")
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.axspos(100,1825)
    dislin.axslen(1900,250)
#    dislin.nograf()
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )
    #
    dislin.height(30)
    yoffset = 0.92
    dislin.rlmess("Average closure phase per triplet (t3phi $\pm$ std), in degrees ==> page 4",xoffset_1,yoffset)

    dislin.height(36)
    xoffset = 0.02
    yoffset_1 = 0.72
    yoffset_2 = 0.52
    dislin.rlmess("Triplet",xoffset,yoffset_1)
    dislin.rlmess("Phi-center",xoffset,yoffset_2)
    dislin.rlmess("Phi-edge",xoffset,yoffset_2-0.2)
    if (flagM == True):
        dislin.rlmess("Phi-M",xoffset,yoffset_2-0.4)
    les_posx = [0.19 , 0.38 , 0.59 , 0.81 ]
    dislin.fixspc(0.5)

    for triplet in range(4) :
        dislin.rlmess(triplet_list[triplet],les_posx[triplet],yoffset_1)
        z1 = numpy.average(oifits.t3phi_sc[triplet,ixwave],weights=oifits.t3phierr_sc[triplet,ixwave]) * (180.0/numpy.pi)
        z2 = numpy.median(oifits.t3phierr_sc[triplet,ixwave]) * (180.0/numpy.pi)
        philog[triplet,0]=z1
        dislin.rlmess("%+.2f"%(z1)+" $\pm$ "+"%.2f"%(z2),les_posx[triplet],yoffset_2)
    for triplet in range(4) :
        dislin.rlmess(triplet_list[triplet],les_posx[triplet],yoffset_1)
        z1 = numpy.average(oifits.t3phi_sc[triplet,ixwaveedge],weights=oifits.t3phierr_sc[triplet,ixwaveedge]) * (180.0/numpy.pi)
        philog[triplet,1]=z1
        z2 = numpy.median(oifits.t3phierr_sc[triplet,ixwaveedge]) * (180.0/numpy.pi)
        dislin.rlmess("%+.2f"%(z1)+" $\pm$ "+"%.2f"%(z2),les_posx[triplet],yoffset_2-0.2)
    if (flagM == True):
        for triplet in range(4) :
            dislin.rlmess(triplet_list[triplet],les_posx[triplet],yoffset_1)
            z1 = numpy.average(oifits.t3phi_sc[triplet,ixwaveM],weights=oifits.t3phierr_sc[triplet,ixwaveM]) * (180.0/numpy.pi)
            philog[triplet,2]=z1
            z2 = numpy.median(oifits.t3phierr_sc[triplet,ixwaveM]) * (180.0/numpy.pi)
            dislin.rlmess("%+.2f"%(z1)+" $\pm$ "+"%.2f"%(z2),les_posx[triplet],yoffset_2-0.4)
        
    
    dislin.reset("fixspc")
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.reset("nograf")
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.axspos(100,2150)
    dislin.axslen(1900,250)

    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.graf( 0.0,1.0,0.0,1.0, -0.05,1.0,0.0,1.0 )

    dislin.height(30)
    xoffset = 0.02
    yoffset = 0.92
    dislin.rlmess("Average photometric flux ("+"%.1e"%(scale_flux)+" photo-e-/s/sp.channel $\pm$ std) ==> page 6",xoffset,yoffset)

    dislin.height(36)
    xoffset = 0.02
    yoffset_1 = 0.55
    yoffset_2 = 0.20
    y_scale = 1.0 / scale_flux
    dislin.rlmess("Telescope",xoffset,yoffset_1)
    dislin.rlmess("Flux",xoffset,yoffset_2)
    #
#    les_posx = [0.20 , 0.55 , 0.31 , 0.66 ]
#    for tel in range(0,2) :
#        dislin.rlmess(telescope_list[tel],les_posx[tel],yoffset_1)
#        z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,tel,:])
#        z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,tel,:])
##        dislin.rlmess("%.3e"%(z1)+" $\pm$ "+"%.3e"%(z2),les_posx[tel],yoffset_2)
#        dislin.rlmess("%.3f"%(z1)+" $\pm$ "+"%.3f"%(z2),les_posx[tel],yoffset_2)
#    for tel in range(2,4) :
#        dislin.rlmess(telescope_list[tel],les_posx[tel],yoffset_3)
#        z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,tel,:])
#        z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,tel,:])
##        dislin.rlmess("%.3e"%(z1)+" $\pm$ "+"%.3e"%(z2),les_posx[tel],yoffset_4)
#        dislin.rlmess("%.3f"%(z1)+" $\pm$ "+"%.3f"%(z2),les_posx[tel],yoffset_4)
    les_posx = [0.20 , 0.40 , 0.60 , 0.80 ]
    for tel in range(4) :
        dislin.rlmess(telescope_list[tel],les_posx[tel],yoffset_1)
        z1 = y_scale * numpy.average(oifits.oi_flux_sc[:,tel,:])
        z2 = y_scale * numpy.median(oifits.oi_flux_err_sc[:,tel,:])
        dislin.rlmess("%.3f"%(z1)+" $\pm$ "+"%.3f"%(z2),les_posx[tel],yoffset_2)
    #
    dislin.height(30)
    page = 1
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    dislin.sendbf()
    dislin.endgrf()
    dislin.newpag()
    
##    #==============================================================================
##    # Time averaged spectrum of VIS2
##    #==============================================================================
    #
#    dislin.height(24)
    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
#    dislin.titlin("Squared visibility vs. wavelength (in microns) ==> VIS2DATA",1)
#    dislin.titlin("Squared visibility vs. wavelength (microns)",2)
    dislin.titlin("Squared visibility vs wavelength (in microns) ==> VIS2DATA",2)
    dislin.titlin("",1)   # to erase previous call
#    dislin.titlin("",2)   # to erase previous call
    dislin.titlin("",3)   # to erase previous call
    dislin.titlin("",4)   # to erase previous call
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
#    dislin.height(24)
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
    posx =  120      # 150
    posy =  580      # 600
    lenx = 1850      # 1900
    leny =  330      # 350
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    yminval = -0.1
    ymaxval = +1.2
    tickx = band[3]
    ticky = 0.1
    posy0 = posy
    posy1 = 440    # 450
##    dislin.grid(1,1)
    #
    for baseline in range(6) :
        posy = posy0 + posy1 * baseline
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if baseline == 0 :
            dislin.setgrf("name","name","labels","labels")    # restore default
        else :
            dislin.setgrf("name","name","ticks","labels")    # restore default
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
        xx = [ float(elem) for elem in oifits.wave_sc ]
        yy = [ float(elem) for elem in oifits.oi_vis2_sc_vis2data[baseline,:] ]
        dislin.color("RED")
        dislin.curve( xx , yy , len(xx) )
        dislin.color("BLUE")
        err = [ float(elem) for elem in oifits.oi_vis2_sc_vis2err[baseline,:] ]
        xx = xx[::stride_err]
        yy = yy[::stride_err]
        err = err[::stride_err]
        dislin.hsymbl(1)
        dislin.marker(3)
        dislin.errbar( xx , yy , err , err , len(xx) )
        dislin.hsymbl(35)
        dislin.color("FORE")
        xpos = xmin + 0.05 * ( xmax - xmin )
        ypos = ymaxval - 0.15 * ( ymaxval - yminval )
        dislin.height(36)
#        dislin.rlmess("VIS2DATA "+base_list[baseline],xpos,ypos)
        dislin.rlmess(base_list[baseline],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page = 2
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    dislin.newpag()

##    #==============================================================================
##    # Time averaged spectrum of VISAMP
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
    dislin.titlin("",1)   # to erase previous call
    dislin.titlin("Time averaged visibility amp. vs wavelength (in microns) ==> VISAMP",2)
    dislin.titlin("",3)   # to erase previous call
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
    dislin.setgrf("name","name","ticks","ticks")    # restore default
    dislin.axends("none","XY")                      # restore defaults
    dislin.ticks(1,"XY")
    #
    dislin.height(16)
    posx =  120     # 150
    posy =  580     # 600
    lenx = 1850     # 1900
    leny =  330     # 350

    xmin = minwaveSC
    xmax = maxwaveSC
    yminval = -0.1
    ymaxval = +1.6
    tickx = band[3]
    ticky = 0.2
    posy0 = posy
    posy1 = 440     # 450
##    dislin.grid(1,1)
    #
    for baseline in range(6) :
        posy = posy0 + posy1 * baseline
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if baseline == 0 :
            dislin.setgrf("name","name","labels","labels")    # restore default
        else :
            dislin.setgrf("name","name","ticks","labels")    # restore default
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,0.0,ticky)
        xx = [ float(elem) for elem in oifits.wave_sc ]
        yy = [ float(elem) for elem in oifits.oi_vis_sc_visamp[baseline,:] ]
        dislin.color("RED")
        dislin.curve( xx , yy , len(xx) )
        dislin.color("BLUE")
        err = [ float(elem) for elem in oifits.oi_vis_sc_visamperr[baseline,:] ]
        xx = xx[::stride_err]
        yy = yy[::stride_err]
        err = err[::stride_err]
        dislin.hsymbl(1)
        dislin.marker(3)
        dislin.errbar( xx , yy , err , err , len(xx) )
        dislin.hsymbl(35)
        dislin.color("FORE")
        xpos = xmin + 0.05 * ( xmax - xmin )
        ypos = ymaxval - 0.15 * ( ymaxval - yminval )
        dislin.height(36)
#        dislin.rlmess("VISAMP "+base_list[baseline],xpos,ypos)
        dislin.rlmess(base_list[baseline],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page = 3
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    dislin.newpag()

##    #==============================================================================
##    # T3PHI spectrum 
##    #==============================================================================

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
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
#    dislin.height(24)
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
    ymaxval = numpy.nanmax(numpy.abs(oifits.t3phi_sc))*(180./numpy.pi)      # NANMAX
    yminval = -ymaxval
    tickx = band[3]
    if ymaxval <= 5.0 :
        ticky = 0.5
    elif ymaxval <= 10.0 :
        ticky = 1.0
    else :
        ticky = 2.0
    ymaxval = ticky * numpy.ceil(ymaxval/ticky)
    yminval = - ymaxval
    posy0 = posy
    posy1 = 650     # 450
    #
    for triplet in range(4) :
        posy = posy0 + posy1 * triplet
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if triplet == 0 :
            dislin.setgrf("name","name","labels","labels")
        else :
            dislin.setgrf("name","name","ticks","labels")
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
        xx = [ float(elem) for elem in oifits.wave_sc ]
        yy = [ float(elem) * (180./numpy.pi) for elem in oifits.t3phi_sc[triplet,:] ]
        dislin.color("RED")
        dislin.curve( xx , yy , len(xx) )
        dislin.color("BLUE")
        err = [ float(elem) * (180./numpy.pi) for elem in oifits.t3phierr_sc[triplet,:] ]
        xx = xx[::stride_err]
        yy = yy[::stride_err]
        err = err[::stride_err]
        dislin.hsymbl(1)
        dislin.marker(3)
        dislin.errbar( xx , yy , err , err , len(xx) )
        dislin.hsymbl(35)
        dislin.color("FORE")
        xpos = xmin + 0.1 * ( xmax - xmin )
        ypos = ymaxval - 0.1 * ( ymaxval - yminval )
        dislin.height(36)
        dislin.rlmess(triplet_list[triplet],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page = 4
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    dislin.newpag()

##    #==============================================================================
##    # VISPHI spectrum
##    #==============================================================================

    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
    dislin.titlin("",1)   # to erase previous call
    dislin.titlin("Differential closure phase (in degres) vs wavelength (in microns) ==> VISPHI",2)
    dislin.titlin("",3)   # to erase previous call
    dislin.titlin("",4)   # to erase previous call
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
#    dislin.height(24)
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
    ymaxval = numpy.nanmax(numpy.abs(oifits.oi_vis_sc_visphi))*(180./numpy.pi)    # NANMAX
    yminval = -ymaxval
    tickx = 2.0 * band[3]
    if ymaxval <= 30.0 :
        ticky = 5.0
    elif ymaxval <= 100.0 :
        ticky = 10.0
    else :
        ticky = 20.0
    ymaxval = ticky * numpy.ceil(ymaxval/ticky)
    yminval = - ymaxval
    #
    for baseline in range(nbase) : 
        posx = posx0 + posx1 * (baseline/3)
        posy = posy0 + posy1 * (baseline%3)
        dislin.axspos(posx,posy)
        dislin.axslen(lenx,leny)
        if baseline == 0 :
            dislin.setgrf("name","name","labels","ticks")
        elif baseline in (1,2) :
            dislin.setgrf("name","name","ticks","ticks")
        elif baseline == 3 :
            dislin.setgrf("name","name","labels","labels")
        elif baseline in (4,5) :
            dislin.setgrf("name","name","ticks","labels")
        else : 
            print("PAS PREVU")
        dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)
        xx = [ float(elem) for elem in oifits.wave_sc ]
        yy = [ float(elem) * (180./numpy.pi) for elem in oifits.oi_vis_sc_visphi[baseline,:] ]
        dislin.color("RED")
        dislin.curve( xx , yy , len(xx) )
        dislin.color("BLUE")
        err = [ float(elem) * (180./numpy.pi) for elem in oifits.oi_vis_sc_visphierr[triplet,:] ]
        xx = xx[::stride_err]
        yy = yy[::stride_err]
        err = err[::stride_err]
        dislin.hsymbl(1)
        dislin.marker(3)
        dislin.errbar( xx , yy , err , err , len(xx) )
        dislin.hsymbl(35)
        dislin.color("FORE")
        xpos = xmin + 0.05 * ( xmax - xmin )
        ypos = ymaxval - 0.10 * ( ymaxval - yminval )
        dislin.height(36)
        dislin.rlmess(base_list[baseline],xpos,ypos)
        dislin.height(16)
        dislin.sendbf()
        dislin.endgrf()
    #
    dislin.height(30)
    page = 5
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #
    dislin.newpag()


##    #==============================================================================
##    # Flux signals vs wavelength
##    #==============================================================================
##
    dislin.height(36)
    dislin.axspos(150,600)
    dislin.axslen(1900,350)
    dislin.titlin("",1)   # to erase previous call
#    dislin.titlin("Average spectrum (in "+"%.1e"%(scale_flux)+" photo-e/DIT) vs wavelength (in microns) ==> OI_FLUX",2)
#    dislin.titlin("Tel1 = red, Tel2 = orange, Tel3 = blue, Tel4 = green",3)
    dislin.titlin("Average spectrum (in "+"%.1e"%(scale_flux)+" photo-e/DIT) vs wavelength (in microns)",2)
    dislin.titlin("==> OI_FLUX  ;  Tel1 = red, Tel2 = orange, Tel3 = blue, Tel4 = green",3)
    dislin.titlin("",4)   # to erase previous call
    dislin.setgrf("none","none","none","none")
    dislin.axends("noends","XY")
    dislin.nograf()
    dislin.graf( 0.0,1.0,0.0,1.0, 0.0,1.0,0.0,1.0 )    # dummy for title only
#    dislin.height(24)
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
    posx =  200
    posy = 1300
    lenx = 1700
    leny = 1000
    #
    xmin = minwaveSC
    xmax = maxwaveSC
    avgflux = numpy.zeros((4,oifits.nwave_sc),dtype='d')
    for tel in range(4):
        avgflux[tel,:] = numpy.nanmean(oifits.oi_flux_sc[:,tel,:],axis=0) # average flux over time
    y_scale = 1.0 / scale_flux
    ymaxval = y_scale * numpy.nanmax(avgflux)    # NANMAX
    yminval = 0.0
    tickx = band[3]
    if ymaxval <= 1.0 :
        ticky = 0.05
    elif ymaxval <= 5.0 :
        ticky = 0.2
    else :
        ticky = 0.5
    ymaxval = ticky * numpy.ceil(ymaxval/ticky)
    #
    dislin.axspos(posx,posy)
    dislin.axslen(lenx,leny)
    dislin.setgrf("name","name","labels","labels")
    dislin.graf( xmin,xmax,xmin,tickx, yminval,ymaxval,yminval,ticky)

    xx = [ float(elem) for elem in oifits.wave_sc ]
    tel = 0
    yy = [ y_scale * float(elem) for elem in avgflux[tel,:] ]
    dislin.color("RED")
    dislin.curve( xx , yy , len(xx) )
    tel = 1
    yy = [ y_scale * float(elem) for elem in avgflux[tel,:] ]
    dislin.color("ORANGE")
    dislin.curve( xx , yy , len(xx) )
    tel = 2
    yy = [ y_scale * float(elem) for elem in avgflux[tel,:] ]
    dislin.color("BLUE")
    dislin.curve( xx , yy , len(xx) )
    tel = 3
    yy = [ y_scale * float(elem) for elem in avgflux[tel,:] ]
    dislin.color("GREEN")
    dislin.curve( xx , yy , len(xx) )
    dislin.color("FORE")
    #
    dislin.sendbf()
    dislin.endgrf()
    #
    dislin.height(30)
    page = 6
    dislin.messag(filename,nxoffset_1,nyoffset)
    dislin.messag("Page : "+str(page),nxoffset_2,nyoffset)
    #

##    #==============================================================================
##    # Create PDF report from Story
##    #==============================================================================

    lamp = band[4]
    if lamp:
        dislin.disfin()
        if device == "pdf" :
           reportname = filename+".pdf"
           print("Create the PDF report : "+reportname)
        else : 
            print(" DONE : "+filename) 
        print ("==========================================================================")
    else:        
        dislin.newpag()


    if band[0] == "LM":
        oml1=oifits.header['HIERARCH ESO INS OML1 VOLT']
        oml2=oifits.header['HIERARCH ESO INS OML2 VOLT']
        oml3=oifits.header['HIERARCH ESO INS OML3 VOLT']
        oml4=oifits.header['HIERARCH ESO INS OML4 VOLT']
        print oml1,oml2,oml3,oml4
        if (numpy.abs(oml1-4) > 0.15 or numpy.abs(oml2-4) > 0.15 or numpy.abs(oml3-4) > 0.15 or numpy.abs(oml4-4) > 0.15):
            opdmod="YES"
        else:
            opdmod="NO"
        print filename+".fits",";",oifits.header['ESO TPL START'],";",oifits.header['HIERARCH ESO DET CHIP NAME'],";",oifits.header['HIERARCH ESO INS DIL ID'],";",oifits.header['HIERARCH ESO INS PIL ID'],";",oifits.header['HIERARCH ESO INS POL ID'],";",oifits.header['HIERARCH ESO INS FIL ID'],";",oifits.header['HIERARCH ESO INS SFL ID'],";",oifits.header['HIERARCH ESO INS BCD1 ID'],";",oifits.header['HIERARCH ESO INS BCD2 ID'],";",opdmod,";",oifits.header['HIERARCH ESO DET SEQ1 DIT'],";",vis2log[0,0],";",vis2log[1,0],";",vis2log[2,0],";",vis2log[3,0],";",vis2log[4,0],";",vis2log[5,0],";",vis2log[0,1],";",vis2log[1,1],";",vis2log[2,1],";",vis2log[3,1],";",vis2log[4,1],";",vis2log[5,1],";",vis2log[0,2],";",vis2log[1,2],";",vis2log[2,2],";",vis2log[3,2],";",vis2log[4,2],";",vis2log[5,2],";",philog[0,0],";",philog[1,0],";",philog[2,0],";",philog[3,0],";",philog[0,1],";",philog[1,1],";",philog[2,1],";",philog[3,1],";",philog[0,2],";",philog[1,2],";",philog[2,2],";",philog[3,2],";",visdifflog[0],";",visdifflog[1],";",visdifflog[2],";",visdifflog[3],";",visdifflog[4],";",visdifflog[5],";",phidifflog[0],";",phidifflog[1],";",phidifflog[2],";",phidifflog[3],";",phidifflog[4],";",phidifflog[5]

    else:
        omn1=oifits.header['HIERARCH ESO INS OMN1 VOLT']
        omn2=oifits.header['HIERARCH ESO INS OMN2 VOLT']
        omn3=oifits.header['HIERARCH ESO INS OMN3 VOLT']
        omn4=oifits.header['HIERARCH ESO INS OMN4 VOLT']
        print omn1,omn2,omn3,omn4
        if (numpy.abs(omn1-5) > 0.15 or numpy.abs(omn2-5) > 0.15 or numpy.abs(omn3-5) > 0.15 or numpy.abs(omn4-5) > 0.15):
            opdmod="YES"
        else:
            opdmod="NO"
        print filename+".fits",";",oifits.header['ESO TPL START'],";",oifits.header['HIERARCH ESO DET CHIP NAME'],";",oifits.header['HIERARCH ESO INS DIN ID'],";",oifits.header['HIERARCH ESO INS PIN ID'],";",oifits.header['HIERARCH ESO INS PON ID'],";",oifits.header['HIERARCH ESO INS FIN ID'],";",oifits.header['HIERARCH ESO INS SFN ID'],";",oifits.header['HIERARCH ESO INS BCD1 ID'],";",oifits.header['HIERARCH ESO INS BCD2 ID'],";",opdmod,";",oifits.header['HIERARCH ESO DET SEQ1 DIT'],";",vis2log[0,0],";",vis2log[1,0],";",vis2log[2,0],";",vis2log[3,0],";",vis2log[4,0],";",vis2log[5,0],";",vis2log[0,1],";",vis2log[1,1],";",vis2log[2,1],";",vis2log[3,1],";",vis2log[4,1],";",vis2log[5,1],";",vis2log[0,2],";",vis2log[1,2],";",vis2log[2,2],";",vis2log[3,2],";",vis2log[4,2],";",vis2log[5,2],";",philog[0,0],";",philog[1,0],";",philog[2,0],";",philog[3,0],";",philog[0,1],";",philog[1,1],";",philog[2,1],";",philog[3,1],";",philog[0,2],";",philog[1,2],";",philog[2,2],";",philog[3,2],";",visdifflog[0],";",visdifflog[1],";",visdifflog[2],";",visdifflog[3],";",visdifflog[4],";",visdifflog[5],";",phidifflog[0],";",phidifflog[1],";",phidifflog[2],";",phidifflog[3],";",phidifflog[4],";",phidifflog[5]


 

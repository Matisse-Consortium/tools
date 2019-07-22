#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 12:28:03 2019

@author: ame
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
import sys
import mat_show_oifits as msoi
import os
from matplotlib.backends.backend_pdf import PdfPages

inch=1/2.54



sta_name_list=np.array(['A0','A1','B0','B1','B2','B3','B4','B5','C0','C1','C2','C3','D0','D1','D2','E0','G0','G1','G2','H0','I1','J1','J2','J3','J4','J5','J6','K0','L0','M0','U1','U2','U3','U4'])

sta_pos_list=np.array([[-32.001,-48.013,-14.642,-55.812],[-32.001,-64.021, -9.434,-70.949],[-23.991,-48.019, -7.065,-53.212],[-23.991,-64.011, -1.863,-68.334],[-23.991,-72.011,  0.739,-75.899],
         [-23.991,-80.029,  3.348,-83.481],[-23.991,-88.013,  5.945,-91.030],[-23.991,-96.012,  8.547,-98.594],[-16.002,-48.013,  0.487,-50.607],[-16.002,-64.011,  5.691,-65.735],
         [-16.002,-72.019,  8.296,-73.307],[-16.002,-80.010, 10.896,-80.864],[  0.010,-48.012, 15.628,-45.397],[  0.010,-80.015, 26.039,-75.660],[  0.010,-96.012, 31.243,-90.787],
         [ 16.011,-48.016, 30.760,-40.196],[ 32.017,-48.0172,45.896,-34.990],[ 32.020,-112.010,66.716,-95.501],[ 31.995,-24.003, 38.063,-12.289],[ 64.015,-48.007, 76.150,-24.572],
         [ 72.001,-87.997, 96.711,-59.789],[ 88.016,-71.992,106.648,-39.444],[ 88.016,-96.005,114.460,-62.151],[ 88.016,  7.996,  80.628,36.193],[ 88.016,23.993,  75.424, 51.320],
         [ 88.016,47.987,  67.618, 74.009],[ 88.016,71.990,  59.810, 96.706],[ 96.002,-48.006,106.397,-14.165],[104.021,-47.998,113.977,-11.549],[112.013,-48.000,121.535, -8.951],
         [-16.000,-16.000, -9.925,-20.335],[ 24.000,24.000,  14.887, 30.502],[64.0013,47.9725, 44.915, 66.183],[112.000,8.000,  103.306,43.999]])

def _vltiplot(tels=np.array([]),baselines=np.array([]),symsize=2,color='k',tcolor=['k','k','k','k'],bcolor=['r','g','b','m','c','y'],labels=False,axe=None):

    UTs=sta_pos_list[-4:,:]
    ATs=sta_pos_list[:-4,:]

    if not(axe):
        axe=plt
    #plt.axis([np.min(sta_pos[:,2])-20,np.max(sta_pos[:,2])+20,np.min(sta_pos[:,3])-20,np.max(sta_pos[:,3])+20])
    axe.plot(UTs[:,2],UTs[:,3],'ro',markersize=3*symsize,color=color)
    axe.plot(ATs[:,2],ATs[:,3],'ro',markersize=symsize,color=color)

    if labels:
        for i in range(len(sta_name_list)):
            axe.text(sta_pos_list[i,2]+symsize,sta_pos_list[i,3]+symsize,sta_name_list[i],fontsize=2.5*symsize,color=color,horizontalalignment='center',zorder=1)


    if len(baselines)>0:
        for i in range(np.shape(baselines)[0]):
            tel1=np.where(sta_name_list == baselines[i,0])
            tel1x=sta_pos_list[tel1,2]
            tel1y=sta_pos_list[tel1,3]
            tel2=np.where(sta_name_list == baselines[i,1])
            tel2x=sta_pos_list[tel2,2]
            tel2y=sta_pos_list[tel2,3]
            axe.plot([tel1x[0,0],tel2x[0,0]],[tel1y[0,0],tel2y[0,0]],color=bcolor[i],zorder=2)

    if len(tels)>0:
         for i in range(len(tels)):
            tel=np.where(sta_name_list == tels[i])
            telx=sta_pos_list[tel,2]
            tely=sta_pos_list[tel,3]
            #print(telx)
            #print(tely)
            axe.scatter(telx[0,0],tely[0,0],c=tcolor[i],s=40*symsize,zorder=3)

def mat_showOiData(filename,wlRange=None,showErr=False):

    try:
        dic=msoi.open_oi(filename)
    except:
        print("Unable to open {0}".format(filename))
        return

    u=dic['VIS2']['U']
    v=dic['VIS2']['V']
    v2=dic['VIS2']['VIS2']
    errv2=dic['VIS2']['VIS2ERR']
    phi=dic['VIS']['DPHI']
    errphi=dic['VIS']['DPHIERR']
    cp=dic['T3']['CLOS']
    cperr=dic['T3']['CLOSERR']

    try:
        flux=dic['FLUX']['FLUX']
        corrFlux=False
    except:
        corrFlux=True
    wl=dic['WLEN']*1e6

    nwl=np.size(wl)
    c=15
    dash=[[3*c,6*c],[1e-3,3*c,3*c,3*c],[1e-3,6*c,3*c,1e-3]]

    sta_index= dic['VIS2']['STA_INDEX']
    sta_name=[
             [dic['STA_NAME'][np.where(dic['STA_INDEX'] == sti[0])[0][0]],
              dic['STA_NAME'][np.where(dic['STA_INDEX'] == sti[1])[0][0]]]
             for sti in sta_index]

    fig=plt.figure(figsize=(29.7*inch,21*inch))

    plts={}

    cols=['r','g','b','m','c','y']

    colT3=['brm','gbc','gmy','cry']


    colTel=['yellow','magenta','cyan','grey']
    cols=['darkcyan','red','blue','darkmagenta','limegreen','y']
    colT3=['darkblue','k','darkred','darkgreen']

    pltv2=[]
    for i in range(6):
        if i==0:
            plts['VIS2_{0}'.format(i)]=fig.add_axes([0.075, 0.1*(i+0.5), 0.22,0.1])
        else:
            plts['VIS2_{0}'.format(i)]=fig.add_axes([0.075, 0.1*(i+0.5), 0.22,0.1],sharex=plts['VIS2_0'],sharey=plts['VIS2_0'])
            plts['VIS2_{0}'.format(i)].get_xaxis().set_visible(False)
        pltv2.append(plts['VIS2_{0}'.format(i)])
    msoi.show_oi_vs_wlen(dic,key='VIS2',datatype="VIS2",plot_errorbars=showErr,subplotList=pltv2,colorList=cols)
        #plts['VIS2_{0}'.format(i)].plot(wl,v2[i,:],color=cols[i])
    pltphi=[]
    for i in range(6):
        if i==0:
            plts['PHI_{0}'.format(i)]=fig.add_axes([0.375, 0.1*(i+0.5), 0.22,0.1],sharex=plts['VIS2_0'])
        else:
            plts['PHI_{0}'.format(i)]=fig.add_axes([0.375, 0.1*(i+0.5), 0.22,0.1],sharex=plts['VIS2_0'],sharey=plts['PHI_0'])
            plts['PHI_{0}'.format(i)].get_xaxis().set_visible(False)
        pltphi.append(plts['PHI_{0}'.format(i)])
        #plts['PHI_{0}'.format(i)].plot(wl,phi[i,:],color=cols[i])

    msoi.show_oi_vs_wlen(dic,key='VIS',datatype="DPHI",plot_errorbars=showErr,subplotList=pltphi,colorList=cols)

    pltcp=[]
    for i in range(4):
        if i==0:
            plts['CP_{0}'.format(i)]=fig.add_axes([0.675, 0.1*(i+0.5), 0.22,0.1],sharex=plts['VIS2_0'])
        else:
            plts['CP_{0}'.format(i)]=fig.add_axes([0.675, 0.1*(i+0.5), 0.22,0.1],sharex=plts['VIS2_0'],sharey=plts['CP_0'])
            plts['CP_{0}'.format(i)].get_xaxis().set_visible(False)
        pltcp.append(plts['CP_{0}'.format(i)])
    msoi.show_oi_vs_wlen(dic,key='T3',datatype="CLOS",plot_errorbars=showErr,subplotList=pltcp,colorList=colT3)


    plts['FLUX']=fig.add_axes([0.075, 0.70, 0.22,0.15],sharex=plts['VIS2_0'])
    pltflx=[plts['FLUX']]*4
    if not(corrFlux):
        msoi.show_oi_vs_wlen(dic,key='FLUX',datatype="FLUX",plot_errorbars=showErr,subplotList=pltflx,colorList=colTel)


    plts['UV']=fig.add_axes([0.675, 0.53, 0.22,0.22*1.41])
    plts['UV'].set_aspect('equal', 'box')
    plts['UV'].scatter(u,v,c=cols,s=2)
    plts['UV'].scatter(-u,-v,c=cols,s=2)
    plts['UV'].set_xlabel('U (m)')
    plts['UV'].set_ylabel('V (m)')

    plts['VLTI']=fig.add_axes([0.305,0.65,0.32,0.32])
    plts['VLTI'].set_axis_off()
    plts['VLTI'].set_aspect('equal', 'box')
    _vltiplot(tels=dic['STA_NAME'],baselines=np.array(sta_name[0:6]),symsize=2,color='k',bcolor=cols,tcolor=colTel,labels=False,axe=plts['VLTI'])

    if (wlRange):
        plts['VIS2_0'].set_xlim(wlRange)

    plts['UV'].set_xlim(150,-150)
    plts['UV'].set_ylim(-150,150)


    if not(corrFlux):
        label="Square Visibility"
    else:
        label="Correlated flux"
    fig.text(0.02,0.35,label,rotation=90,horizontalalignment='center',verticalalignment='center')
    fig.text(0.32,0.35,"Differential Phase ($^o$)",rotation=90,horizontalalignment='center',verticalalignment='center')
    fig.text(0.62,0.25,"Closure Phase ($^o$)",rotation=90,horizontalalignment='center',verticalalignment='center')
    fig.text(0.185,0.01,"$\lambda$ ($\mu$m)",horizontalalignment='center',verticalalignment='center')
    fig.text(0.485,0.01,"$\lambda$ ($\mu$m)",horizontalalignment='center',verticalalignment='center')
    fig.text(0.785,0.01,"$\lambda$ ($\mu$m)",horizontalalignment='center',verticalalignment='center')

    filename=filename.split("/")[-1]
    fig.text(0.01,0.97,filename)

    telconf=""
    for el in dic['STA_NAME']:
        telconf+=el
        telconf+="-"
    telconf=telconf[:-1]

    fig.text(0.71,0.95,telconf,fontsize=20)

    seeing="{0:.2f}".format(dic['SEEING'])
    coherence="{0:.2f}".format(1000*dic['TAU0'])
    airmass="{0:.2f}".format(dic['HDR']['HIERARCH ESO ISS AIRM START'])
    fig.text(0.71,0.92,"Seeing={0}\"".format(seeing))
    fig.text(0.71,0.90,"Coherence={0}ms".format(coherence))
    fig.text(0.71,0.88,"Airmass={0}".format(airmass))

    return fig


if __name__ == '__main__':
    filename=None
    wlRange=None
    showErr=False
    pdf=False
    merged=False

    listArg = sys.argv
    print(listArg)
    if len(listArg)>1:
        if  listArg[1][0]!="-":
            filename=listArg[1]
    else :
        listArg.append("--help")

    for elt in listArg:
        if ('--help' in elt):
            print("Usage: mat_showOiData.py filenameOrDirectory [--wlRange=[min,max] in micron] [--showErr=FALSE] [--pdf] [--mergedpdf]")
            sys.exit(0)
        if ('--wlRange' in elt):
            wlRange=elt.split("=")[1]
        if ('--showErr' in elt):
            showErr=elt.split("=")[1]
        if ('--pdf' in elt):
            pdf=True
        if ('--mergedpdf' in elt):
            merged=True




    if os.path.isdir(filename):
        os.chdir(filename)
        dir0=os.path.basename(os.getcwd())
        filename=[fi for fi in os.listdir(filename) if ".fits" in fi]
    else:
        filename=[filename]
    nfiles=len(filename)

    if merged:
        pdfname=dir0+'_plots.pdf'
        pdf = PdfPages(pdfname)
        print("saving to {0}".format(pdfname))


    for filei in filename:
        fig=mat_showOiData(filei,wlRange=wlRange,showErr=showErr)
        if pdf and not(merged):
            pdfname=filei.split(".fits")[0]+".pdf"
            plt.savefig(pdfname)
            print("saving to {0}".format(pdfname))
        if merged:
            pdf.savefig(fig)
            plt.clf()
   	if not(pdf or merged):
            plt.show(block=False)

    if merged:
        pdf.close()





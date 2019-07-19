import sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np



def mat_time_flux_plot(filename, skyfilename="",unchop=False):
    
    
    doSky = (skyfilename!="" or unchop==True)
    
    d=fits.open(filename)

    if doSky:
        dsky=fits.open(skyfilename)

    tartype=d['IMAGING_DATA'].data['TARTYP']
    tt=(tartype=='T')*1+(tartype=='U')*0.5
    targetframes=np.argwhere(tt==1)
    skyframes=np.argwhere(tt==0)
    

    colors=[]
    symsize=[]
    for i in range(len(tt)):
        if (tt[i]==0):
            colors.append('red')
            symsize.append(20)
        elif (tt[i]==0.5):
            colors.append('#00FF00')
            symsize.append(20)
        else:
            colors.append('blue')
            symsize.append(20)

    phot=[]
    flux=[]
    stddev=[]
    meanflux=[]
    
    h=d[0].header
    
    mode = h['ESO INS MODE']
    if (mode.strip() =="HYBRID" ):
        if (h['ESO DET NAME'].strip()=="MATISSE-LM") :
            mode = "SIPHOT"
        else:
            mode = "HISENS"
    
    print(h['ESO INS MODE'])
    print(h['ESO DET NAME'])   
    print(mode)
    
    if mode== "HISENS":
        regions=[5]
    else :
        regions = [9,10,11,12,13]
    
    nregions=len(regions)
    
    for iregion in regions:
        photi=d['IMAGING_DATA'].data['DATA{0}'.format(iregion)]*1.0
        if doSky:
            skyi=dsky['IMAGING_DATA'].data['DATA{0}'.format(iregion)]*1.0
            meamskyi=np.mean(skyi,axis=0)
            photi=photi-meamskyi
        phot.append(photi)
        
        fluxi=np.sum(np.sum(photi,axis=1),axis=1)
        if unchop:
            fluxi= np.take(fluxi,targetframes)
        flux.append(fluxi)
        stddev.append(np.std(fluxi))
        meanflux.append(np.mean(fluxi))       
    print(np.array(stddev)/np.array(meanflux))
    print(meanflux)
    
    col=["violet","skyblue","salmon","peru","olive"]
    if mode == "SIPHOT":
        #renormalization of the interferometric channel
        mult=[1,1,0.23,1,1]
    else:
        mult=[1]
    for i in range(nregions):
        flux[i]=flux[i]*mult[i]
        plt.plot(flux[i],marker=".",color=col[i])  
        
        
    mini=np.min(flux)
    maxi=np.max(flux)
    x=np.arange(len(tt))
    y=tt*(maxi-mini)+mini
    plt.scatter(x,y, color=colors)

    plt.show()     
    
    

if  __name__== '__main__' :


    arg=sys.argv

    if (arg[1]=="--help" or arg[1]== "-h"):
        print "mat_fringe_flux_plot script to visualize flux vs time on fringes exposures"
        print "Usage : filename [-options]"
        print "options :"
        print " --sky skyfilename"
        print " --unchop"   
    else:   
        filename=arg[1]
        

        narg=len(arg)
        unchop=False
        skyfilename=""
        #print arg
        for i in range(2,narg):
            if (arg[i] == '--sky' ):
                skyfilename=arg[i+1]
            if (arg[i] == '--unchop' ):
                unchop=True    
                
        mat_time_flux_plot(filename,skyfilename,unchop)
                

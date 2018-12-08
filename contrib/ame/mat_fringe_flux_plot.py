import sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np



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
    doSky=False
    unchop=False
    #print arg
    for i in range(2,narg):
        if (arg[i] == '--sky' ):
            skyfilename=arg[i+1]
            doSky=True  
        if (arg[i] == '--unchop' ):
            unchop=True
            doSky=True      
    #print "Analyzing flux in fringe file={0} ".format(filename)
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
    for iregion in [9,10,12,13]:
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
    
    for i in range(4):
        plt.plot(flux[i])   
    plt.show()
    
	    


    """for i in range(4):
		print(fluxsum[i,:])
		plt.plot(fluxsum[i,:],color="k")
		plt.scatter(range(ni),fluxsum[i,:],marker="o",color=colorssum,)"""
    
    #i0 = 26
    #ni = 52
    
    i0 = 1
    ni = 14
    fluxsum=np.zeros([4,ni])
    colorssum=range(ni)
    for i in range(len(flux[0])):
        i2 = (i-i0) % ni
        for j in range(4):
            fluxsum[j,i2]+=flux[j][i]
            colorssum[i2]=colors[i]
    
    
    plt.plot(np.sum(fluxsum,axis=0),color="k")
    plt.scatter(range(ni),np.sum(fluxsum,axis=0),marker="o",color=colorssum)
    
    
    
    plt.show()
    
    

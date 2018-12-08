import sys
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import numpy as np



def mat_fringe_flux(filename,skyfilename):

    d=fits.open(filename)

   
    dsky=fits.open(skyfilename)

    tartype=d['IMAGING_DATA'].data['TARTYP']
    tt=(tartype=='T')*1+(tartype=='U')*0.5
    targetframes=np.argwhere(tt==1)
    skyframes=np.argwhere(tt==0)
    

    phot=[]
    flux=[]
    stddev=[]
    meanflux=[]
    for iregion in [9,10,12,13]:
        photi=d['IMAGING_DATA'].data['DATA{0}'.format(iregion)]*1.0
        skyi=dsky['IMAGING_DATA'].data['DATA{0}'.format(iregion)]*1.0
        meamskyi=np.mean(skyi,axis=0)
        photi=photi-meamskyi
        phot.append(photi)
        
        fluxi=np.sum(np.sum(photi,axis=1),axis=1)
        
        fluxi= np.take(fluxi,targetframes)
        flux.append(fluxi)
        stddev.append(np.std(fluxi))
        meanflux.append(np.mean(fluxi))       
    
    
    
    return {"mean" : meanflux, "var" : np.array(stddev)/np.array(meanflux)}
    


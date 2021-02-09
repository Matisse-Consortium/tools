# -*- coding: utf-8 -*-
"""
Created on Thu Jun 11 07:21:49 2020

@author: ame
"""

from astroquery.vizier import Vizier
from astroquery.simbad import Simbad
import numpy as np
from astropy.io import fits
#import stellarTools as st
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.modeling import models
mas2rad=np.pi/180./3600./1000.



###############################################################################
def mag2Flux(mag,band,unit,struct=True,typ="Johnson",zeropoint=None,wl=None):

    c=3e8

    if band=="Ks":
        typ="2MASS"
    elif band in ['Qs','N0','N1','N2','N3','N4','N5','Q0','Q1','Q2','Q3','Q4','Q5']:
        typ="MERLIN"
    elif band=="G":
        typ="GAIA"

    if zeropoint and wl:
        F=1.*zeropoint*10**(-mag/2.5)
    else:

        if typ=="Johnson":
            Bands=np.array(['U','B','V','R','I','J','H','K','L','M','N','Q'])
            wls=np.array([.36,.44,.55,.71,.97,1.25,1.6,2.22,3.54,4.8,10.6,21.])*1e-6
            F0=np.array([1823,4130,3781,2941,2635,1603,1075,667,288,170,36,9.4])
        elif typ=="2MASS":
            Bands=np.array(['J','H','Ks'])
            wls=np.array([1.235,1.662,2.159])*1e-6
            F0=np.array([1594,1024,666.7])
        elif typ== 'UKIRT':
            Bands=np.array(['V','I','J','H','K','L',"L'",'M','N','Q'])
            wls=np.array([0.5556,0.9,1.25,1.65,2.2,3.45,3.8,4.8,10.1,20])*1e-6
            F0=np.array([3540,2250,1600,1020,657,290,252,163,39.8,10.4])
        elif typ=='MIRLIN':
            Bands=np.array(['K','M','N0','N1','N2','N3','N4','N5','N','Qs','Q0','Q1','Q2','Q3','Q4','Q5'])
            wls=np.array([2.2,4.68,7.91,8.81,9.69,10.27,11.7,12.49,10.79,17.9,17.2,17.93,18.64,20.81,22.81,24.48])*1e-6
            F0=np.array([650,165,60.9,49.4,41.1,36.7,28.5,25.1,33.4,12.4,13.4,12.3,11.4,9.2,7.7,6.7])
        elif typ=='GAIA':
            Bands=np.array(['G'])
            wls=np.array([0.5857])*1e-6
            F0=np.array([2861])

        i=np.where(Bands==band)[0]
        if len(i)!=0:
            i=i[0]
            F=F0[i]*10**(-mag/2.5)
            wl=wls[i]
        else:
            print("Error could not find band {0}".format(band))
            return np.nan

    if unit=='W/m2/Hz':
        F=F*1e-26
    elif unit=='W/m2/m':
        F=F*1e-26*c/wl**2
    elif unit=='W/m2/micron':
        F=F*1e-26*c/wl**2/1e6
    elif unit=='W/cm2/micron':
        F=F*1e-26*c/wl**2/1e6/1e4
    elif unit!="Jy":
        return None

    if struct:
        return {"wl":wl,"f":F}
    else:
        return F


###############################################################################

# Kurucz not implememented as it uses Anthony's stellarTools library
"""
def _kuruczFitloglog(wl,Teff,diam):
    k=st.kurucz.star(diam,Teff,geff=4)
    flx=(k.f*u.W/u.m**2/u.micron).to(u.Jy,equivalencies=u.spectral_density(k.wl*u.m)).value
    #flx=st.fluxConverter(k.f,"W/m2/micron","Jy",wl=k.wl)["f"]
    with np.errstate(divide='ignore'):
        ll=np.interp(np.log10(wl),np.log10(k.wl),np.log10(flx))
    return 10**ll
"""


def _blackbodyFit(wl,Teff,diam):
    #bb=st.blackBodyStarAngular(wl,Teff,diam)
    #flx=(bb*u.W/u.m**3).to(u.Jy,equivalencies=u.spectral_density(wl*u.m)).value
    #flx=st.fluxConverter(bb,"W/m2/m","Jy",wl=wl)["f"]    
    bbfunc = models.BlackBody(temperature=Teff*u.K)
    flx=(bbfunc(wl*u.m)*np.pi*(diam/2*u.mas)**2).to(u.Jy).value  
    return flx

###############################################################################

def mat_queryCruzalebes(star,where="cds"):
    viz = Vizier(columns=['med-Lflux','disp-Lflux','med-Mflux','disp-Mflux',
                          'med-Nflux','disp-Nflux','Teff-GAIA','UDDL-est','e_LDD-est'])
    res2=viz.query_object(star,catalog='II/361/mdfc-v10')
    wl=[3.54e-6,4.8e-6,10.6e-6]
    try:
        fJy=[res2[0]['med-Lflux'][0],res2[0]['med-Mflux'][0],res2[0]['med-Nflux'][0]]
        fJy_err=[res2[0]['disp-Lflux'][0],res2[0]['disp-Mflux'][0],res2[0]['disp-Nflux'][0]]
        uddl=res2[0]['UDDL-est'][0]
        ldd_err=res2[0]['e_LDD-est'][0]
        TeffGaia=res2[0]['Teff-GAIA'][0]

    except:
        fJy=[np.nan,np.nan,np.nan]
        fJy_err=[0,0,0]

    return {"wl":wl,"flx":fJy,"flx_err":fJy_err,"uddl":uddl,"ldd_err":ldd_err,"TeffGaia":TeffGaia}

###############################################################################

def mat_querySimbad(star):
    simbd=Simbad()

    filts=['U','B','V','G','R','I','J','H','K']
    for filt in filts:
        simbd.add_votable_fields('flux({0})'.format(filt))
    simbd.add_votable_fields('sptype')
    simbd.add_votable_fields('plx')
    simbd.add_votable_fields('plx_error')
    res=simbd.query_object(star)
    spType=res['SP_TYPE'][0]
    plx=res['PLX_VALUE'][0]
    plx_err=res['PLX_ERROR'][0]
    dist=1000./plx
    dist_err=1000/plx**2*plx_err
    mag=[]
    for filt in filts:
        mag.append(res['FLUX_{0}'.format(filt)][0])
    return {"filters":filts,"mag":mag,"sptype":spType,"plx":plx,"plx_err":plx_err,"dist":dist,"dist_err":dist_err}



"""
mat_createTheoreticalSED :
    Function to create a theoretical SED for MATISSE calibrators.
    Return the star theoretical/fitted SED in Jy
    star      : The name of the star for which the theoretical SED should be created
    fluxModel : can be "bb" for blackbody or "kz" for kurucz models
    fit       : Should the Theoretical SED be fitted using magnitude or flux
                measurements from SIMBAD and Cruzalebes
    Teff      : Value for the Teff or source : Can be a number in Kelvin, "GAIA"
                or "spType"
    diam      : Either a number for the angular diameter in mas or "JSDC"
    filename  : if given Name of the created fits file fo
    query     : local for quering "local" catalogs or "cds" for quering CDS

"""

def mat_createTheoreticalSED(star,fluxModel="bb",fit=False,Teff="Gaia",
                             diam="JSDC",wl=None,filename=None,query="cds"):

    doQueryCruz=False
    doQuerySimbad=False
    if diam=="JSDC":
        doQueryCruz=True

    if Teff=="Gaia":
        doQueryCruz=True
    elif Teff=="spType":
        doQuerySimbad=True

    if fit==True:
        doQueryCruz=True
        doQuerySimbad=True

    if doQueryCruz==True:
        cruzData=mat_queryCruzalebes(star,where=query)

    if doQuerySimbad==True:
        simbadData=mat_querySimbad(star)

    if Teff=="Gaia":
        Teff=cruzData["TeffGaia"]
    elif Teff=="spType":
      print("Error : estimation of temperature from stellar Type not implemented yet")
      return None
      #starInfo=st.typicalStar(simbadData["sptype"])
      #Teff=starInfo.Teff

    if diam=="JSDC":
        diam=cruzData["uddl"]

    """
    if type(wl)==type(None):
        data=st.kurucz.star(diam,Teff,4.5)
        wl=data.wl
        flx=data.f*1e6 #flux given by kurucz is in W/m2/micron => W/m2/m    
        flx=(flx*u.W/u.m**3).to(u.Jy,equivalencies=u.spectral_density(wl*u.m)).value
        #flx=st.fluxConverter(flx,"W/m2/m","Jy",wl=wl)["f"]
    else:
        flx=_kuruczFitloglog(wl,Teff,diam)
        
    """
    if type(wl)==type(None):
        wl=np.linspace(0.1,20, num=1000)*1e-6
    
    if fluxModel=="bb":
        flx=_blackbodyFit(wl,Teff,diam)
    else:
        print("Error : Only blackbody stellar flux implemented yet")
        return None    


    col1 = fits.Column(name='WAVELENGTH', format='D', array=wl,unit='m')
    col2 = fits.Column(name='FLUX', format='D', array=flx,unit='Jy')
    cols = fits.ColDefs([col1, col2])
    hdu = fits.BinTableHDU.from_columns(cols)
    hdu.header['TEFF']=Teff
    hdu.header['DIAMETER']=diam
    hdu.header['MODEL']=fluxModel
    hdu.name="FLUX_THEORETICAL"
    primary_hdu = fits.PrimaryHDU()
    hdul = fits.HDUList([primary_hdu,hdu])

    if fit==True:
        mes_wl=cruzData["wl"]
        mes_flx=cruzData["flx"]
        mes_flx_err=cruzData["flx_err"]

        for i in range(len(simbadData["filters"])):
            if(not(np.ma.is_masked(simbadData["mag"][i]))):
                r=mag2Flux(simbadData["mag"][i],simbadData["filters"][i],"Jy")
                mes_wl.append(r["wl"])
                mes_flx.append(r["f"])
                if np.isnan(r["f"]):
                    mes_flx_err.append(0)
                else:
                    mes_flx_err.append(r["f"]/10.)


        col1 = fits.Column(name='WAVELENGTH',format='D',array=mes_wl,unit="m")
        col2 = fits.Column(name='FLUX', format='D', array=mes_flx,unit='Jy')
        col3 = fits.Column(name='FLUXERR', format='D', array=mes_flx_err,unit='Jy')
        cols = fits.ColDefs([col1, col2,col3])
        hdu2 = fits.BinTableHDU.from_columns(cols)
        hdu2.name="FLUX_MEASURED"
        hdul.append(hdu2)



        if fluxModel=="bb":
            fitloglogFunct=_blackbodyFit
        else:
            print("Error : Only blackbody stellar flux implemented yet")
            return None
            #fitloglogFunct= _kuruczFitloglog

        popt, pcov = curve_fit(fitloglogFunct, mes_wl, mes_flx,sigma=mes_flx_err,p0=[Teff,diam])
        diam_fitted=popt[1]
        diam_fitted_err=np.sqrt(pcov[1,1])
        Teff_fitted=popt[0]
        Teff_fitted_err=np.sqrt(pcov[0,0])

        flx_fitted=fitloglogFunct(wl,Teff_fitted,diam_fitted)
        col1 = fits.Column(name='WAVELENGTH',format='D',array=wl,unit="m")
        col2 = fits.Column(name='FLUX', format='D', array=flx_fitted,unit='Jy')
        cols = fits.ColDefs([col1, col2])
        hdu3 = fits.BinTableHDU.from_columns(cols)
        hdu3.header['TEFF']=Teff_fitted
        hdu3.header['TEFF_ERR']=Teff_fitted_err
        hdu3.header['DIAMETER']=diam_fitted
        hdu3.header['DIAMERR']=diam_fitted_err
        hdu3.header['MODEL']=fluxModel
        hdu3.name="FLUX_FITTED"
        hdul.append(hdu3)


    if filename!=None:
        hdul.writeto(filename,overwrite=True)

    return hdul

"""
mat_calibrateCorrFlux :
    Function to create a theoretical SED for MATISSE calibrators.
    Return the modified oifits with corrflux calibrated in Jy
        
    oifits    : the oifits for which corr flux should be calibrated
    fluxModel : can be "bb" for blackbody or "kz" for kurucz models
    fit       : Should the Theoretical SED be fitted using magnitude or flux
                measurements from SIMBAD and Cruzalebes
    Teff      : Value for the Teff or source : Can be a number in Kelvin, "GAIA"
                or "spType"
    diam      : Either a number for the angular diameter of the calibrator in mas or "JSDC"
    filename  : if given save the result to a fits file with this name
    query     : local for quering "local" catalogs or "cds" for quering CDS
"""


def mat_calibrateCorrFlux(oifits,fluxModel="bb",fit=False,Teff="Gaia",diam="JSDC",filename=None,filenameThFlux=None,query="cds"):

    if type(oifits)==type(""):
        oifits=fits.open(oifits)

    if "ESO PRO JSDC CAL1" in oifits[0].header.keys():
        calname=oifits[0].header["ESO PRO JSDC CAL1"].split("---")[0]
        print(calname)
    else:
        print("Not a calibrated MATISSE OIfits file abort")
        return

    wl=oifits["OI_WAVELENGTH"].data["EFF_WAVE"]

    dfluxCal=mat_createTheoreticalSED(calname,fluxModel=fluxModel,fit=fit,Teff=Teff,diam=diam,wl=wl,filename=filenameThFlux,query=query)

    if fit==False:
        fCal=dfluxCal["FLUX_THEORETICAL"].data["FLUX"]
    else:
        fCal=dfluxCal["FLUX_FITTED"].data["FLUX"]

    extnames= [hdui.name for hdui in oifits]
        

    if "OI_VIS" in extnames:
        if oifits["OI_VIS"].header["AMPTYP"]=="correlated flux":
            nB=np.shape(oifits["OI_VIS"].data["VISAMP"])[0]
            for iB in range(nB):
                oifits["OI_VIS"].data["VISAMP"][iB,:]*=fCal
                oifits["OI_VIS"].data["VISAMPERR"][iB,:]*=fCal  #TODO add error of flux calibration
            oifits["OI_VIS"].columns["VISAMP"].unit="Jy"
            oifits["OI_VIS"].columns["VISAMPERR"].unit="Jy"            
            
            for hdu in dfluxCal:
                append=False
                if hdu.name=="FLUX_THEORETICAL":
                    hdu.name="FCAL_T"
                    append=True
                elif  hdu.name=="FLUX_MEASURED":
                    hdu.name="FCAL_M"
                    append=True
                elif  hdu.name=="FLUX_FITTED":
                    hdu.name="FCAL_F"
                    append=True
                if append==True:
                    oifits.append(hdu)
              
    if filename!=None:
        oifits.writeto(filename,overwrite=True)

    return oifits

################################################################################


def mat_calibrateTotalFlux(oifitsSci,oifitsCal,fluxModel="bb",fit=False,Teff="Gaia",diam="JSDC",filename=None,filenameThFlux=None,query="cds",avgTel=False):

    if type(oifitsSci)==type(""):
        oifitsSci=fits.open(oifitsSci)
        
    if type(oifitsCal)==type(""):
        oifitsCal=fits.open(oifitsCal)
        
        
    calname=oifitsCal[0].header["ESO OBS TARG NAME"]
   
    
    wl=oifitsSci["OI_WAVELENGTH"].data["EFF_WAVE"]

    dfluxCal=mat_createTheoreticalSED(calname,fluxModel=fluxModel,fit=fit,Teff=Teff,diam=diam,wl=wl,filename=filenameThFlux,query=query)

    if fit==False:
        fCal=dfluxCal["FLUX_THEORETICAL"].data["FLUX"]
    else:
        fCal=dfluxCal["FLUX_FITTED"].data["FLUX"]

    extnamesSci= [hdui.name for hdui in oifitsSci]
    extnamesCal= [hdui.name for hdui in oifitsCal]
    
    if "OI_FLUX" in extnamesSci and "OI_FLUX" in extnamesCal:
        nTel=np.shape(oifitsSci["OI_FLUX"].data["FLUXDATA"])[0]  
        calflux=oifitsCal["OI_FLUX"].data["FLUXDATA"] 
        #calfluxErr=oifitsCal["OI_FLUX"].data["FLUXERR"] 
        for iTel in range(nTel):
            oifitsSci["OI_FLUX"].data["FLUXDATA"][iTel,:] = oifitsSci["OI_FLUX"].data["FLUXDATA"][iTel,:]/calflux[iTel,:]*fCal
            #TODO add error of flux calibration
            oifitsSci["OI_FLUX"].data["FLUXERR"][iTel,:] = oifitsSci["OI_FLUX"].data["FLUXDATA"][iTel,:]/calflux[iTel,:]*fCal 
    
        oifitsSci["OI_FLUX"].columns["FLUXDATA"].unit="Jy"
        oifitsSci["OI_FLUX"].columns["FLUXERR"].unit="Jy"           
   
    
        if avgTel==True:
            flxerr=np.sqrt(np.std(oifitsSci["OI_FLUX"].data["FLUXDATA"],axis=0)**2+np.mean(oifitsSci["OI_FLUX"].data["FLUXERR"],axis=0)**2)
            flx=np.median( oifitsSci["OI_FLUX"].data["FLUXDATA"],axis=0)

            cols=[]
            for col in oifitsSci["OI_FLUX"].columns:
                if col.name=="FLUXDATA":
                    newcol=fits.Column(name='FLUXDATA',format='D',array=flx,unit="Jy")
                    cols.append(newcol)
                elif col.name=="FLUXERR":
                    newcol=fits.Column(name='FLUXERR',format='D',array=flxerr,unit="Jy")
                    cols.append(newcol)
                else:
                    cols.append(col)
            newhdu = fits.BinTableHDU.from_columns(fits.ColDefs(cols)) 
            newhdu.header= oifitsSci["OI_FLUX"].header
            oifitsSci["OI_FLUX"]=newhdu
                

        
        for hdu in dfluxCal:
            append=False
            if hdu.name=="FLUX_THEORETICAL":
                hdu.name="FCAL_T"
                append=True
            elif  hdu.name=="FLUX_MEASURED":
                hdu.name="FCAL_M"
                append=True
            elif  hdu.name=="FLUX_FITTED":
                hdu.name="FCAL_F"
                append=True
            if append==True:
                oifitsSci.append(hdu)
              
    if filename!=None:
        oifitsSci.writeto(filename,overwrite=True)

    return oifitsSci


################################################################################

def mat_plotTheoreticalSED(hdul,ax=None,wlmin=0.35,wlmax=15,fmin=None,fmax=None):

    if not(ax):
        fig,ax=plt.subplots()

    wl_th=hdul["FLUX_THEORETICAL"].data["WAVELENGTH"]*1e6
    f_th=hdul["FLUX_THEORETICAL"].data["FLUX"]
    TEFF_th=hdul["FLUX_THEORETICAL"].header['TEFF']
    DIAM_th=hdul["FLUX_THEORETICAL"].header['DIAMETER']

    ax.loglog(wl_th,f_th,label="Theoretical : {0:.0f}K {1:.2f}mas".format(TEFF_th,DIAM_th))

    if  "FLUX_MEASURED" in [ hdu.name for hdu in hdul]:
        wl_meas=hdul["FLUX_MEASURED"].data["WAVELENGTH"]*1e6
        f_meas=hdul["FLUX_MEASURED"].data["FLUX"]
        f_meas_err=hdul["FLUX_MEASURED"].data["FLUXERR"]

        ax.errorbar(wl_meas,f_meas,yerr=f_meas_err,linestyle="",label="Measurements")

    if  "FLUX_FITTED" in [ hdu.name for hdu in hdul]:
        wl_fit=hdul["FLUX_FITTED"].data["WAVELENGTH"]*1e6
        f_fit=hdul["FLUX_FITTED"].data["FLUX"]
        TEFF_fit=hdul["FLUX_FITTED"].header['TEFF']
        TEFF_fit_err=hdul["FLUX_FITTED"].header['TEFF_ERR']
        DIAM_fit=hdul["FLUX_FITTED"].header['DIAMETER']
        DIAM_fit_err=hdul["FLUX_FITTED"].header['DIAMETER_ERR']

        ax.plot(wl_fit,f_fit,label="Fitted : {0:.0f}$\pm${1:.0f}K {2:.2f}$\pm${3:.2f}mas".format(TEFF_fit,TEFF_fit_err,DIAM_fit,DIAM_fit_err))

    ax.set_xlabel("$\lambda$ ($\mu$m)")
    ax.set_ylabel("Flux (Jy)")

    ax.set_xlim([wlmin,wlmax])

    if (fmin!=None) & (fmax!=None):
         ax.set_ylim([fmin,fmax])

    ax.legend()
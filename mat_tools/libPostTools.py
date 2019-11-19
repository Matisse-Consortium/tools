# -*- coding: utf-8 -*-
"""
Created on Tue Nov 19 13:50:34 2019

@author: ame
"""

from astropy.io import fits
import numpy as np
import os


BCD=[[0,1,2,3,4,5],  #OUT-OUT (0)
     [0,1,4,5,2,3],  #OUT-IN  (1)
     [0,1,3,2,5,4],  #IN-OUT  (2)
     [0,1,5,4,3,2]]  #IN-IN  (3)

BCDsign=[[1,1,1,1,1,1],   #OUT-OUT (0)
         [1,-1,1,1,1,1],  #OUT-IN  (1)
         [-1,-1,1,1,1,1], #IN-OUT  (2) => Pourquoi?
         [-1,-1,1,1,1,1]] #IN-IN  (3)

BCDcp=[[0, 1, 2, 3], #OUT-OUT (0)
       [3, 1, 2, 0], #OUT-IN  (1)
       [0, 2, 1, 3], #IN-OUT  (2)
       [3, 2, 1, 0]] #IN-IN  (3)


sta_index_cp=[[1,2,3],
           [0,1,2],
           [0,1,3],
           [0,2,3]]

uv1=[2,1,1,4]
uv2=[0,2,3,0]

BCDcpsign=[[ 1,  1,  1,  1], #OUT-OUT (0)
           [1,  -1,  -1, 1], #OUT-IN  (1)
           [ -1, 1, 1,  -1], #IN-OUT  (2)
           [-1, -1, -1, -1]] #IN-IN  (3)






################### Sorting oifits files by TPLSTART ##########################


def mat_sortByTplStart(oifitsList):
    data=[]
    if type(oifitsList[0])==type(""):
        data=[fits.open(oifitsi) for oifitsi in oifitsList]
    else:
        data=oifitsList

    tplstart=[d[0].header["ESO TPL START"] for d in data]

    nfiles=len(data)
    tplstartList=[]
    sortData=[]
    for i in range(nfiles):

        if not(tplstart[i] in tplstartList):
            sortData.append([])
            tplstartList.append(tplstart[i] )
        idx=np.where(np.array(tplstartList)==tplstart[i])[0][0]
        sortData[idx].append(data[i])
    return tplstartList,sortData




################### mMerging a list of oifits files ###########################


def mat_mergeOifits(oifitsList):
    data=[]
    if type(oifitsList[0])==type(""):
        data=[fits.open(oifitsi) for oifitsi in oifitsList]
    else:
        data=oifitsList

    nfile=len(data)

    for datai in data:
        mat_removeBCD(datai)
    #------------------- preparing merged container-------------

    nB=np.array([len(datai["OI_VIS2"].data) for datai in data])
    nBmin=np.min(nB)
    idx0=np.where(nB == nBmin)[0][0]
    avgFits=fits.HDUList([hdu.copy() for hdu in data[idx0]])

    #------------------------------OI_VIS2------------------------------------
    nB=np.array([len(datai["OI_VIS2"].data) for datai in data])
    nBmin=np.min(nB)
    idx0=np.where(nB == nBmin)[0][0]
    temp=data[idx0]["OI_VIS2"].copy()

    #mean of the square of vis2 to compute std(vis2) = sqrt(|<vis2>^2-<vis2^2>|)
    vis22=temp.data["VIS2DATA"]**2
    norm=1
    for ifile in range(nfile):
        nmod=nB[ifile]/nBmin
        if ifile!=idx0:
            for imod in range(nmod):
                for key in ["VIS2DATA","UCOORD","VCOORD","TIME","MJD","INT_TIME"]:
                    if len(np.shape(temp.data[key]))==2:
                        temp.data[key]= (temp.data[key]*norm + data[ifile]["OI_VIS2"].data[key][imod*nBmin:(imod+1)*nBmin,:])/(norm+1)
                    else:
                        temp.data[key]= (temp.data[key]*norm + data[ifile]["OI_VIS2"].data[key][imod*nBmin:(imod+1)*nBmin])/(norm+1)
                vis22 = (vis22*norm + data[ifile]["OI_VIS2"].data["VIS2DATA"][imod*nBmin:(imod+1)*nBmin,:]**2)/(norm+1)
                norm+=1
    temp.data["VIS2ERR"]=np.sqrt(np.abs(vis22- temp.data["VIS2DATA"]**2))/np.sqrt(norm)

    avgFits["OI_VIS2"]=temp

    #----------------------------OI_VIS---------------------------------------

    nB=np.array([len(datai["OI_VIS"].data) for datai in data])
    nBmin=np.min(nB)
    idx0=np.where(nB == nBmin)[0][0]
    temp=data[idx0]["OI_VIS"].copy()

    #mean of the square of visamp and visphi to compute the std

    viscompl=temp.data["VISAMP"]*np.exp(np.complex(0,1)*np.deg2rad(temp.data["VISPHI"]))

    #visamp2=temp.data["VISAMP"]**2
    #visphi2=temp.data["VISPHI"]**2

    norm=1
    for ifile in range(nfile):
        nmod=nB[ifile]/nBmin
        if ifile!=idx0:
            for imod in range(nmod):
                for key in ["VISAMPERR","VISPHIERR","UCOORD","VCOORD","TIME","MJD","INT_TIME"]:
                    if len(np.shape(temp.data[key]))==2:
                        temp.data[key]= (temp.data[key]*norm + data[ifile]["OI_VIS"].data[key][imod*nBmin:(imod+1)*nBmin,:])/(norm+1)
                    else:
                        temp.data[key]= (temp.data[key]*norm + data[ifile]["OI_VIS"].data[key][imod*nBmin:(imod+1)*nBmin])/(norm+1)
                #visamp2 = (visamp2*norm + data[ifile]["OI_VIS"].data["VISAMP"][imod*nBmin:(imod+1)*nBmin,:]**2)/(norm+1)
                #visphi2 = (visphi2*norm + data[ifile]["OI_VIS"].data["VISPHI"][imod*nBmin:(imod+1)*nBmin,:]**2)/(norm+1)
                visampi =data[ifile]["OI_VIS"].data["VISAMP"][imod*nBmin:(imod+1)*nBmin,:]
                visphii =data[ifile]["OI_VIS"].data["VISPHI"][imod*nBmin:(imod+1)*nBmin,:]
                viscompl = (viscompl*norm + visampi*np.exp(np.complex(0,1)*visphii))/(norm+1)

                norm+=1
    temp.data["VISAMP"]=np.abs(viscompl)
    temp.data["VISPHI"]=np.rad2deg(np.angle(viscompl))
    temp.data["VISPHIERR"]/=np.sqrt(norm)  # no better estimation than that for now
    temp.data["VISAMPERR"]/=np.sqrt(norm)  # no better estimation than that for now

    #temp.data["VISAMPERR"]=np.sqrt(np.abs(visamp2- temp.data["VISAMP"]**2))/np.sqrt(norm)
    #temp.data["VISPHIERR"]=np.sqrt(np.abs(visphi2- temp.data["VISPHI"]**2))/np.sqrt(norm)

    avgFits["OI_VIS"]=temp

    #-----------------------------OI_T3----------------------------------------

    nB=np.array([len(datai["OI_T3"].data) for datai in data])
    nBmin=np.min(nB)
    idx0=np.where(nB == nBmin)[0][0]
    temp=data[idx0]["OI_T3"].copy()

    expt3phi=np.exp(np.complex(0,1)*np.deg2rad(temp.data["T3PHI"]))


    norm=1
    for ifile in range(nfile):
        nmod=nB[ifile]/nBmin
        if ifile!=idx0:
            for imod in range(nmod):
                for key in ["T3PHIERR","U1COORD","V1COORD","U2COORD","V2COORD","TIME","MJD","INT_TIME"]:
                    if len(np.shape(temp.data[key]))==2:
                        temp.data[key]= (temp.data[key]*norm + data[ifile]["OI_T3"].data[key][imod*nBmin:(imod+1)*nBmin,:])/(norm+1)
                    else:
                        temp.data[key]= (temp.data[key]*norm + data[ifile]["OI_T3"].data[key][imod*nBmin:(imod+1)*nBmin])/(norm+1)
                expt3phi +=  np.exp(complex(0,1)*np.deg2rad(data[ifile]["OI_T3"].data["T3PHI"][imod*nBmin:(imod+1)*nBmin,:]))
                norm+=1
    temp.data["T3PHI"]=np.rad2deg(np.angle(expt3phi))
    temp.data["T3PHIERR"]/=np.sqrt(norm)  # no better estimation than that for now


    avgFits["OI_T3"]=temp

    #-----------------------------OI_FLUX--------------------------------------

    #TODO TODO TODO TODO

    return avgFits




####################### removing BCD in an oifits files ########################



def mat_removeBCD(oifits,saveFits=False):
    if type(oifits)==type(""):
        data=fits.open(oifits)
    else:
        data=oifits
    bcd1=data[0].header["ESO INS BCD1 NAME"]
    bcd2=data[0].header["ESO INS BCD2 NAME"]
    bcd=(bcd2 == "IN") + 2*(bcd1 == "IN")
    if bcd==0:
        #print("no bcd")
        return
    #------------------OI_VIS2-------------------------
    temp=data["OI_VIS2"].copy()

    nB=len(temp.data)
    for iB in range(nB):
        iB2=iB % 6
        shift0 = (iB / 6)*6
        temp.data[iB2+shift0]=data["OI_VIS2"].data[BCD[bcd][iB2]+shift0]
        temp.data[iB2+shift0]["UCOORD"]*=BCDsign[bcd][iB2]
        temp.data[iB2+shift0]["VCOORD"]*=BCDsign[bcd][iB2]
        if BCDsign[bcd][iB2]==-1:
            temp.data[iB2+shift0]["STA_INDEX"]=np.flip(temp.data[iB2+shift0]["STA_INDEX"],axis=0)

    data["OI_VIS2"]=temp
    #------------------OI_VIS--------------------------
    temp=data["OI_VIS"].copy()

    nB=len(temp.data)
    for iB in range(nB):
        iB2=iB % 6
        shift0 = (iB / 6)*6
        temp.data[iB2+shift0]=data["OI_VIS"].data[BCD[bcd][iB2]+shift0]
        temp.data[iB2+shift0]["VISPHI"]*=BCDsign[bcd][iB2]
        temp.data[iB2+shift0]["UCOORD"]*=BCDsign[bcd][iB2]
        temp.data[iB2+shift0]["VCOORD"]*=BCDsign[bcd][iB2]
        if BCDsign[bcd][iB2]==-1:
            temp.data[iB2+shift0]["STA_INDEX"]=np.flip(temp.data[iB2+shift0]["STA_INDEX"],axis=0)

    data["OI_VIS"]=temp

    #------------------OI_T3---------------------------
    sta_index=data["OI_ARRAY"].data["STA_INDEX"]
    temp=data["OI_T3"].copy()

    nB=len(temp.data)
    for iB in range(nB):
        iB2=iB % 4
        shift0 = (iB / 4)*4
        data["OI_T3"].data[BCDcp[bcd][iB2]+shift0]["T3PHI"]*=BCDcpsign[bcd][iB2]
        temp.data[iB2+shift0]=data["OI_T3"].data[BCDcp[bcd][iB2]+shift0]
        temp.data[iB2+shift0]["U1COORD"]=data["OI_VIS2"].data["UCOORD"][uv1[iB2]]
        temp.data[iB2+shift0]["V1COORD"]=data["OI_VIS2"].data["VCOORD"][uv1[iB2]]
        temp.data[iB2+shift0]["U2COORD"]=data["OI_VIS2"].data["UCOORD"][uv2[iB2]]
        temp.data[iB2+shift0]["V2COORD"]=data["OI_VIS2"].data["VCOORD"][uv2[iB2]]
        temp.data[iB2+shift0]["STA_INDEX"]=np.array([sta_index[sta_index_cp[iB2][i]] for i in range(3)])

    data["OI_T3"]=temp

    #------------------OI_FLUX-------------------------

    #TODO TODO TODO

    #------------------END------------------------

    data[0].header["ESO INS BCD1 NAME"]="OUT"
    data[0].header["ESO INS BCD2 NAME"]="OUT"
    if saveFits==True:

        #filenamein=data.filename()
        #filenameout= filenamein.split(".fits")[0]+"_noBCD.fits"
        data.writeto(data.filename())


#=============================================================================


def mat_mergeByTplStart(something,save=False,verbose=True,dirOut="MERGED"):
    data=[]
    currentDir= os.path.abspath("./")
    if type(something)==type(""):
        #should be q directory but check first
        if not(os.path.isdir(something)):
            print("Error : {0} is not a directory".format(something))
            return
        currentDir=something
        something=[something+"/"+fi for fi in os.listdir(something) if ".fits" in fi]
    if type(something[0])==type(""):
        data=[fits.open(oifitsi) for oifitsi in something]
    else:
        data=something

    tplstart,sortedData=mat_sortByTplStart(data)
    ntpl=len(tplstart)

    if verbose:
        print("Number of TPLSTART : {0}".format(ntpl))

    mergedData=[]

    for itpl in range(ntpl):
        band=np.array([d[0].header["ESO DET NAME"] for d in sortedData[itpl]])
        idxL=np.where(band=="MATISSE-LM")[0]
        idxN=np.where(band=="MATISSE-N")[0]
        if verbose:
            print("******({0}/{1}) TPLSTART={2}******".format(itpl+1,ntpl,tplstart[itpl]))
            print("number of files to merge : {0} for LM and {1} for N".format(len(idxL),len(idxN)))
        for idxi in [idxL,idxN]:
            if len(idxi!=0):
                datai=[sortedData[itpl][idata] for idata in idxi]
                #filenames=[dataii.filename() for dataii in datai]
                #print(filenames)
                mergedi=mat_mergeOifits(datai)
                mergedData.append(mergedi)
                if save:
                    if not(os.path.exists(currentDir+"/"+dirOut)):
                        os.mkdir(currentDir+"/"+dirOut)

                    fileout=currentDir+"/"+dirOut+"/"+os.path.basename(sortedData[itpl][idxi[0]].filename()).replace("_OUT","").replace("_IN","").replace("_noChop","").replace("_Chop","")
                    if verbose:
                        print("Saving merged file to {0}".format(fileout))
                    mergedi.writeto(fileout,overwrite=True)




    return mergedData,data
























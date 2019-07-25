# -*- coding: utf-8 -*-
"""
General use routines for the OPD reduction of MATISSE
Created on Mon Aug 01 11:09:12 2016
@author: ame
"""

from   matplotlib import gridspec
import matplotlib.pyplot as plt
from   astropy.io import fits
import numpy as np
from   scipy.optimize import curve_fit
import glob,os

def fit_func(x, a0, a1, a2, a3):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 /2) + a3 
    return y

def MatisseLamFit(y_in,lam_in,y0,dimy):
        a   = np.polyfit(y_in,lam_in,2);
        y   = np.arange(dimy)+y0;
        lam = a[2]+a[1]*y+a[0]*y**2;
        return lam

class MatisseCubeData():
    def __init__(self,filename=None, band=None):
        if not filename==None:
            self.set_band(band);
            self.load(filename);
            
    def set_band(self,band):
        self.band = band;
        if band==1:
            self.sNx  = 16  # size noise x
            self.sPx  = 150 # size photo x
            self.sIx  = 625 # size interfx
            self.sNy  = 8   # size noise y
        else:
            self.sNx = 16  # size noise x
            self.sPx = 78  # size photo x
            self.sIx = 468 # size interfx
            self.sNy = 8   # size noise y

    def load(self,filename):
        try:
            fitsStruct  = fits.open(filename)
            self.data   = fitsStruct[0].data
            fitsStruct.close();
            self.nframe = self.data.shape[0]
            self.dimy   = self.data.shape[1]-2*self.sNy
            self.interf = self.data[:,self.sNy:self.sNy+self.dimy,self.sNx+2*self.sPx:self.sNx+2*self.sPx+self.sIx]
   
            # Debug plot of the data frame
            #plt.figure();
            #plt.imshow(self.interf[1,:,:],vmin=0);
            
            self.dimx   = self.interf.shape[2]
        except:
            print("Failed to open {0}".format(filename))  

    def prepare(self,lam,dimfx,dimfy,P0):
        if lam.size==self.dimy:
            self.dimfx        = dimfx
            self.dimfy        = dimfy
            self.P0           = P0
            self.currentFrame = 0
            #interpolation grid in the X and lam orientation
            self.lam       = lam.reshape(self.dimy,1).repeat(self.dimx,axis=1)
            self.sigma     = 1./self.lam          
            self.dsigm_min = self.sigma[1,0]- self.sigma[0,0]   # For L band 
            self.lsigm     = self.sigma[self.dimy-1,0]-self.sigma[0,0]
            self.nsigm     = int(2*np.fix(self.lsigm/self.dsigm_min/2))
            self.sigmap    = self.sigma[0,0]+self.dsigm_min*np.arange(self.nsigm)
            self.x         = np.arange(self.dimx)-(self.dimx-1)/2.
            self.x         = self.x.reshape(1,self.dimx).repeat(self.dimy,axis=0)
            self.xp        = self.x/self.lam
            self.dx0       = self.xp[0,1]-self.xp[0,0]
            self.dx1       = self.xp[self.dimy-1,1]-self.xp[self.dimy-1,0]             
            if(self.band==1):
                self.maxi  = max(self.xp[self.dimy-1,:])   #For L band 
            else:
                self.maxi  = max(self.xp[0,:])   #For N band 
            self.n         = np.fix(self.maxi*2/self.dx0)               #For L band
            self.xn        = (np.arange(self.n)/(self.n-1.)-0.5)*2*self.maxi
            
            #Conversion between pixels and piston in microns. d'apres koechlin et al
            self.lam_mean  = (np.max(self.lam)+np.min(self.lam))/2
            self.dlam      = (np.max(self.lam)-np.min(self.lam))
            self.pixToPist = self.lam_mean**2/self.dlam/self.dimfy*self.nsigm
            self.scale     = (np.arange(self.dimfy)-(self.dimfy-1.)/2.-0.5)*self.pixToPist
            
            #set the interpolation,fft, and fftCuts arrays
            self.interfpX     = np.zeros([self.nframe,self.dimy,self.n])
            self.interfpXY    = np.zeros([self.nframe,self.nsigm,self.n])       
            self.interfpFinal = np.zeros([self.nframe,self.dimfy,self.dimfx])
            self.fft2D        = np.zeros([self.nframe,self.dimfy,self.dimfx],dtype=np.complex)
            self.peakCut      = np.zeros([self.nframe,7,self.dimfy])
            self.piston       = np.zeros([self.nframe,7])
            self.pistonAvg    = np.zeros([self.nframe,7])
            self.pistonErr    = np.zeros([self.nframe,7])
            
            #compute the 7 peak position (photo + 6 interf)
            self.P            = self.dimfx/2+self.P0*(np.arange(7))*self.dimfx
        
    def interpolateXY(self):
        fringe = self.interf[self.currentFrame,:]
        for i in range(self.dimy):
            self.interfpX[self.currentFrame,i,:]  = np.interp(self.xn,self.xp[i,:],fringe[i,:],left=0,right=0)
            
        # Debug plot of the interpolated data frame
        #plt.figure();
        #plt.imshow(self.interfpX[self.currentFrame,:,:],vmin=0);
        
        
        for i in range(self.n):
            if(self.band==1): # L band
                self.interfpXY[self.currentFrame,:,i] = np.interp(self.sigmap,self.sigma[:,0],self.interfpX[self.currentFrame,:,i])
            else:             # N band
                self.interfpXY[self.currentFrame,::-1,i] = np.interp(self.sigmap[::-1],self.sigma[::-1,0],self.interfpX[self.currentFrame,::-1,i])
                        
        # Debug plot of the interpolated data frame
        #plt.figure();
        #plt.imshow(self.interfpXY[self.currentFrame,:,:],vmin=0,vmax=100);
        
        
        if (self.dimfx >= self.n):
            self.interfpFinal[self.currentFrame,self.dimfy/2-self.nsigm/2:self.dimfy/2+self.nsigm/2,self.dimfx/2-self.n/2:self.dimfx/2+self.n/2]=self.interfpXY[self.currentFrame,:,:]
        else:
            self.interfpFinal[self.currentFrame,self.dimfy/2-self.nsigm/2:self.dimfy/2+self.nsigm/2,:]=self.interfpXY[self.currentFrame,:,self.n/2-self.dimfx/2:self.n/2+self.dimfx/2]
  
        # Debug plot of the interpolated data frame
        #plt.figure();
        #plt.imshow(self.interfpFinal[self.currentFrame,:,:],vmin=0,vmax=100);
        
    def computeFft2D(self):
        # Compute a 2D FFT of the frame
        self.fft2D[self.currentFrame,:,:] = np.fft.fftshift(np.fft.fft2(self.interfpFinal[self.currentFrame,:,:]))
        self.fft2D[self.currentFrame,:,:] = self.fft2D[self.currentFrame,:,:]/np.max(abs(self.fft2D[self.currentFrame,:,:]))                    
                     
    def extractPeakCuts(self):
        for i in range(7):
            self.peakCut[self.currentFrame,i,:] = abs(self.fft2D[self.currentFrame,:,self.P[i]])                  
      
    def computePistons(self):
        for i in range(7):
            peak = self.peakCut[self.currentFrame,i,:]
            maxi = np.max(peak)
            w    = np.where(peak==maxi)[0][0]
            p0   = [maxi,self.scale[w],8,0.]
            try:
                parameters, covariance = curve_fit(fit_func, self.scale, peak,p0)
            except RuntimeError:
                print("Error - curve_fit failed")
                
            self.piston[self.currentFrame,i] = parameters[1]
    
    def compute(self):
        self.interpolateXY()
        self.computeFft2D()
        self.extractPeakCuts()
        self.computePistons()      
        
    def computeAll(self):
        for i in range(self.nframe):   
            self.currentFrame = i
            self.compute()
            
            # Plot the FFT2D
            #plt.figure();
            #plt.clf()
            #plt.imshow(np.log(np.abs(self.fft2D[1,:,:])))
            #plt.plot([self.P,self.P], [0, self.dimfx], 'r-', lw=1) # Red straight line
            #plt.plot([self.P[1],self.P[1]], [0, self.dimfx], 'g-', lw=1) # Red straight line
            #plt.plot([self.P[2],self.P[2]], [0, self.dimfx], 'b-', lw=1) # Red straight line
            
        self.pistonAvg = np.mean(self.piston,axis=0)
        self.pistonErr = np.std(self.piston,axis=0)

    def subtractDark(self,darkFile):
        # Load the dark file
        fitsStruct    = fits.open(darkFile)
        self.dark_data = fitsStruct[0].data
        fitsStruct.close();
        self.dark_nframe = self.dark_data.shape[0]
        self.dark_dimy   = self.dark_data.shape[1]-2*self.sNy
        self.dark_interf = self.dark_data[:,self.sNy:self.sNy+self.dimy,self.sNx+2*self.sPx:self.sNx+2*self.sPx+self.sIx]
        # Subtract the dark            
        self.interf = self.interf - np.median(self.dark_interf, axis=0);
        self.data   = self.data   - np.median(self.dark_data,   axis=0);
           
class MatisseMotorCalib():
    def __init__(self,dir,lam,dimfx,dimfy,P0,peakNum,band,darkFile=None):
        self.darkFile= darkFile
        self.dir     = dir
        self.peakNum = peakNum
        self.dimfx   = dimfx        
        self.dimfy   = dimfy
        self.P0      = P0
        self.lam     = lam
        div          = dir.split("/")
        self.dstep   = div[-2].split("_")[-2]
        self.device  = div[-3]
        self.pos     = div[-2].split("_")[-1]
    
        self.files   = glob.glob(self.dir+"*.fits")
    
        self.stepUp        = []
        self.pistonUp      = []
        self.pistonUpErr   = []
        self.stepDown      = []
        self.pistonDown    = []
        self.pistonDownErr = []
        
        self.data     = MatisseCubeData(self.files[0],band)
        
        self.data.prepare(self.lam,self.dimfx,self.dimfy,self.P0)
        dirsave       = os.path.dirname(self.dir[:-2])+'/'
        self.fileDown = dirsave+"MATISSE_MotorCalib_{0}_{1}_{2}_Down.dat".format(self.device,self.dstep,self.pos)
        self.fileUp   = dirsave+"MATISSE_MotorCalib_{0}_{1}_{2}_Up.dat".format(self.device,self.dstep,self.pos)
        
    def compute(self,num):
        self.files[num]
        self.data.load(self.files[num])
        
        # Remove the dark frames from the fringes frames
        if not self.darkFile==None:
           self.data.subtractDark(self.darkFile)
        
        # Compute the FFT of each file
        self.data.computeAll()
        
        file0 = self.files[num].split("/")[-1]
        step  = (file0.split("_")[-1]).split(".")[0]
        if (step.find("b")==-1):
            step = int(step)
            way  = "up"                        
        else:   
            step = int(step[:-1])
            way  = "down"  
                
        print("Piston({0}{1}) = {2:2.2f} +- {3:2.2f} microns".format(step,way,self.data.pistonAvg[self.peakNum],self.data.pistonErr[self.peakNum]))
        return (step,way,self.data.pistonAvg[self.peakNum],self.data.pistonErr[self.peakNum])
      
    def computeAll(self):
         for ifile in range(len(self.files)):
             step,way,pistonAvg,pistonErr=self.compute(ifile)
             if way=="up":
                self.stepUp.append(step)
                self.pistonUp.append(pistonAvg)
                self.pistonUpErr.append(pistonErr)            
                
             else:       
                self.stepDown.append(step)            
                self.pistonDown.append(pistonAvg)
                self.pistonDownErr.append(pistonErr)  
        
         #on ajoute les pts de départ et de fin qui sontpar défaut dans pistonUp à pistonDown
         mini = min(self.stepUp)
         maxi = max(self.stepUp)
         wmin = next(i for i in range(len(self.stepUp)) if self.stepUp[i]==mini)
         wmax = next(i for i in range(len(self.stepUp)) if self.stepUp[i]==maxi)
         self.stepDown.append(self.stepUp[wmin])            
         self.pistonDown.append(self.pistonUp[wmin])
         self.pistonDownErr.append(self.pistonUpErr[wmin])   
         self.stepDown.append(self.stepUp[wmax])            
         self.pistonDown.append(self.pistonUp[wmax])
         self.pistonDownErr.append(self.pistonUpErr[wmax])  
        
         #Conversion des listes en tableaux numpy
         self.stepDown      = np.array(self.stepDown)
         self.pistonDown    = np.array(self.pistonDown)
         self.pistonDownErr = np.array(self.pistonDownErr)
         self.stepUp        = np.array(self.stepUp)
         self.pistonUp      = np.array(self.pistonUp)
         self.pistonUpErr   = np.array(self.pistonUpErr)
         
    def save(self):
        #Ecriture des pistons dans un fichier
        tabDown = np.array([self.stepDown,self.pistonDown,self.pistonDownErr])
        np.savetxt(self.fileDown,tabDown)
        tabUp   = np.array([self.stepUp,self.pistonUp,self.pistonUpErr])
        np.savetxt(self.fileUp,tabUp)  
        print("Saving piston data for up to {0}".format(self.fileUp))
        print("Saving piston data for Down to {0}".format(self.fileDown))
        
    def plot(self,dirout):
        fileeps = MatisseMotorCalibPlot(self.fileUp,self.fileDown,dirout)
        return fileeps
        
def MatisseMotorCalibPlot(fileUp,fileDown,dirout):
    
    stepUp,pistonUp,pistonUpErr       = np.loadtxt(fileUp)
    stepDown,pistonDown,pistonDownErr = np.loadtxt(fileDown)
        
    div    = fileUp.split("/")[-1].split("_")
    dstep  = div[-3]
    device = div[-4]
    pos    = div[-2]    
    
    #Fit de la loi step/piston
    xUp           = stepUp#np.array([mini,maxi])
    aUp,covUp     = np.polyfit(stepUp,pistonUp,1,w=pistonUpErr,cov=True)
    yUp           = aUp[0]*xUp+aUp[1]
    xDown         = stepDown
    aDown,covDown = np.polyfit(stepDown,pistonDown,1,w=pistonDownErr,cov=True)
    yDown         = aDown[0]*xDown+aDown[1]
    
    fig  = plt.figure()
    gs   = gridspec.GridSpec(3, 1)
    plt1 = fig.add_subplot(gs[0:2,:])
    plt.title("Calibration of {0} with {1} steps around position {2}".format(device,dstep,pos)) 
    plt.ylabel(r"Measured piston ($\mu$m)")
    dot_up,    = plt.plot(stepUp,pistonUp, 'ro')
    line_up,   = plt.plot(xUp,yUp,'r')
    dot_down,  = plt.plot(stepDown,pistonDown, 'go')
    line_down, = plt.plot(xDown,yDown,'g')
    plt.legend([(line_up,dot_up), (line_down,dot_down)], ['Up $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/step'.format(aUp[0],np.sqrt(covUp[0,0])), 'Down $\Rightarrow${0:2.4f}$\pm${1:2.4f} $\mu$m/step'.format(aDown[0],np.sqrt(covDown[0,0]))],loc=3)
    plt.tick_params(axis='x',which='both', bottom='off',top='off',labelbottom='off')
    
    #plt1.set_xticklabels([]);
    fig.add_subplot(gs[2,:],sharex=plt1)
    plt.plot(xUp,yUp-pistonUp,'ro')
    plt.plot(xDown,yDown-pistonDown,'go')
    plt.plot([xDown[0],xDown[-1]],[0,0],'k--')
    plt.xlabel(r"Motor position (steps)")
    plt.ylabel(r"Residual ($\mu$m)")
   
    fileeps = dirout+"MATISSE_MotorCalib_{0}_{1}_{2}_Down.eps".format(device,dstep,pos)
    plt.savefig(fileeps)
    filepng = dirout+"MATISSE_MotorCalib_{0}_{1}_{2}_Down.png".format(device,dstep,pos)
    plt.savefig(filepng)
    return fileeps
    
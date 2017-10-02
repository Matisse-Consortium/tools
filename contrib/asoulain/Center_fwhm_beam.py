# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 16:48:42 2016

@author: asoulain
"""

import numpy as np
from matplotlib import pyplot as plt
import scipy.ndimage.measurements as mes
import warnings
from skimage import measure
import scipy.interpolate as ip
from Display_beam import Open_beam
from dpfit_red import leastsqFit
import pyfits
from mpl_toolkits.axes_grid1 import ImageGrid
from matplotlib.colors import PowerNorm

warnings.filterwarnings('ignore')

def gauss_2D_asym((x, y), param):
    """
    Creates 2D oriented gaussian with an asymmetrical grid.
    """
    amplitude = param['amplitude']
    xo = param['xo']
    yo = param['yo']
    sigma_x = param['sigma_x']
    sigma_y = param['sigma_y']
    theta = np.deg2rad(param['theta'])
    size_x = len(x)
    size_y = len(y)
    im = np.zeros([size_y, size_x])
    xo = float(xo)
    yo = float(yo)    
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    im =  amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) 
                            + c*((y-yo)**2)))
                            
    return im.flatten()
    
def gauss_1D(x, param):
    A = param['amplitude']
    xo = param['xo']
    sigma = param['sigma']
    g = A*np.exp(-((x-xo)**2)/(2*sigma**2))
    return g
    
    
plt.close('all')

#------------------------
#       Open files
#------------------------
band = 'L'
bcd = 'IN'
nbeam = 2

dic = Open_beam(band, nbeam, bcd)
beam = dic['beam']

datadir = '/Users/asoulain/Documents/These/DRS_MATISSE/DATA/focusIR_15112016/'

hdu = pyfits.open(datadir + 'CB%sDiaPhLmPoOp.fits'%nbeam)
hdu_dark = pyfits.open(datadir + 'Cdark.fits')

i = 3
mode = 'INTERF'
if mode == 'INTERF':
    xmin = 700
    xmax = 1420
    aspect = 24
else:
    if nbeam > 2:
        xmin = 1420
        xmax = -1
        aspect = 4
    else:
        xmin = 0
        xmax = 700
        aspect = 4

if band == 'L':
    raw_im = hdu[0].data[i,:,:][:, xmin:xmax]
    dark = hdu_dark[0].data[i,:,:][:, xmin:xmax]
    beam = raw_im - dark
    beam[beam < 0] = 0
else:
    raw_im = dic['lamp']

#--------------------------------------------------
#   FIT 1D (x and y) to initiate parameters for 2D    
#--------------------------------------------------

c_max = np.where(beam == np.max(beam))

x_prf = beam[c_max[0][0],:] #extract xprofil
y_prf = beam[:,c_max[1][0]] #extract yprofil
size_x = np.shape(beam)[1]
size_y = np.shape(beam)[0]
x = np.arange(size_x)
y = np.arange(size_y)

#Init parameters
param_x = {
'amplitude' : np.max(x_prf),
'xo'        : c_max[1][0],
'sigma'   : 3,
}

param_y = {
'amplitude' : np.max(y_prf),
'xo'        : c_max[0][0],
'sigma'   : 10,
}

#Compute uncertainties
err  = np.ones(len(x_prf))*np.std(x_prf)
err2 = np.ones(len(y_prf))*np.std(y_prf)

#Compute fit
fitx = leastsqFit(gauss_1D, x, param_x, x_prf, err = err, verbose = False, epsfcn=1e-12, ftol=1e-10 )
fity = leastsqFit(gauss_1D, y, param_y, y_prf, err = err2, verbose = False, epsfcn=1e-12, ftol=1e-10 )

#Save best value for 2D
sigx = fitx['best']['sigma']
sigy = fity['best']['sigma']

f_anamorph_1D = sigx/sigy

size_x = np.shape(beam)[1]
size_y = np.shape(beam)[0]

x = np.arange(size_x)
y = np.arange(size_y)
tup = np.meshgrid(x,y)

#New init parameters for 2D
param = {
'amplitude' : np.max(beam),
'xo'        : c_max[1][0],
'yo'        : c_max[0][0],
'sigma_x'   : sigx,
'sigma_y'   : sigy,
'theta'     : 5,
}

fit_2D = True #Compute or not 2D fit
if fit_2D == True:
    fit = leastsqFit(gauss_2D_asym, tup, param, beam.ravel(), err = np.std(beam[0:10, 0:10]).ravel(), verbose = False, epsfcn=1e-5, ftol=1e-7,
                     normalizedUncer=True)
    f_anamorph_2D = fit['best']['sigma_x']/fit['best']['sigma_y']
    err_f_anamorph = f_anamorph_2D * (fit['uncer']['sigma_x']/fit['best']['sigma_x'] + fit['uncer']['sigma_y']/fit['best']['sigma_y']) #Error come from covariance matrix

#Compute gaussian model for the best 2D fit
fit_2D = gauss_2D_asym(tup, fit['best']).reshape([size_y, size_x])

chi2_r = fit['chi2']

print '\n-----------------------'
print 'Gauss2D method:'
print '-----------------------'
print 'Centroid : [%2.1f, %2.1f], fwhmx = %2.1f, fwhmy = %2.1f, f = %2.2f ± %2.2f'%(fit['best']['yo'], fit['best']['xo'], 2.355*fit['best']['sigma_x'], 2.355*fit['best']['sigma_y'], f_anamorph_2D, err_f_anamorph)
print 'Orient = %2.3f deg'%np.rad2deg(fit['best']['theta'])
print 'Chi2 = %2.3f'%(chi2_r)


#--------------------------------------------------
#             1D interpolation method
#--------------------------------------------------
def fwhm_ip_1d(x, prf):
    aa = np.where(prf == np.max(prf))
    f_x1 = ip.interp1d(prf[0:int(aa[0])], x[0:int(aa[0])])
    f_x2 = ip.interp1d(prf[int(aa[0]):-1], x[int(aa[0]):-1])

    fwhm_ip = f_x2(.5*np.max(prf)) - f_x1(.5*np.max(prf))
    return fwhm_ip

p_min = 3
p_maj = 20
top = beam[c_max[0][0]-3:c_max[0][0]+3,c_max[1][0]-20:c_max[1][0]+20]
cm = mes.center_of_mass(top)

c_ip = [c_max[0][0] - (p_min - cm[0]), c_max[1][0] - (p_maj - cm[1])]

print '\n-----------------------'
print 'Interpolation method:'
print '-----------------------'
print 'Centroid : [%2.1f, %2.1f], fwhmx = %2.1f, fwhmy = %2.1f, f = %2.1f'%(c_ip[0],  c_ip[1], fwhm_ip_1d(x, x_prf), fwhm_ip_1d(y, y_prf), fwhm_ip_1d(x, x_prf)/fwhm_ip_1d(y, y_prf))

#--------------------------------------------------
#             2D moment method
#--------------------------------------------------
cond_seuil = beam >= .5*np.max(beam) #3%

im = np.zeros([size_y, size_x])

im[cond_seuil] = 1

beam_thres = beam*im

#Compute moment on threshold image
m = measure.moments(beam_thres, order = 3)
cr = m[0, 1]/m[0, 0]
cc = m[1, 0]/m[0, 0]

#Compute central moment on threshold image
mu = measure.moments_central(beam_thres, cr, cc, order = 3)
mu20 = mu[2,0]/mu[0,0]
mu02 = mu[0,2]/mu[0,0]
mu11 = mu[1,1]/mu[0,0]

theta = -0.5*np.arctan(2*mu11/(mu20-mu02))

l1 = ((mu20+mu02) + np.sqrt((mu20-mu02)**2 + 4*mu11**2.))/2.
l2 = ((mu20+mu02) - np.sqrt((mu20-mu02)**2 + 4*mu11**2.))/2.

a = 4*np.sqrt(l1)
b = 4*np.sqrt(l2)

fit_mom_x = gauss_1D(x, {'amplitude':np.max(beam), 'xo':cc, 'sigma':a/2.355})
fit_mom_y = gauss_1D(x, {'amplitude':np.max(beam), 'xo':cr, 'sigma':b/2.355})
fit_mom_2D = gauss_2D_asym(tup, {'amplitude':np.max(beam), 'xo':cc, 'yo':cr, 'sigma_x':a/2.355, 'sigma_y':b/2.355, 'theta':theta}).reshape([size_y, size_x])

print '\n-----------------------'
print 'Moments method :'
print '-----------------------'
print 'Centroid : [%2.1f, %2.1f], fwhmx = %2.1f, fwhmy = %2.1f, f = %2.1f'%(cr,  cc, a, b, a/b)
print 'Orient = %2.1f deg'%np.rad2deg(theta)

#-------------------------------------------------------------------------------
#  Display 1D profil to check if the fit are good (moment and gaussian)
#-------------------------------------------------------------------------------
if True:
    c_fitx = gauss_1D(x, fitx['best'])
    cr = int(cr)
    cc = int(cc)
    x = np.arange(size_x)
    xx = x[c_max[1][0]-200:c_max[1][0]+200]
    y = np.arange(size_y)
    yy = y[cr-10:cr+10]
    major_prf = beam[c_max[0][0],c_max[1][0]-200:c_max[1][0]+200]
    minor_prf = beam[cr-10:cr+10, cc]
    major_res = (beam[c_max[0][0],c_max[1][0]-200:c_max[1][0]+200]-fit_2D[c_max[0][0], c_max[1][0]-200:c_max[1][0]+200])
    minor_res = beam[cr-10:cr+10, cc] - fit_2D[cr-10:cr+10, cc]
    
    major_res_mom = major_prf-fit_mom_x[c_max[1][0]-200:c_max[1][0]+200]
    minor_res_mom = minor_prf-fit_mom_y[cr-10:cr+10]

    fig = plt.figure(0,figsize = (14,8))
    ax1 = plt.subplot2grid((3,4), (0,0), rowspan=2, colspan = 2)
    ax2 = plt.subplot2grid((3,4), (0,2), rowspan=2, colspan = 2)
    ax3 = plt.subplot2grid((3,4), (2,0), colspan = 2)
    ax4 = plt.subplot2grid((3,4), (2,2), colspan = 2)

    fig.subplots_adjust(left=0.07,bottom=0.13,right=0.97,top=0.9, wspace=0.44, hspace=0.13) 
    
    ax1.set_title('Major axis', fontsize=12,color = 'grey', weight = 'bold')
    ax1.plot(xx, major_prf, 'k', linewidth = 1, label = 'Beam %s'%nbeam)
    ax1.plot(xx, fit_2D[c_max[0][0], c_max[1][0]-200:c_max[1][0]+200], '--', label = 'Gauss')
    ax1.plot(xx, fit_mom_x[c_max[1][0]-200:c_max[1][0]+200], '--', label = 'Moment')
    ax1.set_xlim([xx.min(),xx.max()])
    ax1.set_ylabel('Intensity [count]')
    ax1.legend()
    ax1.grid()

    ax2.set_title('Minor axis', fontsize=12,color = 'grey', weight = 'bold')
    ax2.plot(yy,minor_prf, 'k', linewidth = 1, label = 'Beam %s'%nbeam)
    ax2.plot(yy,fit_2D[cr-10:cr+10, cc], '--', label = 'Gauss')
    ax2.plot(yy,fit_mom_y[cr-10:cr+10], '--', label = 'Moment')
    ax2.legend()
    ax2.grid()
    
    ax3.plot(xx, major_res/np.std(major_prf), '.')
    ax3.plot(xx, major_res_mom/np.std(major_prf), '.')
    ax3.set_ylim([-1,1])
    ax3.set_xlim([xx.min(),xx.max()])
    ax3.set_xlabel('Size [pixel]')
    ax3.set_ylabel('Residual [$\sigma$]')
    ax3.grid()
    
    ax4.plot(yy, minor_res/np.std(minor_prf), '.')
    ax4.plot(yy, minor_res_mom/np.std(minor_prf), '.')
    
    ax4.set_ylim([-1,1])
    ax4.grid()
    plt.tight_layout()
    
    
#-------------------------------------------------
#  Display beam, gaussian model and moment model
#-------------------------------------------------

p = 1 #Power law for the colorbar normalization (i.e.: .5 for sqrt display)

vmax = np.max(fit_2D)

# Set up figure and image grid
fig = plt.figure(figsize=(9.5, 3))
fig.subplots_adjust(right=0.95, left = 0.05)
grid = ImageGrid(fig, 111,          # as in plt.subplot(111)
                 nrows_ncols=(1,4),
                 axes_pad=0.1,
                 share_all=False,
                 add_all  =True,
                 cbar_location="right",
                 cbar_mode="single",
                 cbar_size="5%",
                 cbar_pad=0.15,
                 )

im2 = grid[0].imshow(raw_im, interpolation = 'nearest', cmap = 'gist_stern', aspect = aspect, norm = PowerNorm(p))#, vmin = 0, vmax = vmax)
grid[0].axis([cc-200, cc+200, cr-10,cr+10])
grid[0].set_title('Raw data')
grid[1].imshow(beam, interpolation = 'nearest', cmap = 'gist_stern', aspect = aspect, vmin = 0, vmax = vmax, norm = PowerNorm(p))
grid[1].axis([cc-200, cc+200, cr-10,cr+10])
grid[1].set_title('Raw data - Dark')
grid[2].imshow(fit_2D, interpolation = 'nearest', cmap = 'gist_stern', aspect = aspect, vmin = 0, vmax = vmax, norm = PowerNorm(p))
grid[2].contour(fit_2D, [.5*np.max(fit_2D)], colors = 'orange', linewidths=1, linestyles = '--')
grid[2].axis([cc-200, cc+200, cr-10,cr+10])
grid[2].set_title('Gaussian model')
grid[3].imshow(beam_thres, interpolation = 'nearest', cmap = 'gist_stern', aspect = aspect, vmin = 0, vmax = vmax, norm = PowerNorm(p))
grid[3].contour(fit_mom_2D, [.5*np.max(fit_mom_2D)], colors = 'orange', linewidths=1, linestyles = '--')
grid[3].axis([cc-200, cc+200, cr-10,cr+10])
grid[3].set_title('Thresholded image')
grid[3].cax.colorbar(im2)
grid[3].cax.toggle_label(True)
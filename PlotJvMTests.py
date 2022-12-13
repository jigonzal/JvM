from matplotlib import rcParams
rcParams.update({'figure.autolayout': True})
import matplotlib.pyplot as plt
from astropy.coordinates import ICRS
from astropy import units as u
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from photutils import CircularAperture
from photutils import SkyCircularAperture
from photutils import aperture_photometry
from astropy.modeling import models, fitting
import seaborn as sns
from astropy.io import fits
import numpy as np
sns.set_style("white", {'legend.frameon': True})
sns.set_style("ticks", {'legend.frameon': True})
sns.set_context("talk")
sns.set_palette('Dark2',desat=1)
cc = sns.color_palette()
from spectral_cube import SpectralCube
from astropy.convolution import convolve,convolve_fft,Kernel2D
from astropy.modeling.models import Gaussian2D
from essence import *

def GetBeam(hdulist):
    head = hdulist[0].header
    pix_size = np.abs(head['CDELT2'])
    return head['BMAJ']/pix_size,head['BMIN']/pix_size,head['BPA']

def GetKernel(CubePath):
    hdulist =   fits.open(CubePath,memmap=True)
    head = hdulist[0].header
    data = hdulist[0].data[0]

    try:
        BMAJ = hdulist[1].data.field('BMAJ')
        BMIN = hdulist[1].data.field('BMIN')
        BPA = hdulist[1].data.field('BPA')
    except:
        BMAJ = []
        BMIN = []
        BPA = []
        for i in range(len(data)):
            BMAJ.append(head['BMAJ']*3600.0)
            BMIN.append(head['BMIN']*3600.0)
            BPA.append(head['BPA'])
        BMAJ = np.array(BMAJ)
        BMIN = np.array(BMIN)
        BPA = np.array(BPA)

    pix_size = head['CDELT2']*3600.0
    factor = 2*(np.pi*BMAJ*BMIN/(8.0*np.log(2)))/(pix_size**2)
    factor = 1.0/factor
    FractionBeam = 1.0/np.sqrt(2.0)
    FractionBeam = 1.0
    KernelList = []

    for i in range(len(BMAJ)):
        SigmaPixel = int((BMAJ[i]*FractionBeam/2.355)/pix_size)+1
        x = np.arange(-(3*SigmaPixel), (3*SigmaPixel))
        y = np.arange(-(3*SigmaPixel), (3*SigmaPixel))        
        x, y = np.meshgrid(x, y)
        arr = models.Gaussian2D(amplitude=1.0,x_mean=0,y_mean=0,
                                x_stddev=(BMAJ[i]*FractionBeam/2.355)/pix_size,
                                y_stddev=(BMIN[i]*FractionBeam/2.355)/pix_size,
                                theta=(BPA[i]*2.0*np.pi/360.0)+np.pi/2)(x,y)
        kernel = Kernel2D(model=models.Gaussian2D(amplitude=1.0,x_mean=0,y_mean=0,
                            x_stddev=(BMAJ[i]*FractionBeam/2.355)/pix_size,
                            y_stddev=(BMIN[i]*FractionBeam/2.355)/pix_size,
                            theta=(BPA[i]*2.0*np.pi/360.0)+np.pi/2),
                            array=arr,width=len(x))
        KernelList.append(kernel)

    return KernelList[0],pix_size,1.0/factor

#Epsilon: 0.1435

CleanImageShallowPath = 'Point_Imaging_ShallowClean.fits'
CleanImageShallowJvMPath = 'Point_Imaging_ShallowClean_JvM.fits'
CleanImageDeepPath = 'Point_Imaging_DeepClean.fits'
CleanImageDeepJvMPath = 'Point_Imaging_DeepClean_JvM.fits'



CleanImageShallowImage = fits.open(CleanImageShallowPath)[0].data[0][0]
CleanImageShallowJvMImage = fits.open(CleanImageShallowJvMPath)[0].data[0][0]
CleanImageDeepImage = fits.open(CleanImageDeepPath)[0].data[0][0]
CleanImageDeepJvMImage = fits.open(CleanImageDeepJvMPath)[0].data[0][0]


kernel,pix_size,factor = GetKernel(CleanImageShallowPath)


wcs = WCS(fits.open(CleanImageShallowPath)[0].header)
wcs2 = wcs.celestial

position = SkyCoord(ra='00:00:00.0',dec='-23:00:00.0', unit=(u.hourangle, u.deg))

Radii = np.arange(0.02,1.02,0.02)

FluxCleanShallow = []
FluxCleanShallowJvM = []
FluxCleanDeep = []
FluxCleanDeepJvM = []

for i in Radii:
    aperture = SkyCircularAperture(position, r=i * u.arcsec)
    pix_aperture = aperture.to_pixel(wcs2)

    phot_table = aperture_photometry(CleanImageShallowImage, pix_aperture)
    FluxCleanShallow.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageShallowJvMImage, pix_aperture)
    FluxCleanShallowJvM.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageDeepImage, pix_aperture)
    FluxCleanDeep.append(phot_table['aperture_sum']/factor[0])

    phot_table = aperture_photometry(CleanImageDeepJvMImage, pix_aperture)
    FluxCleanDeepJvM.append(phot_table['aperture_sum']/factor[0])    


FluxCleanShallow = np.array(FluxCleanShallow)
FluxCleanShallowJvM = np.array(FluxCleanShallowJvM)
FluxCleanDeep = np.array(FluxCleanDeep)
FluxCleanDeepJvM = np.array(FluxCleanDeepJvM)


w, h = 1.5*plt.figaspect(0.9)
fig1 = plt.figure(figsize=(w,h))
plt.plot(Radii,1e3*FluxCleanShallow,label='Shallow clean')
plt.plot(Radii,1e3*FluxCleanShallowJvM,label='Shallow clean + JvM')
plt.plot(Radii,1e3*FluxCleanDeep,label='Deep clean')
plt.plot(Radii,1e3*FluxCleanDeepJvM,label='Deep clean + JvM')
plt.axhline(0.99,ls='--',color='gray',label='Flux density in uv-plane')
plt.axvline(0.13,ls=':',color='gray',label='(Bmaj+Bmin)/2')

plt.legend(loc=0)
plt.fill_between(Radii,(0.99-5e-02),(0.99+5e-02),alpha=0.5,color='gray',label='Flux density in uv-plane')
plt.xlabel('Radius [Arcsec]')
plt.ylabel(r'Flux density [$\mathrm{mJy}$]')
plt.ylim(0.5,1.9)
plt.savefig('FluxMeasurementsWithJvM_Point.pdf')

if True:
    Size = 12.0

    kernel,pix_size,factor = GetKernel(CleanImageShallowPath)
    ShallowCleanCube = SpectralCube.read(CleanImageShallowPath)[0]
    OriginalCube_pb = SpectralCube.read('Disk_Imaging_Dirty_pb.fits')
    NoiseShallowClean = mk_noisemap(ShallowCleanCube,OriginalCube_pb,pbfrac=0.2,sigma=6) 

    kernel,pix_size,factor = GetKernel(CleanImageDeepPath)
    DeepCleanCube = SpectralCube.read(CleanImageDeepPath)[0]
    NoiseDeepClean = mk_noisemap(DeepCleanCube,OriginalCube_pb,pbfrac=0.2,sigma=6) 


    kernel,pix_size,factor = GetKernel(CleanImageShallowJvMPath)
    ShallowJvMCube = SpectralCube.read(CleanImageShallowJvMPath)[0]
    NoiseShallowJvM = mk_noisemap(ShallowJvMCube,OriginalCube_pb,pbfrac=0.2,sigma=6) 

    RMSShallow = np.nanstd(NoiseShallowClean)
    RMSDeep = np.nanstd(NoiseDeepClean)
    RMSShallowJvM = np.nanstd(NoiseShallowJvM)

    SigmaShallow = []
    SigmaShallowJvM = []
    SigmaDeep = []
    for index,i in enumerate( Radii):
        FluxDensityShallow = []
        FluxDensityShallow = np.array(FluxDensityShallow)
     

        FluxDensityShallowJvM = []
        FluxDensityShallowJvM = np.array(FluxDensityShallowJvM)

        FluxDensityDeep = []
        FluxDensityDeep = np.array(FluxDensityDeep)

        for j in np.arange(-1.0*Size,Size,2*i):
            for k in np.arange(-1.0*Size,Size,2*i):
                NewPosition = SkyCoord(ra=position.ra.deg+j*np.cos(position.dec.radian)/3600.0,dec=position.dec.deg+k/3600.0, unit=(u.deg, u.deg))
                aperture = SkyCircularAperture(NewPosition, r=i * u.arcsec)
                pix_aperture = aperture.to_pixel(wcs2)

                phot_table = aperture_photometry(NoiseShallowClean, pix_aperture)

                FluxDensityShallow = np.append(FluxDensityShallow, phot_table['aperture_sum'][0]/factor[0])
                
                phot_table = aperture_photometry(NoiseShallowJvM, pix_aperture)

                FluxDensityShallowJvM = np.append(FluxDensityShallowJvM, phot_table['aperture_sum'][0]/factor[0])    

                phot_table = aperture_photometry(NoiseDeepClean, pix_aperture)

                FluxDensityDeep = np.append(FluxDensityDeep, phot_table['aperture_sum'][0]/factor[0])                     

                if len(FluxDensityShallow[np.isfinite(FluxDensityShallow)])>200:
                    break


        SigmaShallow.append(np.nanstd(FluxDensityShallow))
        SigmaShallowJvM.append(np.nanstd(FluxDensityShallowJvM))
        SigmaDeep.append(np.nanstd(FluxDensityDeep))



    w, h = 1.5*plt.figaspect(0.9)
    fig1 = plt.figure(figsize=(w,h))

    SigmaShallowJvM = np.array(SigmaShallowJvM)
    SigmaDeep = np.array(SigmaDeep)

    plt.plot(Radii,np.array(SigmaShallowJvM)*1e6,label='Deep clean + JvM Error')
    plt.plot(Radii,np.array(SigmaDeep)*1e6,label='Deep clean Error')
    plt.semilogy()
    plt.axhline(1e6*5e-05,ls='--',color=cc[2],label='Error Point source in uv-plane')
    plt.axhline(0.155*1000.0,ls='--',color=cc[3],label='Error Disk source in uv-plane')
    plt.axvline(0.5,ls=':',color='red',label='Disk size')
    plt.axvline(0.13,ls=':',color='gray',label='(Bmaj+Bmin)/2')
    plt.legend(loc=0)
    plt.xlabel('Radius [Arcsec]')
    plt.ylabel(r'Error [$\mu \mathrm{Jy}$]')
    plt.savefig('ErrorMeasurementsWithJvM_Point.pdf')




############################################################
############################################################
############################################################


CleanImageShallowPath = 'Disk_Imaging_ShallowClean.fits'
CleanImageShallowJvMPath = 'Disk_Imaging_ShallowClean_JvM.fits'
CleanImageDeepPath = 'Disk_Imaging_DeepClean.fits'
CleanImageDeepJvMPath = 'Disk_Imaging_DeepClean_JvM.fits'




CleanImageShallowImage = fits.open(CleanImageShallowPath)[0].data[0][0]
CleanImageShallowJvMImage = fits.open(CleanImageShallowJvMPath)[0].data[0][0]
CleanImageDeepImage = fits.open(CleanImageDeepPath)[0].data[0][0]
CleanImageDeepJvMImage = fits.open(CleanImageDeepJvMPath)[0].data[0][0]


kernel,pix_size,factor = GetKernel(CleanImageShallowPath)


wcs = WCS(fits.open(CleanImageShallowPath)[0].header)
wcs2 = wcs.celestial

position = SkyCoord(ra='00:00:00.0',dec='-23:00:00.0', unit=(u.hourangle, u.deg))


FluxCleanShallow = []
FluxCleanShallowJvM = []
FluxCleanDeep = []
FluxCleanDeepJvM = []

for i in Radii:
    aperture = SkyCircularAperture(position, r=i * u.arcsec)
    pix_aperture = aperture.to_pixel(wcs2)

    phot_table = aperture_photometry(CleanImageShallowImage, pix_aperture)
    FluxCleanShallow.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageShallowJvMImage, pix_aperture)
    FluxCleanShallowJvM.append(phot_table['aperture_sum'][0]/factor[0])
       

    phot_table = aperture_photometry(CleanImageDeepImage, pix_aperture)
    FluxCleanDeep.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageDeepJvMImage, pix_aperture)
    FluxCleanDeepJvM.append(phot_table['aperture_sum'][0]/factor[0])    

   


FluxCleanShallow = np.array(FluxCleanShallow)
FluxCleanShallowJvM = np.array(FluxCleanShallowJvM)
FluxCleanDeep = np.array(FluxCleanDeep)
FluxCleanDeepJvM = np.array(FluxCleanDeepJvM)


w, h = 1.5*plt.figaspect(0.9)
fig1 = plt.figure(figsize=(w,h))
plt.plot(Radii,1e3*FluxCleanShallow,label='Shallow clean')
plt.plot(Radii,1e3*FluxCleanShallowJvM,label='Shallow clean + JvM')
plt.plot(Radii,1e3*FluxCleanDeep,label='Deep Clean')
plt.plot(Radii,1e3*FluxCleanDeepJvM,label='Deep Clean + JvM')

ModelFlux = 10*np.power(Radii,2)/np.power(0.5,2)
ModelFlux[Radii>=0.5] = 10.0
plt.plot(Radii,ModelFlux,label='Model')
plt.axhline(9.85,ls='--',color='gray',label='Flux density in uv-plane')
plt.axvline(0.5,ls=':',color='red',label='Disk size')
plt.legend(loc=0)
plt.fill_between(Radii,(9.85-0.27),(9.85+0.27),alpha=0.5,color='gray',label='Flux density in uv-plane')

plt.xlabel('Radius [Arcsec]')
plt.ylabel(r'Flux density [$\mathrm{mJy}$]')
plt.ylim(5,13)
plt.savefig('FluxMeasurementsWithJvM_Disk.pdf')

plt.close('all')
w, h = 1.5*plt.figaspect(0.9)
fig1 = plt.figure(figsize=(w,h))

Area = np.pi*np.power(Radii,2)

plt.plot(Radii,1e3*FluxCleanShallow/Area,label='Shallow clean',color=cc[0])
plt.plot(Radii,1e3*FluxCleanShallowJvM/Area,label='Shallow clean + JvM',color=cc[1])
plt.plot(Radii,1e3*FluxCleanDeep/Area,label='Deep clean',color=cc[2])
plt.plot(Radii,1e3*FluxCleanDeepJvM/Area,label='Deep clean + JvM',color=cc[3])

plt.fill_between(Radii,1e3*FluxCleanShallow/Area - 1e3*SigmaDeep/Area,1e3*FluxCleanShallow/Area + 1e3*SigmaDeep/Area,alpha=0.5,color=cc[0])
plt.fill_between(Radii,1e3*FluxCleanShallowJvM/Area - 1e3*SigmaShallowJvM/Area,1e3*FluxCleanShallowJvM/Area + 1e3*SigmaShallowJvM/Area,alpha=0.5,color=cc[1])
plt.fill_between(Radii,1e3*FluxCleanDeep/Area - 1e3*SigmaDeep/Area,1e3*FluxCleanDeep/Area + 1e3*SigmaDeep/Area,alpha=0.5,color=cc[2])
plt.fill_between(Radii,1e3*FluxCleanDeepJvM/Area - 1e3*SigmaShallowJvM/Area,1e3*FluxCleanDeepJvM/Area + 1e3*SigmaShallowJvM/Area,alpha=0.5,color=cc[3])

plt.axhline(9.8/(np.power(0.5,2)*np.pi),label='Model',ls='--',color='gray')
plt.axvline(0.5,ls=':',color='red',label='Disk size')
plt.legend(loc=0)
plt.xlabel('Radius [Arcsec]')
plt.ylabel(r'Surface Flux density [$\mathrm{mJy}/arcsec^2$]')
plt.savefig('SurfaceFluxMeasurementsWithJvM_Disk.pdf')



############################################################
############################################################
############################################################


CleanImageShallowPath = 'Disk_Imaging_ShallowClean_MS.fits'
CleanImageShallowJvMPath = 'Disk_Imaging_ShallowClean_MS_JvM.fits'
CleanImageDeepPath = 'Disk_Imaging_DeepClean_MS.fits'
CleanImageDeepJvMPath = 'Disk_Imaging_DeepClean_MS_JvM.fits'



CleanImageShallowImage = fits.open(CleanImageShallowPath)[0].data[0][0]
CleanImageShallowJvMImage = fits.open(CleanImageShallowJvMPath)[0].data[0][0]
CleanImageDeepImage = fits.open(CleanImageDeepPath)[0].data[0][0]
CleanImageDeepJvMImage = fits.open(CleanImageDeepJvMPath)[0].data[0][0]


kernel,pix_size,factor = GetKernel(CleanImageShallowPath)


wcs = WCS(fits.open(CleanImageShallowPath)[0].header)
wcs2 = wcs.celestial

position = SkyCoord(ra='00:00:00.0',dec='-23:00:00.0', unit=(u.hourangle, u.deg))

FluxCleanShallow = []
FluxCleanShallowJvM = []
FluxCleanDeep = []
FluxCleanDeepJvM = []

for i in Radii:
    aperture = SkyCircularAperture(position, r=i * u.arcsec)
    pix_aperture = aperture.to_pixel(wcs2)

    phot_table = aperture_photometry(CleanImageShallowImage, pix_aperture)
    FluxCleanShallow.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageShallowJvMImage, pix_aperture)
    FluxCleanShallowJvM.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageDeepImage, pix_aperture)
    FluxCleanDeep.append(phot_table['aperture_sum'][0]/factor[0])

    phot_table = aperture_photometry(CleanImageDeepJvMImage, pix_aperture)
    FluxCleanDeepJvM.append(phot_table['aperture_sum'][0]/factor[0])    


FluxCleanShallow = np.array(FluxCleanShallow)
FluxCleanShallowJvM = np.array(FluxCleanShallowJvM)
FluxCleanDeep = np.array(FluxCleanDeep)
FluxCleanDeepJvM = np.array(FluxCleanDeepJvM)


w, h = 1.5*plt.figaspect(0.9)
fig1 = plt.figure(figsize=(w,h))
plt.plot(Radii,1e3*FluxCleanShallow,label='Shallow clean')
plt.plot(Radii,1e3*FluxCleanShallowJvM,label='Shallow clean + JvM')
plt.plot(Radii,1e3*FluxCleanDeep,label='Deep clean')
plt.plot(Radii,1e3*FluxCleanDeepJvM,label='Deep clean + JvM')

ModelFlux = 10*np.power(Radii,2)/np.power(0.5,2)
ModelFlux[Radii>=0.5] = 10.0
plt.plot(Radii,ModelFlux,label='Model')
plt.axhline(9.85,ls='--',color='gray',label='Flux density in uv-plane')
plt.axvline(0.5,ls=':',color='red',label='Disk size')

plt.legend(loc=0)
plt.fill_between(Radii,(9.85-0.27),(9.85+0.27),alpha=0.5,color='gray',label='Flux density in uv-plane')
plt.xlabel('Radius [Arcsec]')
plt.ylabel(r'Flux density [$\mathrm{mJy}$]')
# plt.semilogy()
plt.ylim(5,13)
plt.savefig('FluxMeasurementsWithJvM_Disk_MS.pdf')

plt.close('all')
w, h = 1.5*plt.figaspect(0.9)
fig1 = plt.figure(figsize=(w,h))

Area = np.pi*np.power(Radii,2)

plt.plot(Radii,1e3*FluxCleanShallow/Area,label='Shallow clean',color=cc[0])
plt.plot(Radii,1e3*FluxCleanShallowJvM/Area,label='Shallow clean + JvM',color=cc[1])
plt.plot(Radii,1e3*FluxCleanDeep/Area,label='Deep clean',color=cc[2])
plt.plot(Radii,1e3*FluxCleanDeepJvM/Area,label='Deep clean + JvM',color=cc[3])

plt.fill_between(Radii,1e3*FluxCleanShallow/Area - 1e3*SigmaDeep/Area,1e3*FluxCleanShallow/Area + 1e3*SigmaDeep/Area,alpha=0.5,color=cc[0])
plt.fill_between(Radii,1e3*FluxCleanShallowJvM/Area - 1e3*SigmaShallowJvM/Area,1e3*FluxCleanShallowJvM/Area + 1e3*SigmaShallowJvM/Area,alpha=0.5,color=cc[1])
plt.fill_between(Radii,1e3*FluxCleanDeep/Area - 1e3*SigmaDeep/Area,1e3*FluxCleanDeep/Area + 1e3*SigmaDeep/Area,alpha=0.5,color=cc[2])
plt.fill_between(Radii,1e3*FluxCleanDeepJvM/Area - 1e3*SigmaShallowJvM/Area,1e3*FluxCleanDeepJvM/Area + 1e3*SigmaShallowJvM/Area,alpha=0.5,color=cc[3])

plt.axhline(9.8/(np.power(0.5,2)*np.pi),label='Model',ls='--',color='gray')
plt.axvline(0.5,ls=':',color='red',label='Disk size')

plt.legend(loc=0)
plt.xlabel('Radius [Arcsec]')
plt.ylabel(r'Surface Flux density [$\mathrm{mJy}/arcsec^2$]')
plt.savefig('SurfaceFluxMeasurementsWithJvM_Disk_MS.pdf')

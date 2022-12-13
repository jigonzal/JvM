import numpy as np
import os
import sys 
sys.path.append("analysis_scripts/")
import analysisUtils as au 
import matplotlib.pyplot as plt
import numpy as np

def ApplyCorrection(ImagePath,MaximumRadius):

	BeamParameters = au.getFitsBeam(ImagePath+'.image')
	os.system('rm -rf CleanBeamImage.image')
	au.makeGaussianForImage(ImagePath+'.image',1,BeamParameters[5]/2,BeamParameters[6]/2,BeamParameters[0]/np.abs(BeamParameters[3]),BeamParameters[1]/np.abs(BeamParameters[3]),BeamParameters[2],'CleanBeamImage.image')

	RangePixels = np.arange(1,MaximumRadius/np.abs(BeamParameters[3]),1)

	CleanBeamVolume = []
	DirtyBeamVolume = []
	for r in RangePixels:
		myoutput = imstat('CleanBeamImage.image',region='circle[['+str(BeamParameters[5]/2)+'pix, '+str(BeamParameters[6]/2)+'pix], '+str(r)+'pix]')
		# print(r,myoutput['sum'])
		CleanBeamVolume.append(myoutput['sum'])

	for r in RangePixels:
		myoutput = imstat(ImagePath+'.psf',region='circle[['+str(BeamParameters[5]/2)+'pix, '+str(BeamParameters[6]/2)+'pix], '+str(r)+'pix]')
		# print(r,myoutput['sum'])	
		DirtyBeamVolume.append(myoutput['sum'])

	CleanBeamVolume = np.array(CleanBeamVolume)
	DirtyBeamVolume = np.array(DirtyBeamVolume)

	epsilon = np.min(CleanBeamVolume/DirtyBeamVolume)
	print('Epsilon:',epsilon)

	plt.figure()
	plt.plot(RangePixels*np.abs(BeamParameters[3]),CleanBeamVolume/DirtyBeamVolume,lw=2)
	plt.xlabel('Radius Arcsec')
	plt.ylabel('Ratio Clean/Dirty Beam Volume')
	plt.axhline(epsilon,ls='--',color='black',lw=2)
	plt.savefig('RatioVolumesBeams.pdf')

	plt.figure()
	plt.plot(RangePixels*np.abs(BeamParameters[3]),CleanBeamVolume,ls='--',label='Clean Beam',lw=2)
	plt.plot(RangePixels*np.abs(BeamParameters[3]),DirtyBeamVolume,ls='-',label='Dirty Beam',lw=2)
	plt.xlabel('Radius Arcsec')
	plt.ylabel('Beam Volume')
	plt.legend(loc=0)
	plt.axvline(np.argmin(CleanBeamVolume/DirtyBeamVolume)*np.abs(BeamParameters[3]),ls='--',lw=2,color='black')
	plt.savefig('VolumeBeams.pdf')
	plt.close('all')


	immath(imagename=[ImagePath+'.image',ImagePath+'.residual'], expr='IM0-IM1',outfile=ImagePath+'.modelConvolved',imagemd=ImagePath+'.image')

	# imsmooth(imagename=ImagePath+'.model', kernel='gauss',major=BeamParameters[0], minor=BeamParameters[1], pa=BeamParameters[2],outfile=ImagePath+'.modelConvolved',overwrite=True)


	os.system('rm -rf '+ImagePath+'.imageJvM')
	immath(imagename=[ImagePath+'.residual',ImagePath+'.modelConvolved'], expr='IM0*'+str(epsilon)+'+IM1',outfile=ImagePath+'.imageJvM',imagemd=ImagePath+'.image')
	exportfits(imagename=ImagePath+'.imageJvM',fitsimage=ImagePath+'_JvM.fits',overwrite=True)
	return 

TimeOnSource = 1.0 #hours
User_pwv = 0.913
T_sky = 39.538
Tau0 = 0.158


direction = "J2000 00h00m00.0s -23d00m00.0s"

os.system('rm -rf Point.cl')
cl.done()
cl.addcomponent(dir=direction, flux=1e-3, fluxunit='Jy', freq='350.0GHz', shape="point")
cl.rename('Point.cl')
cl.done()

Option = 'Point'
ConcatMS = Option+'_concatenated.ms'

simobserve(project=Option,
		complist = Option+'.cl',
		compwidth = '1GHz',
		inwidth='1GHz',
		totaltime=str(3600*TimeOnSource)+'s',
		antennalist='alma.cycle8.7.cfg',
		thermalnoise='tsys-atm',
		seed=int(np.random.uniform(0,100000)),
		incenter='350.0GHz',
		user_pwv = User_pwv,
		t_ground = 273,
		t_sky = T_sky,
		tau0 = Tau0,	
		overwrite=True,
		incell = '0.01arcsec')		

simobserve(project=Option,
		complist = Option+'.cl',
		compwidth = '1GHz',
		inwidth='1GHz',
		totaltime=str(3600*0.5*TimeOnSource)+'s',
		antennalist='alma.cycle8.3.cfg',
		thermalnoise='tsys-atm',
		seed=int(np.random.uniform(0,100000)),
		incenter='350.0GHz',
		user_pwv = User_pwv,
		t_ground = 273,
		t_sky = T_sky,
		tau0 = Tau0,	
		overwrite=True,
		incell = '0.01arcsec')	

os.system('rm -rf '+ConcatMS)

concat(vis=[Option+'/'+Option+'.alma.cycle8.7.noisy.ms',Option+'/'+Option+'.alma.cycle8.3.noisy.ms'],concatvis=ConcatMS)


ImageOutput = Option + '_Imaging_ShallowClean'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=1000,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=3.0,
	weighting='natural',
	robust=2.0,
	deconvolver='hogbom', 
	scales=[0,10,30,50],
	interactive=False,
	smallscalebias=0)

ImageOutput = Option + '_Imaging_DeepClean'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=1000,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=1.0,
	weighting='natural',
	robust=2.0,
	deconvolver='hogbom', 
	scales=[0,10,30,50],
	interactive=False,
	smallscalebias=0)


ImageOutput = Option + '_Imaging_Dirty'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=0,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=3.0,
	weighting='natural',
	robust=2.0,
	deconvolver='hogbom', 
	scales=[0,10,30,50],
	interactive=False,
	smallscalebias=0)



uvmodelfit(vis=ConcatMS,field='0',spw='0', niter=5, comptype='P', sourcepar=[0.001, 0.0, 0.0],varypar=[True,False,False], outfile='Fitted_UV_Point')

# There are 1950480 - 1 = 1950479 degrees of freedom.
#  iter=0:   reduced chi2=0.00244689:  I=0.001,  dir=[0, 0] arcsec
#  iter=1:   reduced chi2=0.00244689:  I=0.000998618,  dir=[0, 0] arcsec
#  iter=2:   reduced chi2=0.00244689:  I=0.000998617,  dir=[0, 0] arcsec
#  iter=3:   reduced chi2=0.00244689:  I=0.000998617,  dir=[0, 0] arcsec
#  iter=4:   reduced chi2=0.00244689:  I=0.000998617,  dir=[0, 0] arcsec
#  iter=5:   reduced chi2=0.00244689:  I=0.000998617,  dir=[0, 0] arcsec
# If data weights are arbitrarily scaled, the following formal errors
#  will be underestimated by at least a factor sqrt(reduced chi2). If 
#  the fit is systematically poor, the errors are much worse.
# I = 0.000998617 +/- 0.00101261
# x = 0 +/- 0 arcsec
# y = 0 +/- 0 arcsec
# Writing componentlist to file: /Users/jgonzal/Desktop/Projects/JvM/Simulations/Fitted_UV_Point

ApplyCorrection(Option+'_Imaging_ShallowClean',2)
ApplyCorrection(Option+'_Imaging_DeepClean',2)
ApplyCorrection(Option+'_Imaging_Dirty',2)



exportfits(imagename=Option+'_Imaging_ShallowClean.image',fitsimage=Option+'_Imaging_ShallowClean.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_ShallowClean.mask',fitsimage=Option+'_Imaging_ShallowClean_mask.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_ShallowClean.residual',fitsimage=Option+'_Imaging_ShallowClean_residual.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean.image',fitsimage=Option+'_Imaging_DeepClean.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean.mask',fitsimage=Option+'_Imaging_DeepClean_mask.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean.residual',fitsimage=Option+'_Imaging_DeepClean_residual.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty.image',fitsimage=Option+'_Imaging_Dirty.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty.psf',fitsimage=Option+'_Imaging_Dirty_psf.fits',overwrite=True)



os.system('rm -rf Disk.cl')
cl.done()
cl.addcomponent(dir=direction, flux=10e-3, fluxunit='Jy', freq='350.0GHz', shape="disk",majoraxis='1arcsec',minoraxis='1arcsec',positionangle='0.0deg')
cl.rename('Disk.cl')
cl.done()

Option = 'Disk'
ConcatMS = Option+'_concatenated.ms'

simobserve(project=Option,
		complist = Option+'.cl',
		compwidth = '1GHz',
		inwidth='1GHz',
		totaltime=str(3600*TimeOnSource)+'s',
		antennalist='alma.cycle8.7.cfg',
		thermalnoise='tsys-atm',
		seed=int(np.random.uniform(0,100000)),
		incenter='350.0GHz',
		user_pwv = User_pwv,
		t_ground = 273,
		t_sky = T_sky,
		tau0 = Tau0,	
		overwrite=True,
		incell = '0.01arcsec')	

simobserve(project=Option,
		complist = Option+'.cl',
		compwidth = '1GHz',
		inwidth='1GHz',
		totaltime=str(3600*0.5*TimeOnSource)+'s',
		antennalist='alma.cycle8.3.cfg',
		thermalnoise='tsys-atm',
		seed=int(np.random.uniform(0,100000)),
		incenter='350.0GHz',
		user_pwv = User_pwv,
		t_ground = 273,
		t_sky = T_sky,
		tau0 = Tau0,	
		overwrite=True,
		incell = '0.01arcsec')	

os.system('rm -rf '+ConcatMS)

concat(vis=[Option+'/'+Option+'.alma.cycle8.7.noisy.ms',Option+'/'+Option+'.alma.cycle8.3.noisy.ms'],concatvis=ConcatMS)
uvmodelfit(vis=ConcatMS,field='0',spw='0', niter=5, comptype='D', sourcepar=[0.01, 0.0, 0.0,1,0.9,0],varypar=[True,False,False,True,True,True], outfile='Fitted_UV_Disk')
uvmodelfit(vis=ConcatMS,field='0',spw='0', niter=5, comptype='D', sourcepar=[0.01, 0.0, 0.0,1,1,0],varypar=[True,False,False,False,False,False], outfile='Fitted_UV_Disk2')

# There are 1950480 - 4 = 1950476 degrees of freedom.
#  iter=0:   reduced chi2=0.00244056:  I=0.01,  dir=[0, 0] arcsec,  shape=[1, 0.9, 0]
#  iter=1:   reduced chi2=0.00244054:  I=0.010084,  dir=[0, 0] arcsec,  shape=[1.03958, 1, 34.2886]
#  iter=2:   reduced chi2=0.00244052:  I=0.00986755,  dir=[0, 0] arcsec,  shape=[1.02817, 0.928359, 34.2886]
#  iter=3:   reduced chi2=0.00244052:  I=0.00985111,  dir=[0, 0] arcsec,  shape=[1.03794, 0.945721, 46.0772]
#  iter=4:   reduced chi2=0.00244052:  I=0.00985111,  dir=[0, 0] arcsec,  shape=[1.03794, 0.945721, 46.0772]
#  iter=5:   reduced chi2=0.00244052:  I=0.00985111,  dir=[0, 0] arcsec,  shape=[1.03794, 0.945721, 46.0772]
# If data weights are arbitrarily scaled, the following formal errors
#  will be underestimated by at least a factor sqrt(reduced chi2). If 
#  the fit is systematically poor, the errors are much worse.
# I = 0.00985111 +/- 0.0055852
# x = 0 +/- 0 arcsec
# y = 0 +/- 0 arcsec
# a = 1.03794 +/- 1.58048 arcsec
# r = 0.945721 +/- 1.83528
# p = 46.0772 +/- 973.001 deg
# Writing componentlist to file: /Users/jgonzal/Desktop/Projects/JvM/Simulations/Fitted_UV_Disk


# There are 1950480 - 1 = 1950479 degrees of freedom.
#  iter=0:   reduced chi2=0.00244052:  I=0.01,  dir=[0, 0] arcsec,  shape=[1, 1, 0]
#  iter=1:   reduced chi2=0.00244052:  I=0.00979313,  dir=[0, 0] arcsec,  shape=[1, 1, 0]
#  iter=2:   reduced chi2=0.00244052:  I=0.00979292,  dir=[0, 0] arcsec,  shape=[1, 1, 0]
#  iter=3:   reduced chi2=0.00244052:  I=0.00979292,  dir=[0, 0] arcsec,  shape=[1, 1, 0]
#  iter=4:   reduced chi2=0.00244052:  I=0.00979292,  dir=[0, 0] arcsec,  shape=[1, 1, 0]
#  iter=5:   reduced chi2=0.00244052:  I=0.00979292,  dir=[0, 0] arcsec,  shape=[1, 1, 0]
# If data weights are arbitrarily scaled, the following formal errors
#  will be underestimated by at least a factor sqrt(reduced chi2). If 
#  the fit is systematically poor, the errors are much worse.
# I = 0.00979292 +/- 0.0031469
# x = 0 +/- 0 arcsec
# y = 0 +/- 0 arcsec
# a = 1 +/- 0 arcsec
# r = 1 +/- 0
# p = 0 +/- 0 deg
# Writing componentlist to file: /Users/jgonzal/Desktop/Projects/JvM/Simulations/Fitted_UV_Disk2

#sigma disk 0.00015546195611080286 Jy


ImageOutput = Option + '_Imaging_ShallowClean'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=1000,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=3.0,
	weighting='natural',
	robust=2.0,
	deconvolver='hogbom', 
	scales=[0,10,30,50],
	interactive=False,
	smallscalebias=0)



ImageOutput = Option + '_Imaging_DeepClean'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=1000,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=1.0,
	weighting='natural',
	robust=2.0,
	deconvolver='hogbom', 
	scales=[0,10,30,50],
	interactive=False,
	smallscalebias=0)

ImageOutput = Option + '_Imaging_Dirty'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=0,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=3.0,
	weighting='natural',
	robust=2.0,
	deconvolver='hogbom', 
	scales=[0,10,30,50],
	interactive=False,
	smallscalebias=0)


ApplyCorrection(Option+'_Imaging_ShallowClean',2)
ApplyCorrection(Option+'_Imaging_DeepClean',2)
ApplyCorrection(Option+'_Imaging_Dirty',2)


#epsilon 0.1435

exportfits(imagename=Option+'_Imaging_ShallowClean.image',fitsimage=Option+'_Imaging_ShallowClean.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_ShallowClean.mask',fitsimage=Option+'_Imaging_ShallowClean_mask.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_ShallowClean.residual',fitsimage=Option+'_Imaging_ShallowClean_residual.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean.image',fitsimage=Option+'_Imaging_DeepClean.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean.mask',fitsimage=Option+'_Imaging_DeepClean_mask.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean.residual',fitsimage=Option+'_Imaging_DeepClean_residual.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty.image',fitsimage=Option+'_Imaging_Dirty.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty.psf',fitsimage=Option+'_Imaging_Dirty_psf.fits',overwrite=True)


############################


ImageOutput = Option + '_Imaging_ShallowClean_MS'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=1000,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=3.0,
	weighting='natural',
	robust=2.0,
	deconvolver='multiscale', 
	scales=[0,13,40],
	interactive=False,
	smallscalebias=0)


ImageOutput = Option + '_Imaging_DeepClean_MS'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=1000,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=1.0,
	weighting='natural',
	robust=2.0,
	deconvolver='multiscale', 
	scales=[0,13,40],
	interactive=False,
	smallscalebias=0)

ImageOutput = Option + '_Imaging_Dirty_MS'
os.system('rm -rf '+ImageOutput+'.*')

tclean(vis=ConcatMS,
	imagename=ImageOutput,
	field='0',
	imsize=2500,
	cell='0.01arcsec',
	niter=0,
	usemask='auto-multithresh',
	smoothfactor=0.5,
	sidelobethreshold=2,
	noisethreshold=5,
	lownoisethreshold=2,
	minbeamfrac=0,
	nsigma=3.0,
	weighting='natural',
	robust=2.0,
	deconvolver='multiscale', 
	scales=[0,13,40],
	interactive=False,
	smallscalebias=0)


ApplyCorrection(Option+'_Imaging_ShallowClean_MS',2)
ApplyCorrection(Option+'_Imaging_DeepClean_MS',2)
ApplyCorrection(Option+'_Imaging_Dirty_MS',2)



exportfits(imagename=Option+'_Imaging_ShallowClean_MS.image',fitsimage=Option+'_Imaging_ShallowClean_MS.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_ShallowClean_MS.residual',fitsimage=Option+'_Imaging_ShallowClean_MS_residual.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_ShallowClean_MS.mask',fitsimage=Option+'_Imaging_ShallowClean_MS_mask.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean_MS.image',fitsimage=Option+'_Imaging_DeepClean_MS.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean_MS.residual',fitsimage=Option+'_Imaging_DeepClean_MS_residual.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_DeepClean_MS.mask',fitsimage=Option+'_Imaging_DeepClean_MS_mask.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty_MS.image',fitsimage=Option+'_Imaging_Dirty_MS.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty_MS.psf',fitsimage=Option+'_Imaging_Dirty_MS_psf.fits',overwrite=True)
exportfits(imagename=Option+'_Imaging_Dirty_MS.pb',fitsimage=Option+'_Imaging_Dirty_MS_pb.fits',overwrite=True)

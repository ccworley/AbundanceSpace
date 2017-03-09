#Plot Abundance data
from load_data import *
from math import *
import numpy as np

import matplotlib.pyplot as plt

lab = len(idr4abuncol)

#Extract GE_MW stars  (also look at GE_SD_GGC and GE_MW_BL
imat = np.where(idr4gestyp == 'GE_MW')

plt.close() 
for i in range(0,lab):
#for i in range(0,0):

	imat = np.where(idr4gesfld == gesfldlist[i])
	#imat = np.where(idr4gesfld == 'NGC1851')
	
	#For 
	teff = idr4teff[imat]
	logg = idr4logg[imat]
	feh = idr4feh[imat]
	vel = idr4vel[imat]
	na = idr4na1[imat]
	o = idr4o1[imat]
	ba = idr4ba2[imat]

	#Finds Vel and FEH limits in file
	ilims = np.where(lim == gesfldlist[i])
	#ilims = np.where(lim == 'NGC1851')

	#Extract indices of member  stars
	ilims = (fehllim[ilims] < feh) & (fehulim[ilims] > feh) & (vradllim[ilims] < vel) & (vradulim[ilims] > vel)

	mnfeh = np.mean(feh[ilims])
	sdfeh = np.std(feh[ilims])
#, r'$\mu=100,\ \sigma=15$'
	plt.ion()
	plt.scatter(o[ilims],na[ilims],c=ba[ilims], cmap='nipy_spectral')
	plt.xlabel('[O/Fe]')
	plt.ylabel('[Na/Fe]')
	#plt.clabel('[Ba/Fe]')
	plt.text(6500,0.5,'<[Fe/H]>=%5.2f'%mnfeh+'+/-%4.2f'%sdfeh)  # +'+/- %5.2f' sdfeh)
	#plt.axis([8000,3000,5,0])
	plt.title(gesfldlist[i])
	plt.colorbar()
	#plt.gca().invert_xaxis()
	#plt.gca().invert_yaxis()
	plt.show()
	plt.savefig('Figures/'+gesfldlist[i]+'_NaOBa.pdf')
	plt.close() 

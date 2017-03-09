#Plot Abundance data
from load_data import *
from math import *
import numpy as np
from astropy.table import Table, Column

import matplotlib.pyplot as plt
label_size = 6
plt.rcParams['xtick.labelsize'] = label_size 
plt.rcParams['ytick.labelsize'] = label_size 

lab = len(idr4abuncol)

#Extract GE_MW stars  (also look at GE_SD_GGC and GE_MW_BL
blin = (idr4gesfld == 'NGC1851') & (idr4recwg == 'WG11')
#blin = (idr4gestyp == 'GE_MW') & (idr4recwg == 'WG11')
imat = np.where(blin==True)
#imat = imat1[imat2]
spgf = idr4gestyp[imat]
lspe = len(spgf)

#Find where abundancecolumns have values
abcolkp = []
for i in range(0,lab):
#for i in range(0,0):

	abvect = idr4data.field(idr4abuncol[i])
	mwabvect = abvect[imat]
	inan = np.where(np.isfinite(mwabvect)==False)
	ireal = np.where(np.isfinite(mwabvect)==True)
	#imat = np.where(idr4gesfld == 'NGC1851')

	lnan = len(mwabvect[inan])
	lvect = len(mwabvect)
	if lnan != lvect:  #there are values in vector so keep
		abcolkp.append(idr4abuncol[i])

#Number of abundance columns with values
labkp = len(abcolkp)

#Create matrix of nrows = abund, ncols = spec
abmat = np.zeros((labkp,lspe))
#aberrmat = np.zeros((labkp,lspe))

mnabmat = np.zeros((labkp))
sdabmat = np.zeros((labkp))

for i in range(0,labkp):
	iss = np.where(idr4abcol==abcolkp[i])
	ssab = idr4abscc[iss]

	abvect = idr4data.field(abcolkp[i])-idr4data.field('FEH')-ssab
	mwabvect = abvect[imat]
	ireal = np.where(np.isfinite(mwabvect)==True)
	rmwabvect = mwabvect[ireal]

	#Convert NaN to zero -> doing later to whole matrix
	#inan = np.where(np.isfinite(mwabvect)==False)
	#mwabvect[inan] = 0

	abmat[i,:] = mwabvect
	#aberrmat[i,:] = mwabvect

	mnabmat[i] = np.mean(rmwabvect)
	sdabmat[i] = np.std(rmwabvect)

#-------------------------------------
##Create Abudnaces FITS table
##defined columns -> abundances
#cna = idr4cname[imat]
#coln = fits.Column(name='CNAME', format='20A', array=cna)
#colist = []
#colist.append(coln)

#for i in range(0,len(mnabmat)):
#	a = abmat[i,:]
#	cola = fits.Column(name=abcolkp[i], format='E', array=a)
#	colist.append(cola)

#cols = fits.ColDefs(colist)
#tbhdu = fits.BinTableHDU.from_columns(cols)
#tbhdu.writeto('GES_NGC1851_Abundances.fits')
#-------------------------------------

#-------------------------------------
#Straight Covariance Matrix
tabmat = abmat.transpose()
covmat = np.dot(abmat,tabmat)

#-------------------------------------
#Remove mean of abundance vector from each element of that abundance
#Size of matrix -> abmat.shape
#Mean of a row(0) or column(1) of a matrix -> abmat.mean(axis=0), but doen'st account for 0,nan

#Remove mean abundance from each element
dabmat = abmat - mnabmat[:, np.newaxis]
zdabmat = dabmat

#Remove nans and put zeros
inan = np.where(np.isfinite(dabmat)==False)
zdabmat[inan] = 0

tdabmat = dabmat.transpose()
dcovmat = np.dot(dabmat,tdabmat)

#-------------------------------------
#Scale by diagonal elements in jk for each element in dcovmat
#Create matrix of nrows = abund, ncols = abund
ndcovmat = np.zeros((labkp,labkp))
asndcovmat = np.zeros((labkp+1,labkp+1))

#Test on specific i and j - vary values
#i=9
#j=10
#nele = dcovmat[i,j]
#di = dcovmat[i,i]
#dj = dcovmat[j,j]
#nval = np.sqrt(np.multiply(di,dj))
#outnele = np.divide(nele,nval)

for i in range(0,len(mnabmat)):
	for j in range(0,len(mnabmat)):
		nele = dcovmat[i,j]
		di = dcovmat[i,i]
		dj = dcovmat[j,j]
		nval = np.sqrt(np.multiply(di,dj))
		ndcovmat[i,j] = np.divide(nele,nval)


#-------------------------------------
##Create FITS table
##defined columns -> abundances
#coln = fits.Column(name='Element', format='20A', array=abcolkp)

#colist = []
#colist.append(coln)
#for i in range(0,len(mnabmat)):
#	a = ndcovmat[i,:]
#	cola = fits.Column(name=abcolkp[i], format='E', array=a)
#	colist.append(cola)

#cols = fits.ColDefs(colist)
#tbhdu = fits.BinTableHDU.from_columns(cols)
#tbhdu.writeto('GES_MW_WG11_NormAbundCov.fits')

#abcolkp.append(ndcovmat)
#data = Table([abcolkp, ndcovmat], names=['Element', abcolkp])
#np.ascii.write(data, 'data_new.txt')
#-------------------------------------
corlim = 0.75
elevect = np.linspace(1,len(mnabmat),len(mnabmat))

plt.close() 
for j in range(0,len(mnabmat)):
	plt.subplot(4,10,j+1) 

	icor = np.where(abs(ndcovmat[:,j]) > corlim)

	plt.ion()
	plt.plot([0,39],[0,0],'-c')
	plt.plot([0,39],[0+corlim,0+corlim],':c')
	plt.plot([0,39],[0-corlim,0-corlim],':c')

	for i in range(0,len(icor)):
		plt.plot([elevect[icor[i]],elevect[icor[i]]],[-1.2,1.2],'--r')

	plt.plot([elevect[j],elevect[j]],[-1.2,1.2],'--b')
	
	plt.plot(elevect,ndcovmat[:,j],'ok',markersize=3) #,'MarkerFaceColor','k')

	plt.title(abcolkp[j], fontsize=8)
	if np.remainder(j,10) == 0:
		plt.ylabel('Normalised Covariance', fontsize=6)

	if j > 29:
		plt.xlabel('GES Elements', fontsize=6)

		#plt.text(6500,0.5,'<[Fe/H]>=%5.2f'%mnfeh+'+/-%4.2f'%sdfeh)  # +'+/- %5.2f' sdfeh)
	plt.axis([0,39,-1.2,1.2])
	plt.xticks(elevect, abcolkp, rotation='vertical', fontsize=4)
	plt.subplots_adjust(left=0.05, bottom=0.05, right=0.98, top=0.97, wspace=0.5, hspace=0.3)

	plt.show()

plt.savefig('Figures/elementsbynormcovariance.pdf')
		#plt.close() 


import numpy as np
from astropy.io import fits
import fnmatch

#Load iDR4 FFFT Columns to get abundance fields
idr4colfile = open('/home/cclare/GES/iDR5/FFFT/GESiDR5_FFFT/IDL/GES_iDR5_WGYY_FFFT_ColumnsDefinition.txt')
idr4coldata = np.loadtxt(idr4colfile, dtype={'names': ('keywrd', 'keywrdcomment','datatype','datacomment','units','unitscomment'),'formats': ('S50','S50','S50','S50','S50','S50')},delimiter=';')
idr4kywrd = idr4coldata['keywrd']
idr4kywrdcmt = idr4coldata['keywrdcomment']
idr4dtyp = idr4coldata['datatype']
idr4dtypcmt = idr4coldata['datacomment']
idr4dunts = idr4coldata['units']
idr4duntscmt = idr4coldata['unitscomment']

idr4abuncol = idr4kywrd[np.where(np.char.find(idr4kywrdcmt,'Abundance') > -1)]

#Load iDR4 recommended
idr4rec = '/home/cclare/GES/Repository/Dropbox/iDR4/GES_iDR4_WG15_Recommended.fits'
hdulist = fits.open(idr4rec)

#header
idr4hd = hdulist[1].header

#Retrieve Data
idr4data = hdulist[1].data

#Length of header
#lenhdr = len(hdulist[1].header)
#keywrd = []
#for i in range(0,lenhdr):
#	keywrd.append(hdulist[1].header[i])

#load idr4 data
idr4cname = idr4data.field('CNAME')
idr4gestyp = idr4data.field('GES_TYPE')
idr4gesfld = idr4data.field('GES_FLD')
idr4recsetup = idr4data.field('REC_SETUP')
idr4recwg = idr4data.field('REC_WG')
idr4teff = idr4data.field('TEFF')
idr4eteff = idr4data.field('E_TEFF')
idr4logg = idr4data.field('LOGG')
idr4elogg = idr4data.field('E_LOGG')
idr4feh = idr4data.field('FEH')
idr4efeh = idr4data.field('E_FEH')
idr4vel= idr4data.field('VRAD')
idr4evel = idr4data.field('E_VRAD')
idr4al1 = idr4data.field('AL1') - 6.46 - idr4feh
idr4eal1 = idr4data.field('E_AL1')
idr4limal1 = idr4data.field('UPPER_COMBINED_AL1')
idr4na1 = idr4data.field('NA1') - 6.3 - idr4feh
idr4ena1 = idr4data.field('E_NA1')
idr4limna1 = idr4data.field('UPPER_COMBINED_NA1')
idr4o1 = idr4data.field('O1') - 8.69 - idr4feh
idr4eo1 = idr4data.field('E_O1')
idr4limo1 = idr4data.field('UPPER_COMBINED_O1')
idr4mg1 = idr4data.field('MG1') - 7.55 - idr4feh
idr4emg1 = idr4data.field('E_MG1')
idr4limmg1 = idr4data.field('UPPER_COMBINED_MG1')
idr4li1= idr4data.field('LI1') - 3.28 - idr4feh
idr4eli1 = idr4data.field('E_LI1')
idr4limli1 = idr4data.field('UPPER_COMBINED_LI1')
idr4sr1 = idr4data.field('SR1') - 2.91 - idr4feh
idr4esr1 = idr4data.field('E_SR1')
idr4limsr1 = idr4data.field('UPPER_COMBINED_SR1')
idr4zr2= idr4data.field('ZR2') - 2.6 - idr4feh
idr4ezr2 = idr4data.field('E_ZR2')
idr4limzr2 = idr4data.field('UPPER_ZR2')
idr4ba2= idr4data.field('BA2') - 2.18 - idr4feh
idr4eba2 = idr4data.field('E_BA2')
idr4limba2 = idr4data.field('UPPER_BA2')
idr4la2= idr4data.field('LA2') - 1.18 - idr4feh
idr4ela2 = idr4data.field('E_LA2')
idr4limla2 = idr4data.field('UPPER_LA2')
idr4eu2= idr4data.field('EU2') - 0.52 - idr4feh
idr4eeu2 = idr4data.field('E_EU2')
idr4limeu2 = idr4data.field('UPPER_EU2')


#Read in Solar Chemical Composition
solf = open('Inputs/Abund_SolChemComp_iDR4.txt')
soldata = np.loadtxt(solf, dtype={'names': ('solcol', 'solz','solabscc'),'formats': ('S4', 'I2','f6')})
idr4abcol = soldata['solcol']
idr4abz = soldata['solz']
idr4abscc = soldata['solabscc']


"""
This is the python implementation of GMBCG using the radially weighted ECGMM.
J. Hao @ Fermilab
5/10/2011
"""

from gmbcgFinder import *



#-----setup catalog directories---
InputCatDir = '/home/jghao/research/data/des_mock/v3.04/obsCat/'
OutputCatDir = 'output dir'

galF = gl.glob(InputCatDir+'*.fit')
NF = len(galF)

#----read in file and prepare the input variable -----
i=0
cat=pf.getdata(galF[i],1)
#central = cat.field('central')
ra=cat.field('ra')
dec=cat.field('dec')
photoz=cat.field('PHOTOZ_GAUSSIAN')
mag=cat.field('mag_z')
gmr=cat.field('mag_g') - cat.field('mag_r')
gmrErr=np.sqrt(cat.field('magerr_g')**2+cat.field('magerr_r')**2)
rmi=cat.field('mag_r') - cat.field('mag_i')
rmiErr=np.sqrt(cat.field('magerr_r')**2+cat.field('magerr_i')**2)
imz=cat.field('mag_i') - cat.field('mag_z')
imzErr=np.sqrt(cat.field('magerr_i')**2+cat.field('magerr_z')**2)
zmy=cat.field('mag_z') - cat.field('mag_y')
zmyErr=np.sqrt(cat.field('magerr_z')**2+cat.field('magerr_y')**2)
objID=cat.field('id')
Idx=np.arange(len(ra))

bcgCandidateIdx = selectBCGcandidate(photoz,gmr,rmi,imz,zmy)

grCandIdx = bcgCandidateIdx[photoz[bcgCandidateIdx] < 0.4]
res= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=gmr,colorErr=gmrErr,mag=mag,bcgCandidateIdx=grCandIdx)

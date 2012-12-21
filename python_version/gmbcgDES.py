"""
This is the python implementation of GMBCG using the radially weighted ECGMM.
J. Hao @ Fermilab
5/10/2011
"""

from gmbcgFinder import *
import cPickle as p

#-----setup catalog directories---
InputCatDir = '/home/jghao/research/data/des_mock/v4.00/obsCat/'
TruCatDir = '/home/jghao/research/data/des_mock/v4.00/truthCat/'
OutputCatDir = '/home/jghao/research/data/des_mock/v4.00/truthCat/clusterCat/'

TruGalF = gl.glob(TruCatDir+'*.fit')
TruGalF.sort()

galF = gl.glob(InputCatDir+'*.fit')
galF.sort()
NF = len(galF)

Truth = 0
#----read in file and prepare the input variable -----
i=0
cat=pf.getdata(galF[i],1)
ok = cat.field('PHOTOZ_GAUSSIAN') <= 0.65
cat = cat[ok]
central = pf.getdata(TruGalF[i],1).field('central')

ra=cat.field('ra')
dec=cat.field('dec')
photoz=cat.field('PHOTOZ_GAUSSIAN')

if Truth == 1:
    mag = cat.field('omag')[:,3]
    gmr = cat.field('omag')[:,0] - cat.field('omag')[:,1]
    gmrErr = np.sqrt(cat.field('omag')[:,0]**2+cat.field('omag')[:,1]**2)
    rmi = cat.field('omag')[:,1] - cat.field('omag')[:,2]
    rmiErr = np.sqrt(cat.field('omag')[:,1]**2+cat.field('omag')[:,2]**2)
    imz = cat.field('omag')[:,2] - cat.field('omag')[:,3]
    imzErr = np.sqrt(cat.field('omag')[:,2]**2+cat.field('omag')[:,3]**2)
    zmy = cat.field('omag')[:,3] - cat.field('omag')[:,4]
    zmyErr = np.sqrt(cat.field('omag')[:,3]**2+cat.field('omag')[:,4]**2)
    central = cat.field('central')
else:    
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

if Truth == 1:
    bcgCandidateIdx = Idx[central == 1]
else:
    bcgCandidateIdx = selectBCGcandidate(photoz,gmr,rmi,imz,zmy)


grCandIdx = bcgCandidateIdx[photoz[bcgCandidateIdx] < 0.4]
resGR= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=gmr,colorErr=gmrErr,mag=mag,bcgCandidateIdx=grCandIdx)


riCandIdx = bcgCandidateIdx[(photoz[bcgCandidateIdx] >= 0.4)*(photoz[bcgCandidateIdx] < 0.75)]
resRI= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=rmi,colorErr=rmiErr,mag=mag,bcgCandidateIdx=riCandIdx)


izCandIdx = bcgCandidateIdx[(photoz[bcgCandidateIdx] >= 0.75)*(photoz[bcgCandidateIdx] < 1.15)]
resIZ= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=imz,colorErr=imzErr,mag=mag,bcgCandidateIdx=izCandIdx)


zyCandIdx = bcgCandidateIdx[(photoz[bcgCandidateIdx] >= 1.15)*(photoz[bcgCandidateIdx] < 1.4)]
resZY= gmbcgFinder(objID=objID,ra=ra, dec=dec, photoz=photoz,color=zmy,colorErr=zmyErr,mag=mag,bcgCandidateIdx=zyCandIdx)


"""
BCGIdx, BCGobjID, BCGalpha0, BCGalpha1,BCGmu0, BCGmu1, BCGsigma0,BCGsigma1,BCGaic1,BCGaic2,BCGamp


"""

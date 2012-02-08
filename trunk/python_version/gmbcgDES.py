"""
This is the python implementation of GMBCG using the radially weighted ECGMM.
J. Hao @ Fermilab
5/10/2011
"""

from gmbcgFinder import *

class Gal:
    def __init__(self):
        self.ID = None
        self.childID = None
        self.alpha = None
        self.mu = None
        self.sigma = None
        self.rich = None
        self.bic1 = None
        self.bic2 = None
 
#-----0.4L* in i-band ---------
def limi(x):
    A=np.exp(3.1638)
    k=0.1428
    lmi=A*x**k
    return(lmi)

#-----0.4L* in z-band ---------
def limz(x):
    """
    corresponding i_absmag <= -20.5 
    """
    A=np.exp(3.17)
    k=0.15
    lmz=A*x**k
    return(lmz)


#-----setup catalog directories---
InputCatDir = '/home/jghao/research/des_mock/collrunJun2011'
OutputCatDir = 'output dir'

galF = gl.glob(InputCatDir+'/*.fit')
NF = len(galF)

i=0
cat=pf.getdata(galF[i],1)
central = cat.field('central')
ra=cat.field('ra')
dec=cat.field('dec')
photoz=cat.field('z')
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
bcgCandidateIdx = Idx[central == 1]

galObj=[objID,gmr,gmrErr,rmi,rmiErr,imz,imzErr,mag,photoz,ra,dec,Idx]

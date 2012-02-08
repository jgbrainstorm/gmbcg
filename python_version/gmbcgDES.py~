"""
This is the python implementation of GMBCG using the radially weighted ECGMM.
J. Hao @ Fermilab
5/10/2011
"""

import numpy as np
import pyfits as pf
import pylab as pl
import esutil as es
import ecgmmPy as gmm
import rwecgmmPy as rwgmm
import scipy.stats as sts
import glob as gl

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

#-----setup catalog directories---
InputCatDir = 'input dir'
OutputCatDir = 'output dir'

galF = gl.glob(InputCatDir+'/*.fit')
NF = len(galF)

catid=1
coadd=pf.getdata('/home/jghao/research/data/coadd10_29_09/gmbcg_input_small_'+str(catid)+'.fit')
ra=coadd.field('ra')
dec=coadd.field('dec')
photoz=coadd.field('photoz')
mag=coadd.field('model_counts')[:,3]
bcgCandidateIdx=np.arange(0,100000)
gmr=coadd.field('gmr')
gmrErr=coadd.field('gmr_err')
rmi=coadd.field('rmi')
rmiErr=coadd.field('rmi_err')
imz=coadd.field('imz')
imzErr=coadd.field('imz_err')
objID=coadd.field('objid')
Idx=np.arange(len(ra0))
childIdx=np.zeros(len(ra))

galObj=[objID,gmr,gmrErr,rmi,rmiErr,imz,imzErr,mag,photoz,ra,dec,Idx,childIdx]

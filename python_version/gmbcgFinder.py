"""
This define the GMBCG finder for individual catalog

"""

import numpy as np
import pyfits as pf
import pylab as pl
import esutil as es
import ecgmmPy as gmm
import rwecgmmPy as rwgmm
import scipy.stats as sts
import glob as gl
import time

class BCG:
    def __init__(self):
        self.Idx = None
        self.childIdx = None
        self.alpha = None
        self.mu = None
        self.sigma = None
        self.rich = None
        self.bic1 = None
        self.bic2 = None

def bgDsty(color=None,mag=None,area=None):
    N=len(color)
    value=np.array(zip(color,mag))
    kde=sts.gaussian_kde(value.T)
    
        
def gmbcgFinder(objID=None,ra=None, dec=None, photoz=None,gmr=None,gmrErr=None,rmi=None,rmiErr=None,imz=None, imzErr=None,mag=None,bcgCandidateIdx=None):
    cra = ra[bcgCandidateIdx]
    cdec = dec[bcgCandidateIdx]
    cmag = mag[bcgCandidateIdx]
    cphotoz = photoz[bcgCandidateIdx]
    depth = 12
    h=es.htm.HTM(depth)
    DA=es.cosmology.Da(0,cphotoz,h=0.7)
    srad=np.rad2deg(1./DA)
    m1,m2,d12 = h.match(cra,cdec,ra,dec,srad,maxmatch=5000)
    r12=np.deg2rad(d12)*DA[m1]
    indices=(mag[m2]<=limz(cphotoz[m1]))*(cmag[m1] < mag[m2])
    m1 = m1[indices]
    m2 = m2[indices]
    h,rev = es.stat.histogram(m1, binsize=1, rev=True)
    BCG=[]
    startTime=time.time()
    for i in range(5000):
        print str(i)
        if rev[i] != rev[i+1]:
            indx = rev[ rev[i]:rev[i+1]]
            if photoz[m1[indx[0]]] < 0.4:
                color=gmr[m2[indx]]
                colorErr=gmrErr[m2[indx]]
            if photoz[m1[indx[0]]] >= 0.4 and photoz[m1[indx[0]]] < 0.75:
                color=rmi[m2[indx]]
                colorErr=rmiErr[m2[indx]]
            if photoz[m1[indx[0]]] >= 0.75 and photoz[m1[indx[0]]] <= 1.0:
                color=imz[m2[indx]]
                colorErr=imzErr[m2[indx]]
            if photoz[m1[indx[0]]] > 1.0:
                continue
            Ntot = len(indx)
            alpha=np.array([0.5,0.5])
            mu=np.array([sts.scoreatpercentile(color,per=70),sts.scoreatpercentile(color,per=40)])
            sigma=np.array([0.04,0.3])
            bic2,alpha,mu,sigma=rwgmm.bic2EM(color,colorErr,r12[indx],alpha,mu,sigma)
            bic1 = rwgmm.bic1EM(color,colorErr,r12[indx])[0]
            if bic2 < bic1:
                bcgi=candidateBCG()
                srt=np.argsort(sigma)
                bcgi.Idx = m1[indx[0]]
                bcgi.alpha=alpha[srt]
                bcgi.mu=mu[srt]
                bcgi.sigma=sigma[srt]
                bcgi.rich = Ntot * alpha[0]
                bcgi.childIdx=m2[indx]
                bcgi.bic1=bic1
                bcgi.bic2=bic2
                BCG.append(bcgi)
        endTime=time.time()
        elapseTime=endTime-startTime
                

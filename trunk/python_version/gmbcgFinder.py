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

#---- color vs redshift, based on mock v3.04---
def zgmr(z):
    gmr = 2.681*z + 0.665
    return gmr

def zrmi(z):
    rmi = 2.063*z -0.266
    return rmi

def zimz(z):
    imz = 1.667*z-0.708
    return imz

def zzmy(z):
    zmy = 0.192 * z + 0.223
    return zmy

#----the inverse from color to redshift, based on mock v3.04---
def gmrz(gmr):
    z = 1./2.681*gmr - 0.665/2.681
    return z

def rmiz(rmi):
    z = 1./2.063*rmi - (-0.266/2.063)
    return z

def imzz(imz):
    z = 1./1.667*imz - (-0.708/1.667)
    return z

def zmyz(zmy):
    z = 1./0.192 * zmy - (0.223/0.192)
    return z
#--------------------------------------------------------------

def selectBCGcandidate(photoz,gmr,rmi,imz,zmy):
    idx = np.arange(len(photoz))
    idxGR = idx[photoz < 0.4]
    idxRI = idx[(photoz >= 0.4)*(photoz < 0.75)]
    idxIZ = idx[(photoz >= 0.75)*(photoz < 1.15)]
    idxZY = idx[(photoz >= 1.15)]
    idxGRin = abs(gmr[idxGR]-zgmr(photoz[idxGR])) <= 0.2
    idxRIin = abs(rmi[idxRI]-zrmi(photoz[idxRI])) <= 0.2
    idxIZin = abs(imz[idxIZ]-zimz(photoz[idxIZ])) <= 0.2
    idxZYin = abs(zmy[idxZY]-zzmy(photoz[idxZY])) <= 0.2
    idxIN = np.concatenate((idxGR[idxGRin],idxRI[idxRIin], idxIZ[idxIZin],idxZY[idxZYin]))
    return idxIN
    
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




"""
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
"""

def bgDsty(color=None,mag=None,area=None):
    N=len(color)
    value=np.array(zip(color,mag))
    kde=sts.gaussian_kde(value.T)
    

def gmbcgFinder(objID=None,ra=None, dec=None, photoz=None,color=None,colorErr=None,mag=None,bcgCandidateIdx=None):
    #--define some quantity to be returned ----
    BCGIdx = []
    BCGobjID = []
    BCGalpha0 = []
    BCGalpha1 = []
    BCGmu0 = []
    BCGmu1 = []
    BCGsigma0 = []
    BCGsigma1 = []
    BCGntot = []
    BCGamp = []
    BCGaic1 = []
    BCGaic2 = []
    #-----------------------------------------
    cra = ra[bcgCandidateIdx]
    cdec = dec[bcgCandidateIdx]
    cmag = mag[bcgCandidateIdx]
    cphotoz = photoz[bcgCandidateIdx]
    Ncandidates = len(cphotoz)
    ridgeZ = np.zeros(Ncandidates)
    depth = 10
    h=es.htm.HTM(depth)
    Cosmo = es.cosmology.Cosmo(h=0.7)
    DA=Cosmo.Da(0,cphotoz)
    srad=np.rad2deg(1./DA)
    m1,m2,d12 = h.match(cra,cdec,ra,dec,srad,maxmatch=5000)
    r12=np.deg2rad(d12)*DA[m1]
    indices=(mag[m2]<=limz(cphotoz[m1]))*(cmag[m1] < mag[m2])
    m1 = m1[indices]
    m2 = m2[indices]
    h,rev = es.stat.histogram(m1, binsize=1, rev=True)
    startTime=time.time()
    for i in range(h.size):
        if rev[i] != rev[i+1]:
            print i
            indx = rev[ rev[i]:rev[i+1]]
            alpha=np.array([0.5,0.5])
            mu=np.array([sts.scoreatpercentile(color[m2[indx]],per=70),sts.scoreatpercentile(color[m2[indx]],per=40)])
            sigma=np.array([0.04,0.3])
            aic2,alpha,mu,sigma=rwgmm.aic2EM(color[m2[indx]],colorErr[m2[indx]],r12[indx],alpha,mu,sigma)
            aic1 = rwgmm.aic1EM(color[m2[indx]],colorErr[m2[indx]],r12[indx])[0]
            if aic2 < aic1:
                srt=np.argsort(sigma)
                BCGIdx.append(m1[indx[0]])
                BCGobjID.append(objID[m1[indx[0]]])
                BCGalpha0.append(alpha[srt[0]])
                BCGalpha1.append(alpha[srt[1]])
                BCGmu0.append(mu[srt[0]])
                BCGmu1.append(mu[srt[1]])
                BCGsigma0.append(sigma[srt[0]])
                BCGsigma1.append(sigma[srt[1]])
                BCGaic1.append(aic1)
                BCGaic2.append(aic2)
                BCGamp.append(len(indx)*alpha[srt[0]]) 
    endTime=time.time()
    elapseTime=endTime-startTime
    print '---elapsed time: ' + str(elapseTime)
    return np.array(BCGIdx), np.array(BCGobjID),np.array(BCGalpha0),np.array(BCGalpha1),np.array(BCGmu0),np.array(BCGmu1),np.array(BCGsigma0),np.array(BCGsigma1),np.array(BCGaic1),np.array(BCGaic2),np.array(BCGamp)


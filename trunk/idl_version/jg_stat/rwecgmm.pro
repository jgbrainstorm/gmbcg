;This is the translaton of rwecgmmPy.py
; J Hao 4/2/2012


function poisson,r
    const=2.
    res=const*r
    return, res
end

function nfw,r
    rhos=104.5
    rs=0.15
    x=r/rs
    nobj = N_elements(x)
    resf=fltarr(nobj)
    lt1= where(x lt 1)
    resf[lt1] = 1-2./sqrt(1-x[lt1]^2)*atanh(sqrt((1-x[lt1])/(x[lt1]+1)))
    gt1= where(x gt 1)
    resf[gt1] = 1-2./sqrt(x[gt1]^2-1)*np.atanh(sqrt((x[gt1]-1)/(x[gt1]+1)))
    eq1= where(x eq 1)
    resf[eq1] = 0.
    res=2*rhos*rs/(x^2-1)*resf*r
    return, res
end


function pcz1,c,delta,r,mu1,sigma1
    res=nfw(r)*exp(-(c-mu1)^2/2./(sigma1^2+delta^2))/sqrt(2.*3.1415926*(sigma1^2+delta^2))
    return, res
end

function pcz2,c,delta,r,mu2,sigma2
    res=poisson(r)*exp(-(c-mu2)^2/2./(sigma2^2+delta^2))/sqrt(2.*3.1415926*(sigma2^2+delta^2))
    return, res
end

function pz1,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    res=pcz1(c,delta,r,mu1,sigma1)*w1/(pcz1(c,delta,r,mu1,sigma1)*w1+pcz2(c,delta,r,mu2,sigma2)*w2)
    return, res
end

function pz2,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    res=pcz2(c,delta,r,mu2,sigma2)*w2/(pcz1(c,delta,r,mu1,sigma1)*w1+pcz2(c,delta,r,mu2,sigma2)*w2)
    return, res
end

function w1new,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    nobj = n_elements(c)
    res=total(pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2))/float(nobj)
    return, res
end

function w2new,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    nobj = n_elements(c)
    res=total(pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2))/float(nobj)
    return, res
end

function mu1new, c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    upper=total(c*pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1^2/(sigma1^2+delta^2))
    lower=total(pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1^2/(sigma1^2+delta^2))
    res=upper/lower
    return, res
end

function mu2new,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    upper=total(c*pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2^2/(sigma2^2+delta^2))
    lower=total(pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2^2/(sigma2^2+delta^2))
    res=upper/lower
    return, res
end

function sigma1new,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    upper=total((c-mu1)^2*pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1^4/(sigma1^2+delta^2)^2)
    lower=total(pz1(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma1^4/(sigma1^2+delta^2)^2)
    res=sqrt(upper/lower)
    return, res
end


function sigma2new,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2):
    upper=total((c-mu2)^2*pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2^2/(sigma2^2+delta^2)^2)
    lower=total(pz2(c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2)*sigma2^2/(sigma2^2+delta^2)^2)
    res=sqrt(upper/lower)
    return, res
end


function lkhood2,c,delta,r,w1,w2,mu1,mu2,sigma1,sigma2
    res=pcz1(c,delta,r,mu1,sigma1)*w1+pcz2(c,delta,r,mu2,sigma2)*w2
    res=total(alog(res))
    return, res
end


pro bic1EM,c,delta,r, BIC, mu, sigma   
    ;Input: color, color errors, radius to the center
    ;Output: BIC,mu,sigma
    nobj = n_elements(c)
    r[where(r eq 0.)]=0.00001
    weighted_stat,c,delta,mu,sigma,wsigma
    lkhood_b=total(alog(pcz2(c,delta,r,mu,sigma)))
    BIC=-2.*lkhood_b + 2.*alog(nobj)
    return 
end
            
pro aic1EM,c,delta,r,AIC,mu,sigma
    
    ;Input: color, color errors, radius to the center
    ;Output: BIC,mu,sigma
    nobj = n_elements(c)
    r[where(r eq 0.)]=0.00001
    weighted_stat,c,delta,mu,sigma,wsigma
    ;mu,sigma,aic,bic=ec.wstat(c,delta)
    lkhood_b=total(alog(pcz2(c,delta,r,mu,sigma)))
    AIC=-2.*lkhood_b + 2.*2.
    return
end
 
pro bic2EM,c,delta,r,alpha,mu,sigma,BIC
    
    ;Input: color, color errors, radius to the center, initial guess of alpha, mu and sigma
    ;Output: BIC,alpha,mu,sigma
    nobj = n_elements(c)
    r[where(r eq 0.)]=0.00001
    w1=alpha[0]
    w2=alpha[1]
    mu1=mu[0]
    mu2=mu[1]
    sigma1=sigma[0]
    sigma2=sigma[1]
    NIter=1000
    acc=0.00001
    for k=0, Niter -1 do begin
        w1_new=w1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        w2_new=w2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        sigma1_new=sigma1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        sigma2_new=sigma2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        mu1_new=mu1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        mu2_new=mu2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        lkhood_a=lkhood2(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        lkhood_b=lkhood2(c,delta,r,w1_new,w2_new,mu1_new,mu2_new,sigma1_new,sigma2_new)
        alpha[0]=w1_new
        alpha[1]=w2_new
        mu[0]=mu1_new
        mu[1]=mu2_new
        sigma[0]=sigma1_new
        sigma[1]=sigma2_new
        if (abs(lkhood_a - lkhood_b) le acc) then break
    endfor
    BIC=-2.*lkhood_b + 5.*alog(nobj)
    return
end

pro aic2EM,c,delta,r,alpha,mu,sigma,AIC
    
    ;Input: color, color errors, radius to the center, initial guess of alpha, mu and sigma
    ;Output: AIC,alpha,mu,sigma
    nobj = n_elements(c)
    r[where(r eq 0.)]=0.00001
    w1=alpha[0]
    w2=alpha[1]
    mu1=mu[0]
    mu2=mu[1]
    sigma1=sigma[0]
    sigma2=sigma[1]
    NIter=1000
    acc=0.00001
    for k=0, Niter -1 do begin
        w1_new=w1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        w2_new=w2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        sigma1_new=sigma1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        sigma2_new=sigma2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        mu1_new=mu1new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        mu2_new=mu2new(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        lkhood_a=lkhood2(c,delta,r,alpha[0],alpha[1],mu[0],mu[1],sigma[0],sigma[1])
        lkhood_b=lkhood2(c,delta,r,w1_new,w2_new,mu1_new,mu2_new,sigma1_new,sigma2_new)
        alpha[0]=w1_new
        alpha[1]=w2_new
        mu[0]=mu1_new
        mu[1]=mu2_new
        sigma[0]=sigma1_new
        sigma[1]=sigma2_new
        if (abs(lkhood_a - lkhood_b) le acc) then break
    endfor
    AIC=-2.*lkhood_b + 5.*2.
    return
end

pro rwecgmm

    print, 'loaded'

end

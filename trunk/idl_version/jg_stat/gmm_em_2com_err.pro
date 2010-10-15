;univariate gaussian mixture using EM algorithm. This is two
;componenet mixture.

;newly modified on 10/22/2007
;new update after remove the bug of weighted stat 9/14/2009
;err: measurement errors
;sigma: intrinsic scatter of the color. 

FUNCTION gauss,x,x_err,mu,sigma
         gau = 1.0/sqrt(2*!pi*(sigma^2.+x_err^2.))*exp(-(x-mu)^2/(2.*(sigma^2+x_err^2)))
         RETURN,gau
END


FUNCTION Loglhood,x,x_err,mu,sigma,alpha
         lh=alpha[0]*gauss(x,x_err,mu[0],sigma[0])+alpha[1]*gauss(x,x_err,mu[1],sigma[1])
         loglh=alog(lh)
         RETURN,loglh
END


FUNCTION p0,x_i,x_err,alpha,mu,sigma
        pp = alpha[0]*gauss(x_i,x_err,mu[0],sigma[0])/(alpha[0]*gauss(x_i,x_err,mu[0],sigma[0])+alpha[1]*gauss(x_i,x_err,mu[1],sigma[1]))
        RETURN,pp   
END

FUNCTION p1,x_i,x_err,alpha,mu,sigma
        pp = alpha[1]*gauss(x_i,x_err,mu[1],sigma[1])/(alpha[0]*gauss(x_i,x_err,mu[0],sigma[0])+alpha[1]*gauss(x_i,x_err,mu[1],sigma[1]))
        RETURN,pp   
END


pro gmm_em_2com_err,xx,xx_err,alpha,mu,sigma,force2=force2,aic=aic,robust=robust
   !except=0           ;to avoid the report of underflow or overflows. 
    acc=0.0000001
    M=2
    iter=500L
        
    if keyword_set(robust) then begin
        tst=where(xx_err gt 0)
        if tst[0] eq -1 then begin
            x=xx
            x_err=xx_err
        endif else begin
            outlier_rem,xx_err,select,rem
            x=xx[select]
            x_err=xx_err[select]
        endelse

    endif else begin
            x=xx
            x_err=xx_err
    endelse
 
    N = N_ELEMENTS(x)

    num=Float(N)
    x=double(x)
    x_err=double(x_err)
    alpha=double(alpha)
    mu=double(mu)
    sigma=double(sigma)
    alpha_new = fltarr(M)
    mu_new = fltarr(M)
    sigma_new = fltarr(M)

    testerr=where(x_err gt 0)
    if testerr[0] ne -1 then begin

        FOR k=0L,iter DO BEGIN

            alpha_new[0] = total(p0(x,x_err,alpha,mu,sigma))/num
            mu_new[0]=total(x*p0(x,x_err,alpha,mu,sigma)*sigma[0]^2./(sigma[0]^2+x_err^2.))/total(p0(x,x_err,alpha,mu,sigma)*sigma[0]^2./(sigma[0]^2.+x_err^2.))
            sigma_new[0] = sqrt(total((x-mu_new[0])^2.*p0(x,x_err,alpha,mu,sigma)*sigma[0]^4./(sigma[0]^2.+x_err^2.)^2.)/total(p0(x,x_err,alpha,mu,sigma)*sigma[0]^2./(sigma[0]^2.+x_err^2.)))

            alpha_new[1] = total(p1(x,x_err,alpha,mu,sigma))/num
            mu_new[1]=total(x*p1(x,x_err,alpha,mu,sigma)*sigma[1]^2/(sigma[1]^2+x_err^2.))/total(p1(x,x_err,alpha,mu,sigma)*sigma[1]^2/(sigma[1]^2+x_err^2.))
            sigma_new[1] = sqrt(total((x-mu_new[1])^2.*p1(x,x_err,alpha,mu,sigma)*sigma[1]^4./(sigma[1]^2.+x_err^2.)^2.)/total(p1(x,x_err,alpha,mu,sigma)*sigma[1]^2./(sigma[1]^2.+x_err^2.)))
            
            loglh=total(loglhood(x,x_err,mu,sigma,alpha))
            loglh_new=total(loglhood(x,x_err,mu_new,sigma_new,alpha_new))

            if finite(loglh_new ) eq 0  then begin

                    BIC=1000
                    continue

             endif
             alpha=alpha_new
             mu=mu_new
             sigma=sigma_new
           
             if keyword_set(aic) then begin
                 BIC=-2.*loglh+5.*2. ;this is AIC
             endif else begin
                 BIC=-2*loglh+5.*alog(num) ;the BIC
             endelse
             
             IF(abs(loglh-loglh_new) le acc) Then begin
                 break
             ENDIF
         ENDFOR
       endif else begin

         FOR k=0L,iter DO BEGIN

            alpha_new[0] = total(p0(x,x_err,alpha,mu,sigma))/num
            mu_new[0]=total(x*p0(x,x_err,alpha,mu,sigma))/total(p0(x,x_err,alpha,mu,sigma))
            sigma_new[0] = sqrt(total((x-mu_new[0])^2.*p0(x,x_err,alpha,mu,sigma))/total(p0(x,x_err,alpha,mu,sigma)))

            alpha_new[1] = total(p1(x,x_err,alpha,mu,sigma))/num
            mu_new[1]=total(x*p1(x,x_err,alpha,mu,sigma))/total(p1(x,x_err,alpha,mu,sigma))
            sigma_new[1] = sqrt(total((x-mu_new[1])^2.*p1(x,x_err,alpha,mu,sigma))/total(p1(x,x_err,alpha,mu,sigma)))
            
            loglh=total(loglhood(x,x_err,mu,sigma,alpha))
            loglh_new=total(loglhood(x,x_err,mu_new,sigma_new,alpha_new))

            if finite(loglh_new ) eq 0 then begin

                    BIC=1000
                    continue

             endif
             alpha=alpha_new
             mu=mu_new
             sigma=sigma_new
           
             if keyword_set(aic) then begin
                 BIC=-2.*loglh+5.*2. ;this is AIC
             endif else begin
                 BIC=-2*loglh+5.*alog(num) ;the BIC
             endelse
             
             IF(abs(loglh-loglh_new) le acc) Then begin
                 break
             ENDIF
         ENDFOR

       endelse


       ;----------------single component---------------------
         weighted_stat,x,x_err,wxm,wsd,wsdm       
         
         if n_elements(x) eq 1 then loglh_1com=total(alog(gauss(x,x_err,x,1))) else loglh_1com=total(alog(gauss(x,x_err,wxm,wsd))) ;for non cov GMM(not ecgmm) 

         if keyword_set(aic) then begin
             BIC_1com=-2*loglh_1com+2.*2 ;this is aic 
         endif else begin
             BIC_1com=-2*loglh_1com+2.*alog(num)
         endelse
         
         if not keyword_set(force2) Then begin    
             IF(BIC_1com le BIC) Then Begin
                 if n_elements(x) eq 1 then begin
                     sigma=1 
                     mu=x
                     aplha=1.
                 endif else begin 
                     sigma=wsd
                     mu=wxm
                     alpha=1.
                 endelse
             Endif
         endif

         Return
 end



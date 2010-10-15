;This pro calculate the weighted statistics. 
;wxm: weighted mean
;wsd: weighted standard deviation
;wsdm: weighted standard deviation to weighted mean
; reference: http://en.wikipedia.org/wiki/Weighted_mean
;            http://www.itl.nist.gov/div898/software/dataplot/refman2/ch2/weightsd.pdf
;J G Hao @ U of Michigan 12/17/2007
;          modified on 3/24/2009 for the standard deviations to the
;          weighted means
;9/8/2009 modified back to the traditional sd. For the one used in
;maxbcg ridgeline paper, it is changed to the name weighted_stat_new.pro

pro weighted_stat,xx,xxerr,wxm,wsd,wsdm
    ok=where(xxerr gt 0)
    if ok[0] ne -1 then begin
        x=xx[ok]
        xerr=xxerr[ok]
        nn=float(n_elements(x))
        if nn eq 1 then begin
            wxm=x
            wsd=xerr
            wsdm=xerr

        endif else begin
            w=1./xerr^2.
            wxm=total(w*x)/total(w) ; weighted mean
                                ;  wsd=sqrt(total(w*(x-wxm)^2)/(total(w)-1))  ; weighed standard deviation by spss
            wsd=sqrt(total(w*(x-wxm)^2)/(total(w)*(1-1./nn))) ; weighed standard deviation by NIST
            wsdm=sqrt(1./total(w))    ; naive sd to the weighted mean
           ; wsdm=sqrt(1./total(w)/(nn-1)*total((x-wxm)^2.*w)) ; modified weighted standard deviation to the mean by wiki
        endelse
    endif else begin
        nn=float(n_elements(xx))
        if nn eq 1 then begin
            wxm=xx    
            wsd=xxerr
            wsdm=xxerr                                                                                                                                     
        endif else begin
            
            wxm=mean(xx)
            wsd=stdev(xx)
            wsdm=stdev(xx)/sqrt(nn-1)
        endelse         

     endelse
     
    return

end

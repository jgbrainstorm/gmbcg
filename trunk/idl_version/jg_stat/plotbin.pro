pro plotbin,x,y,yerr=yerr,binsize=binsize,overplot=overplot,xtitle=xtitle,ytitle=ytitle,xrange=xrange,yrange=yrange,weighted=weighted,color=color,linefit=linefit,xstyle=xstyle,ystyle=ystyle,struct=struct,thick=thick,textcolors=textcolors,charsize=charsize,median=median,xbin=xbin,ybin=ybin,symsize=symsize,scatter=scatter,xlog=xlog,ylog=ylog,title=title

    x=float(x)
    y=float(y)
    

    if keyword_set(yerr) then yerr=float(yerr)
    

    h=histogram(x,binsize=binsize,reverse_indices=ri,OMIN=om)
    nbin=n_elements(h)

    xbin=fltarr(nbin)
    ybin=fltarr(nbin)
    ysdbin=fltarr(nbin)   ; sd
    ysdmbin=fltarr(nbin)  ;sd to the mean
    nnbin=fltarr(nbin)    ; # in each bin
    for i=0,nbin-1 do begin
        
        if ri[i+1] ge ri[i]+2 then begin
           d=ri[ri[i]:ri[i+1]-1]
 
           xbin[i]=mean(x[d])

           if keyword_set(weighted) then begin 

            weighted_stat_new,y[d],yerr[d], wmean, wsig, werr
            if keyword_set(median) then ybin[i]=median(y[d]) else ybin[i]=wmean
             ysdmbin[i]=werr
             ysdbin[i]=wsig
             nnbin[i]=n_elements(d)

           endif else begin

           nnbin[i]=n_elements(d)
           if keyword_set(median) then ybin[i]=median(y[d]) else ybin[i]=mean(y[d])
           ysdmbin[i]=stdev(y[d])/sqrt(float(nnbin[i]))
           ;ysdbin[i]=stdev(y[d])
           ysdbin[i]=robust_std(y[d])
           endelse
        endif
    endfor

       ok=where(xbin ne 0 and ybin ne 0)

       xbin=xbin[ok]
       ybin=ybin[ok]
       ysdbin=ysdbin[ok]
       nnbin=nnbin[ok]
       if keyword_set(scatter) then begin
           ysdmbin=ysdbin[ok]
       endif else begin
           ysdmbin=ysdmbin[ok]
       endelse


       stt=create_struct('x',0.,'y',0.,'y_sd',0.,'y_sdm',0.,'nn',0L)
       struct=replicate(stt,n_elements(xbin))
       struct.x=xbin
       struct.y=ybin
       struct.y_sd=ysdbin
       struct.y_sdm=ysdmbin
       struct.nn=nnbin
 
 



    if keyword_set(overplot) then begin
       oplot,xbin,ybin,psym=8,color=color,symsize=symsize
       errplot, xbin,ybin-ysdmbin,ybin+ysdmbin,color=color
    endif else begin

       plot,xbin,ybin,psym=8,xtitle=xtitle,ytitle=ytitle,color=color,xrange=xrange,yrange=yrange,xstyle=xstyle,ystyle=ystyle,charsize=charsize,symsize=symsize,xlog=xlog,ylog=ylog,title=title

       errplot,xbin,ybin-ysdmbin,ybin+ysdmbin
    endelse

    if keyword_set(linefit) then begin

      ; if keyword_set(weighted) then res = LINFIT(xbin,ybin,measure_errors=ysmdbin,/double,sigma=res_sigma,chisq=gmr_chi2) else res = LINFIT(xbin,ybin,/double,sigma=res_sigma,chisq=gmr_chi2)
       
      ;  oplot,[min(xbin)-0.1,max(xbin)+0.1],[min(xbin)-0.1,max(xbin)+0.1]*res[1]+res[0],col=color,thick=thick
    
      ; if res_sigma[0] eq 0 and res_sigma[1] eq 0 then legend,['Slope: '+ntostr(res[1]),'Intercept: '+ntostr(res[0])],charsize=1,/bottom,/right,textcolors=textcolors else legend,['Slope: '+ntostr(res[1])+textoidl('\pm')+ ntostr(res_sigma[1]),'Intercept: '+ntostr(res[0])+textoidl('\pm')+ ntostr(res_sigma[0])],charsize=1.5,/bottom,box=0,/right,textcolors=textcolors
   
    endif

    return

end

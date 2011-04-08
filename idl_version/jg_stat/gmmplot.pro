;This pro plot the GMM fitting result
; J G Hao 1/10/2008 @ U of Michigan

function gaussfn,x,mu,sigma

         res=exp(-(x-mu)^2./2./sigma^2.)/sqrt(2*!pi)/sigma
         return,res
end

pro gmmplot,alpha,mu,sigma,xrange=xrange,npoint=npoint,overplot=overplot,linestyle=linestyle,xtitle=xtitle,charsize=charsize
      
         x=xrange[0]+lindgen(npoint)*(xrange[1]-xrange[0])/float(npoint)

            fn=fltarr(n_elements(x),n_elements(alpha))
            ff=fltarr(n_elements(x))

            for i=0,n_elements(alpha)-1 do begin
                 
                 fn[*,i]=alpha[i]*gaussfn(x,mu[i],sigma[i])                 
                 ff=ff+fn[*,i]
            endfor
                 peak=float(max(ff))
                 fn=fn/peak
                 ff=ff/peak

         if keyword_set(overplot) then begin

            oplot,x,ff,psym=3,thick=4,linestyle=linestyle,symsize=3

            for i=0,n_elements(alpha)-1 do begin
        
                if i eq 0 then clr=!red
                if i eq 1 then clr=!blue
                if i eq 2 then clr=!green
                oplot,x,fn[*,i],psym=3,col=clr,linestyle=linestyle,symsize=3
               
            endfor

                
 
         endif else begin

                plot,x,ff,psym=3,thick=4,linestyle=linestyle,symsize=3,color=color,xtitle=xtitle,charsize=charsize

            for i=0,n_elements(alpha)-1 do begin
        
                if i eq 0 then clr=!red
                if i eq 1 then clr=!blue
                if i eq 2 then clr=!green
                oplot,x,fn[*,i],psym=3,color=clr,linestyle=linestyle,symsize=3
               
            endfor
                
 
      endelse

         
      

end

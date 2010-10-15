;given an array, remove the outlier and return a clean array

;Jiangang Hao @ FNAL
;9/10/2009

pro outlier_rem,inArray,okidx,badidx,harsh=harsh
    y=inArray
    qt1=hquantile(y,1)
    qt3=hquantile(y,3)
    IRQ=qt3-qt1
    if keyword_set(harsh) then begin
      lowfense=qt1-1.5*IRQ   ; the standard one, but probably too harsh for ecgmm
      highfense=qt3+1.5*IRQ
    endif else begin
      lowfense=qt1-3.*IRQ         ;change from 1.5 to 3.
      highfense=qt3+3.*IRQ
    endelse
    okidx=where(y ge lowfense and y le highfense)
    badidx=where(y lt lowfense or y gt highfense)
    return
end

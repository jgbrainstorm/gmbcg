;This pro return the quantile of an array
;J G Hao @ U of Michigan, 10/3/2006


FUNCTION hquantile,x,ith
        
         ary=x       
         n = N_ELEMENTS(ary)
         if n gt 1 then begin
             srt = sort(ary)
             ary = ary[srt]
             pos = round((n*ith+ith)/4.0-1.0)
         endif else begin
             pos=0
         endelse
         RETURN, ary[pos]
END

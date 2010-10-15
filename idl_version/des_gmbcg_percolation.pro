;This blending process consist of two steps. 1. close blending
;2. far blending

;modified on 2/25/09 by matching the xray clusters.  The current
;version is good. 

function percolate, bg

;----step 1 0.25 mpc merge-----------------
    srad=0.2/angdist_lambda(bg.photoz)
    depth=10
    htm_find_neighbors,depth,bg.ra,bg.dec,bg.ra,bg.dec,srad,ind1,ind2,dist,/output_dist
    dmag=bg[ind1].omag[3]-bg[ind2].omag[3]
    in=where(dmag lt 0)
    in1=ind1[in]
    in2=ind2[in]

    his2=histogram(in1,OMIN=om,reverse_indices=rj,binsize=1,min=0)
    yy=rem_dup(in1)
    inr1=in1[yy]
 
    s=sort(bg[inr1].omag[3])
    inr1=inr1[s]
    num4=n_elements(inr1)

     
    For j=0L,num4-1 do begin
         print,j
  
          if(rj[inr1[j]+1] gt rj[inr1[j]]) then begin 
      
             dd=rj[rj[inr1[j]]:rj[inr1[j]+1]-1]

             If(bg[in1[dd[0]]].used ne 1) then begin
                 
             within= where(abs(bg[in1[dd]].photoz-bg[in2[dd]].photoz) le 0.3 and bg[in1[dd]].nfw_lh ge bg[in2[dd]].nfw_lh*0.2)
             if within[0] ne -1 then bg[in2[dd[within]]].used =1


             Endif

          endif

      Endfor

     bg=bg[where(bg.used ne 1)]
   
 return,bg
 end 

 ;------step 2 ----lr200 merge--

 function percolate_r200,bg
   
    srad=1./angdist_lambda(bg.photoz)

    depth=10
    htm_find_neighbors,depth,bg.ra,bg.dec,bg.ra,bg.dec,srad,ind1,ind2,dist,/output_dist
    dmag=bg[ind1].omag[3]-bg[ind2].omag[3]
    in=where(dmag lt 0)
    in1=ind1[in]
    in2=ind2[in]

    his2=histogram(in1,OMIN=om,reverse_indices=rj,binsize=1,min=0)
    yy=rem_dup(in1)
    inr1=in1[yy]
    s=sort(bg[inr1].omag[3])
    inr1=inr1[s]
    num4=n_elements(inr1)

     
    For j=0L,num4-1 do begin
         print,j
  
          if(rj[inr1[j]+1] gt rj[inr1[j]]) then begin 
      
             dd=rj[rj[inr1[j]]:rj[inr1[j]+1]-1]

             If(bg[in1[dd[0]]].used ne 1) then begin
                     
               within= where(abs(bg[in1[dd]].photoz-bg[in2[dd]].photoz) le 0.05 and bg[in1[dd]].nfw_lh ge bg[in2[dd]].nfw_lh*0.2)
              
                if within[0] ne -1 then bg[in2[dd[within]]].used =1


             Endif

          endif

      Endfor

      
      ok=where(bg.used ne 1)

      bg=bg[ok]
      return,bg         

end



pro des_gmbcg_percolation,catdir,version


    bcg=mrdfits(catdir+'des_mock_v'+ntostr(version,4)+'_BCG_gmbcg_v2.5.fit',1)
    bcg=bcg[where(bcg.ngals le 300)]
    bcg.used = 0

  ;------- step 0 ---remove junks----------

   
   ; bcg=bcg[ok]
 
   ; bcg1com=bcg[where(bcg.gm_nn eq 1 and ((bcg.gm_gmr_wdh gt 0 and bcg.gm_gmr_wdh le 0.12) or (bcg.gm_rmi_wdh gt 0 and bcg.gm_rmi_wdh le 0.12)) and bcg.gm_scaled_ngals ge 15 and bcg.bcglh ge 1 and bcg.nfw_lh ge 10)]

  
  ;  bcg2com=bcg[where(bcg.gm_nn eq 2)]   
 
  ;  b1=percolate(bcg1com)
  ;  b2=percolate(bcg2com)

  ;  b=[b1,b2]
    bcg=percolate(bcg)
    bcg.used=0
    bg=percolate_r200(bcg)

 


    
   mwrfits,bg,catdir+'des_mock_v'+ntostr(version,4)+'_BCG_blended_gmbcg_v2.5.fit',/create
    
  

end

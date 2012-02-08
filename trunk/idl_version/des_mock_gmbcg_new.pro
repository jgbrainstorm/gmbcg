s                    ;------------------------------------------------
;Define arctanh function which is not exist in IDL
;------------------------------------------------

FUNCTION arctanh,x
        
         result=1/2.0*alog((1+x)/(1-x))
         
         RETURN, result

END


;--------------------------------------------------------
;Define the NFW profile function \Sigma(x)
;-------------------------------------------------------


FUNCTION NFW_radial,r

         rho_s=9.91227          ;It is a constant of normalization.see the mathematica fiel NFW.
                                ;Note the normalization is an annular
                                ;area, i.e. int 2\Pi*\sigma(r)*r dr=1   
         r_s=0.15               ;150 kpc. Here, r_s is in the unit of Mpc.
         x=double(r/r_s)        ;x=1 is a singular point, removed 8/28/2007

         fx=fltarr(N_elements(r))
         sig=fltarr(N_elements(r))
         dd=where(x eq 1.)
         d0=where(x le 0.1)
         d1=where(x gt 20.0)
         d2=where(x gt 0.1 and x lt 1.0)
         d3=where(x gt 1.0 and x le 20.0)
         IF(dd[0] ne -1) THEN fx[dd]=0.
         IF(d0[0] ne -1) THEN fx[d0]=6.03237
         IF(d1[0] ne -1) THEN fx[d1]=0.
         IF(d2[0] ne -1) THEN fx[d2]=(2.*rho_s*r_s)/(x[d2]^2-1)*(1-2/sqrt(1-x[d2]^2)*arctanh(sqrt((1-x[d2])/(x[d2]+1))))
         IF(d3[0] ne -1) THEN fx[d3]=(2.*rho_s*r_s)/(x[d3]^2-1)*(1-2/sqrt(x[d3]^2-1)*atan(sqrt((x[d3]-1)/(x[d3]+1))))
           
              
        
         RETURN,fx   
     END




;---------------------------------------------------------------
;Main engine for cluster selection.
;---------------------------------------------------------------



pro des_mock_gmbcg,cat_dir,gal,radius,patch,version

    
    t0=systime(1)  
    gal.gmr_err=0
    gal.rmi_err=0
    gal.imz_err=0
    gal.zmy_err=0

    select_des_mock_red,gal,inok   
    bright=where(gal[inok].omag[2] le gal[inok].lim_i and gal[inok].omag[2] le 22)

    inok=inok[bright]

    num=long(n_elements(inok))
    srad=double(radius/angdist_lambda(gal[inok].photoz))/(1+gal[inok].photoz) 

    depth=10
    htm_match,gal[inok].ra,gal[inok].dec,gal.ra,gal.dec,srad,ind1,ind2,dist,maxmatch=3000,depth=depth ;dist in radian
   
    dmag=gal[inok[ind1]].omag[3]-gal[ind2].omag[3] ; use zmag to compare the magnitude
    in=where(dmag le 0 and abs(gal[inok[ind1]].photoz-gal[ind2].photoz) le 0.2 and gal[ind2].omag[2] le gal[inok[ind1]].lim_i)
      
    in1=ind1[in]
    in2=ind2[in]
    dist=dist[in]

    ind1=0
    ind2=0
    in=0

     
     
;----------------------------------------------------------
;calcualte the bcg likelihood and mag likelihood for all galaxies. The
;bcg lh is the product of g-r lh and r-i lh
;----------------------------------------------------------
    
     his1=histogram(in1,reverse_indices=ri,binsize=1,OMIN=om,min=0)
     nhis1=n_elements(his1)

     For i=0L,nhis1-1 do begin

         print,'------------',i,'-------------------'

  
         If(ri[i+1]-1 ge ri[i]+10) then begin 
                 
            gg=ri[ri[i]:ri[i+1]-1]
            ;-------------------------gr ridgeline--------------------    
            if (gal[inok[in1[gg[0]]]].gr_ridge eq 1) then begin 

                alpha=[0.5,0.5]
              
                mu=[hquantile(gal[in2[gg]].gmr,3.),hquantile(gal[in2[gg]].gmr,0.8)]
                sigma=[0.04,0.3]
                gmm_em_2com_err,gal[in2[gg]].gmr,gal[in2[gg]].gmr_err,alpha,mu,sigma,/robust
                gal[inok[in1[gg[0]]]].Ntot=n_elements(gg) 
                if (n_elements(alpha) eq 2 and n_elements(mu) eq 2) then begin
                                      
                      alpha=[0.5,0.5]
                      mu=[hquantile(gal[in2[gg]].gmr,3.),hquantile(gal[in2[gg]].gmr,0.8)]
                      sigma=[0.04,0.3]
                      gmm_em_2com_err,gal[in2[gg]].gmr,0.,alpha,mu,sigma,/robust,/force2  
                      ss=reverse(sort(alpha*gauss(gal[inok[in1[gg[0]]]].gmr,0.,mu,sigma)))
                      ;ss=sort(sigma)
                      ss=sort()
                      gal[inok[in1[gg[0]]]].GM_mix_gmr_clr=alpha[ss[0]]
                      gal[inok[in1[gg[0]]]].GM_mix_gmr_bgd=alpha[ss[1]]

                      gal[inok[in1[gg[0]]]].GM_gmr=mu[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_gmr_bgd=mu[ss[1]]
 
                      gal[inok[in1[gg[0]]]].GM_gmr_wdh=sigma[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_gmr_wdh_bgd=sigma[ss[1]]
                      gal[inok[in1[gg[0]]]].GM_NN=2 
                      gal[inok[in1[gg[0]]]].gm_ngals_weighted=n_elements(gg)*gal[inok[in1[gg[0]]]].gm_mix_gmr_clr

                  endif else begin
                      gal[inok[in1[gg[0]]]].GM_NN=1
                      sigma=robust_std(gal[in2[gg]].gmr)  
                    
                      if abs(gal[inok[in1[gg[0]]]].gmr - mu) gt 3.*sigma then begin

                          gal[inok[in1[gg[0]]]].GM_mix_gmr_clr=0
                          gal[inok[in1[gg[0]]]].GM_mix_gmr_bgd=alpha

                          gal[inok[in1[gg[0]]]].GM_gmr=0 
                          gal[inok[in1[gg[0]]]].GM_gmr_bgd=mu
 
                          gal[inok[in1[gg[0]]]].GM_gmr_wdh=0
                          gal[inok[in1[gg[0]]]].GM_gmr_wdh_bgd=sigma

                      endif else begin

                          gal[inok[in1[gg[0]]]].GM_mix_gmr_clr=alpha
                          gal[inok[in1[gg[0]]]].GM_mix_gmr_bgd=0

                          gal[inok[in1[gg[0]]]].GM_gmr=mu
                          gal[inok[in1[gg[0]]]].GM_gmr_bgd=0
 
                          gal[inok[in1[gg[0]]]].GM_gmr_wdh=sigma
                          gal[inok[in1[gg[0]]]].GM_gmr_wdh_bgd=0

                      endelse
 
                  endelse

                  
                  within=where(abs(gal[in2[gg]].gmr-gal[inok[in1[gg]]].GM_gmr) le 2.*gal[inok[in1[gg]]].GM_gmr_wdh and gal[in2[gg]].omag[3] gt gal[inok[in1[gg]]].omag[3])
      
                  if within[0] ne -1 then begin
                    xd=(dist[gg[within]]*angdist_lambda(gal[inok[in1[gg[0]]]].photoz))
                    gal[inok[in1[gg[0]]]].bcglh=gauss_err(gal[inok[in1[gg[0]]]].gmr,gal[inok[in1[gg[0]]]].gmr_err,gal[inok[in1[gg[0]]]].GM_gmr,gal[inok[in1[gg[0]]]].GM_gmr_wdh);
                    gal[inok[in1[gg[0]]]].ngals=n_elements(xd)
                    gal[inok[in1[gg[0]]]].NFW_lh=total(NFW_radial(xd))
                    gal[inok[in1[gg[0]]]].lh=gal[inok[in1[gg[0]]]].NFW_lh*gal[inok[in1[gg[0]]]].bcglh
                   
                 endif
           
             endif
             ;------gr ridgeline end-------------------------

             ;-------------------------ri ridgeline--------------------    
             if (gal[inok[in1[gg[0]]]].ri_ridge eq 1) then begin 

                alpha=[0.5,0.5]
                mu=[hquantile(gal[in2[gg]].rmi,3.),hquantile(gal[in2[gg]].rmi,0.8)]     
                sigma=[0.04,0.3]
                gmm_em_2com_err,gal[in2[gg]].rmi,gal[in2[gg]].rmi_err,alpha,mu,sigma,/robust  
                gal[inok[in1[gg[0]]]].Ntot=n_elements(gg) 
                if (n_elements(alpha) eq 2 and n_elements(mu) eq 2) then begin
                    
                    alpha=[0.5,0.5] 
                    mu=[hquantile(gal[in2[gg]].rmi,3.),hquantile(gal[in2[gg]].rmi,0.8)]
                    sigma=[0.04,0.3]
                    gmm_em_2com_err,gal[in2[gg]].rmi,0.,alpha,mu,sigma,/force2,/robust 
                    ss=reverse(sort(alpha*gauss(gal[inok[in1[gg[0]]]].rmi,0.,mu,sigma)))
                      ;ss=sort(sigma)
                      gal[inok[in1[gg[0]]]].GM_mix_rmi_clr=alpha[ss[0]]
                      gal[inok[in1[gg[0]]]].GM_mix_rmi_bgd=alpha[ss[1]]

                      gal[inok[in1[gg[0]]]].GM_rmi=mu[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_rmi_bgd=mu[ss[1]]
 
                      gal[inok[in1[gg[0]]]].GM_rmi_wdh=sigma[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_rmi_wdh_bgd=sigma[ss[1]]
                      gal[inok[in1[gg[0]]]].GM_NN=2 
                      gal[inok[in1[gg[0]]]].gm_ngals_weighted=n_elements(gg)*gal[inok[in1[gg[0]]]].gm_mix_rmi_clr
              
                  endif else begin
                      gal[inok[in1[gg[0]]]].GM_NN=1

                      if abs(gal[inok[in1[gg[0]]]].rmi - mu) gt 3.*sigma then begin

                          gal[inok[in1[gg[0]]]].GM_mix_rmi_clr=0
                          gal[inok[in1[gg[0]]]].GM_mix_rmi_bgd=alpha

                          gal[inok[in1[gg[0]]]].GM_rmi=0 
                          gal[inok[in1[gg[0]]]].GM_rmi_bgd=mu
 
                          gal[inok[in1[gg[0]]]].GM_rmi_wdh=0 
                          gal[inok[in1[gg[0]]]].GM_rmi_wdh_bgd=sigma 

                      endif else begin

                          gal[inok[in1[gg[0]]]].GM_mix_rmi_clr=alpha
                          gal[inok[in1[gg[0]]]].GM_mix_rmi_bgd=0

                          gal[inok[in1[gg[0]]]].GM_rmi=mu
                          gal[inok[in1[gg[0]]]].GM_rmi_bgd=0
 
                          gal[inok[in1[gg[0]]]].GM_rmi_wdh=sigma  
                          gal[inok[in1[gg[0]]]].GM_rmi_wdh_bgd=0

                      endelse
 
                  endelse

                   within=where(abs(gal[in2[gg]].rmi-gal[inok[in1[gg]]].GM_rmi) le 2.*gal[inok[in1[gg]]].GM_rmi_wdh and gal[in2[gg]].omag[3] gt gal[inok[in1[gg]]].omag[3])

       
                 if within[0] ne -1 then begin
                    
                    xd=(dist[gg[within]]*angdist_lambda(gal[inok[in1[gg[0]]]].photoz))
                    gal[inok[in1[gg[0]]]].bcglh=gauss_err(gal[inok[in1[gg[0]]]].rmi,gal[inok[in1[gg[0]]]].rmi_err,gal[inok[in1[gg[0]]]].GM_rmi,gal[inok[in1[gg[0]]]].GM_rmi_wdh) 
                    gal[inok[in1[gg[0]]]].ngals=n_elements(xd)
                    gal[inok[in1[gg[0]]]].NFW_lh=total(NFW_radial(xd))
                    gal[inok[in1[gg[0]]]].lh=gal[inok[in1[gg[0]]]].NFW_lh*gal[inok[in1[gg[0]]]].bcglh
                  
                 endif
           
             endif
             ;------ri ridgeline end-------------------------

            
             ;-------------------------iz ridgeline--------------------    
             if (gal[inok[in1[gg[0]]]].iz_ridge eq 1) then begin 


                alpha=[0.5,0.5]
                mu=[hquantile(gal[in2[gg]].imz,3.),hquantile(gal[in2[gg]].imz,0.8)]     
                sigma=[0.04,0.3]
                gmm_em_2com_err,gal[in2[gg]].imz,gal[in2[gg]].imz_err,alpha,mu,sigma,/robust  
                gal[inok[in1[gg[0]]]].Ntot=n_elements(gg) 
                if (n_elements(alpha) eq 2 and n_elements(mu) eq 2) then begin
              
                     alpha=[0.5,0.5]
                     mu=[hquantile(gal[in2[gg]].imz,3.),hquantile(gal[in2[gg]].imz,0.8)] 
                     sigma=[0.04,0.3]
                     gmm_em_2com_err,gal[in2[gg]].imz,0.,alpha,mu,sigma,/force2,/robust
                      ss=reverse(sort(alpha*gauss(gal[inok[in1[gg[0]]]].imz,0.,mu,sigma)))                                    
                      ;ss=sort(sigma)
                      gal[inok[in1[gg[0]]]].GM_mix_imz_clr=alpha[ss[0]]
                      gal[inok[in1[gg[0]]]].GM_mix_imz_bgd=alpha[ss[1]]

                      gal[inok[in1[gg[0]]]].GM_imz=mu[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_imz_bgd=mu[ss[1]]
 
                      gal[inok[in1[gg[0]]]].GM_imz_wdh=sigma[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_imz_wdh_bgd=sigma[ss[1]]
                      gal[inok[in1[gg[0]]]].GM_NN=2 
                      gal[inok[in1[gg[0]]]].gm_ngals_weighted=n_elements(gg)*gal[inok[in1[gg[0]]]].gm_mix_imz_clr
              
                  endif else begin
                      gal[inok[in1[gg[0]]]].GM_NN=1
                     
                      if abs(gal[inok[in1[gg[0]]]].imz - mu) gt 3.*sigma then begin

                          gal[inok[in1[gg[0]]]].GM_mix_imz_clr=0
                          gal[inok[in1[gg[0]]]].GM_mix_imz_bgd=alpha

                          gal[inok[in1[gg[0]]]].GM_imz=0 
                          gal[inok[in1[gg[0]]]].GM_imz_bgd=mu
 
                          gal[inok[in1[gg[0]]]].GM_imz_wdh=0
                          gal[inok[in1[gg[0]]]].GM_imz_wdh_bgd=sigma 

                      endif else begin

                          gal[inok[in1[gg[0]]]].GM_mix_imz_clr=alpha
                          gal[inok[in1[gg[0]]]].GM_mix_imz_bgd=0

                          gal[inok[in1[gg[0]]]].GM_imz=mu
                          gal[inok[in1[gg[0]]]].GM_imz_bgd=0
 
                          gal[inok[in1[gg[0]]]].GM_imz_wdh=sigma
                          gal[inok[in1[gg[0]]]].GM_imz_wdh_bgd=0

                      endelse
 
                  endelse

                  
                  within=where(abs(gal[in2[gg]].imz-gal[inok[in1[gg]]].GM_imz) le 2.*gal[inok[in1[gg]]].GM_imz_wdh and gal[in2[gg]].omag[3] gt gal[inok[in1[gg]]].omag[3])

               
              
                 
                 if within[0] ne -1 then begin
                    xd=(dist[gg[within]]*angdist_lambda(gal[inok[in1[gg[0]]]].photoz))
                    gal[inok[in1[gg[0]]]].bcglh=gauss_err(gal[inok[in1[gg[0]]]].imz,gal[inok[in1[gg[0]]]].imz_err,gal[inok[in1[gg[0]]]].GM_imz,gal[inok[in1[gg[0]]]].GM_imz_wdh) 
                    gal[inok[in1[gg[0]]]].ngals=n_elements(xd)
                    gal[inok[in1[gg[0]]]].NFW_lh=total(NFW_radial(xd))
                    gal[inok[in1[gg[0]]]].gm_ngals_weighted = n_elements(gg)*gal[inok[in1[gg[0]]]].gm_mix_imz_clr
                    gal[inok[in1[gg[0]]]].lh=gal[inok[in1[gg[0]]]].NFW_lh*gal[inok[in1[gg[0]]]].bcglh
                 
                 endif
           
             endif
             ;------iz ridgeline end-------------------------


             ;-------------------------zy ridgeline--------------------    
             if (gal[inok[in1[gg[0]]]].zy_ridge eq 1) then begin 

                alpha=[0.5,0.5]
                mu=[hquantile(gal[in2[gg]].zmy,3.),hquantile(gal[in2[gg]].zmy,0.8)]     
                sigma=[0.04,0.3]
                gmm_em_2com_err,gal[in2[gg]].zmy,gal[in2[gg]].zmy_err,alpha,mu,sigma,/robust 
                gal[inok[in1[gg[0]]]].Ntot=n_elements(gg) 
                if (n_elements(alpha) eq 2 and n_elements(mu) eq 2) then begin
              
                      alpha=[0.5,0.5] 
                      mu=[hquantile(gal[in2[gg]].zmy,3.),hquantile(gal[in2[gg]].zmy,0.8)]
                      sigma=[0.04,0.3]
                      gmm_em_2com_err,gal[in2[gg]].zmy,0.,alpha,mu,sigma,/robust,/force2
                      ss=reverse(sort(alpha*gauss(gal[inok[in1[gg[0]]]].zmy,0.,mu,sigma)))
                     ; ss=sort(sigma)
                      gal[inok[in1[gg[0]]]].GM_mix_zmy_clr=alpha[ss[0]]
                      gal[inok[in1[gg[0]]]].GM_mix_zmy_bgd=alpha[ss[1]]

                      gal[inok[in1[gg[0]]]].GM_zmy=mu[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_zmy_bgd=mu[ss[1]]
 
                      gal[inok[in1[gg[0]]]].GM_zmy_wdh=sigma[ss[0]] 
                      gal[inok[in1[gg[0]]]].GM_zmy_wdh_bgd=sigma[ss[1]]
                      gal[inok[in1[gg[0]]]].GM_NN=2 
                      gal[inok[in1[gg[0]]]].gm_ngals_weighted=n_elements(gg)*gal[inok[in1[gg[0]]]].gm_mix_zmy_clr
              
              
                  endif else begin
                      gal[inok[in1[gg[0]]]].GM_NN=1
                     
                      if abs(gal[inok[in1[gg[0]]]].zmy - mu) gt 3.*sigma then begin

                          gal[inok[in1[gg[0]]]].GM_mix_zmy_clr=0
                          gal[inok[in1[gg[0]]]].GM_mix_zmy_bgd=alpha

                          gal[inok[in1[gg[0]]]].GM_zmy=0 
                          gal[inok[in1[gg[0]]]].GM_zmy_bgd=mu
 
                          gal[inok[in1[gg[0]]]].GM_zmy_wdh=0 
                          gal[inok[in1[gg[0]]]].GM_zmy_wdh_bgd=sigma 

                      endif else begin

                          gal[inok[in1[gg[0]]]].GM_mix_zmy_clr=alpha
                          gal[inok[in1[gg[0]]]].GM_mix_zmy_bgd=0

                          gal[inok[in1[gg[0]]]].GM_zmy=mu
                          gal[inok[in1[gg[0]]]].GM_zmy_bgd=0
 
                          gal[inok[in1[gg[0]]]].GM_zmy_wdh=sigma  
                          gal[inok[in1[gg[0]]]].GM_zmy_wdh_bgd=0

                      endelse
 
                  endelse

                   within=where(abs(gal[in2[gg]].zmy-gal[inok[in1[gg]]].GM_zmy) le 2.*gal[inok[in1[gg]]].GM_zmy_wdh and gal[in2[gg]].omag[3] gt gal[inok[in1[gg]]].omag[3])
                  
              
                 if within[0] ne -1 then begin
                    xd=(dist[gg[within]]*angdist_lambda(gal[inok[in1[gg[0]]]].photoz))
                    gal[inok[in1[gg[0]]]].bcglh=gauss_err(gal[inok[in1[gg[0]]]].zmy,gal[inok[in1[gg[0]]]].zmy_err,gal[inok[in1[gg[0]]]].GM_zmy,gal[inok[in1[gg[0]]]].GM_zmy_wdh) 
                    gal[inok[in1[gg[0]]]].ngals=n_elements(xd)
                    gal[inok[in1[gg[0]]]].NFW_lh=total(NFW_radial(xd))
                    gal[inok[in1[gg[0]]]].gm_ngals_weighted=n_elements(gg)*gal[inok[in1[gg[0]]]].gm_mix_zmy_clr
                    gal[inok[in1[gg[0]]]].lh=gal[inok[in1[gg[0]]]].NFW_lh*gal[inok[in1[gg[0]]]].bcglh
                 
                 endif
           
             endif
             ;------iz ridgeline end-------------------------

         Endif
     Endfor
        
        
     
        x_gr=where(gal[inok[in1]].gr_ridge eq 1)
        x_ri=where(gal[inok[in1]].ri_ridge eq 1)
        x_iz=where(gal[inok[in1]].iz_ridge eq 1)
      ; x_zy=where(gal[inok[in1]].zy_ridge eq 1)

      
        in_gr=where(abs(gal[in2[x_gr]].gmr-gal[inok[in1[x_gr]]].GM_gmr) le 2.*gal[inok[in1[x_gr]]].GM_gmr_wdh and gal[in2[x_gr]].omag[3] gt gal[inok[in1[x_gr]]].omag[3] and abs(gal[inok[in1[x_gr]]].gmr-gal[inok[in1[x_gr]]].GM_gmr) le 2.*gal[inok[in1[x_gr]]].GM_gmr_wdh);;only 2 sig of nocov 10/1/09

        in_ri=where(abs(gal[in2[x_ri]].rmi-gal[inok[in1[x_ri]]].GM_rmi) le 2.*gal[inok[in1[x_ri]]].GM_rmi_wdh and gal[in2[x_ri]].omag[3] gt gal[inok[in1[x_ri]]].omag[3] and abs(gal[inok[in1[x_ri]]].rmi-gal[inok[in1[x_ri]]].GM_rmi) le 2.*gal[inok[in1[x_ri]]].GM_rmi_wdh) 
        
        in_iz=where(abs(gal[in2[x_iz]].imz-gal[inok[in1[x_iz]]].GM_imz) le 2.*gal[inok[in1[x_iz]]].GM_imz_wdh and gal[in2[x_iz]].omag[3] gt gal[inok[in1[x_iz]]].omag[3] and abs(gal[inok[in1[x_iz]]].imz-gal[inok[in1[x_iz]]].GM_imz) le 2.*gal[inok[in1[x_iz]]].GM_imz_wdh) 


      


        in_all=[x_gr[in_gr],x_ri[in_ri],x_iz[in_iz]]
 
        in1=in1[in_all]   
        in2=in2[in_all] 

        x_ri=0                  ;mem clear
        in_ri=0                 ;mem clear
        x_gr=0                  ;mem clear
        in_gr=0                 ;mem clear
        x_iz=0
        in_iz=0
       ; x_zy=0
       ; in_zy=0
       ; in_all=0

        num2=n_elements(in2)

        his2=histogram(in1,reverse_indices=rj,binsize=1,OMIN=om1,min=0)
        yy=rem_dup(in1)
        inr1=in1[yy]
        s=reverse(sort(gal[inr1].nfw_lh)) ;replaced with nfw_lh 
        inr1=inr1[s]
        num4=n_elements(inr1)
        dup=0.0

       

         bcgmbt=create_struct('objid',long64(0),'bcgid',long64(0),'gmr',0.,'rmi',0.,'imz',0.,'zmy',0.,'bcg_gmr',0.,'bcg_rmi',0.,'bcg_imz',0.,'bcg_zmy',0.,'omag',[0.,0.,0.,0.,0],'bcgomag',[0.,0.,0.,0.,0],'amag',[0.,0.,0.,0.,0],'bcgamag',[0.,0.,0.,0.,0],'ra',-999.D,'dec',-999.D,'bcgra',0.D,'bcgdec',0.D,'bcgstatus',0,'photoz',0.,'photoz_err',0.,'bcg_photoz',0.,'bcg_photoz_err',0.,'ngals',0,'lh',0.,'spz',0.,'bcg_spz',0.,'bcg_limi',0.,'bcg_gmr_err',0.,'bcg_rmi_err',0.,'bcg_imz_err',0.,'bcg_zmy_err',0.,'gmr_err',0.,'rmi_err',0.,'imz_err',0.,'zmy_err',0.,'bcglh',0.,'bcgmag_lh',0.,'NFW_lh',0.,'bcg_gr_lh',0.,'bcg_ri_lh',0.,'bcg_iz_lh',0.,'rcenter',0.,'concentration',0.,'gr_ridge',0.,'ri_ridge',0.,'iz_ridge',0.,'zy_ridge',0.,'bcg_gr_ridge',0.,'bcg_ri_ridge',0.,'bcg_iz_ridge',0.,'bcg_zy_ridge',0.,'ngals_500kpc',0,'GM_gmr',0.,'GM_rmi',0.,'GM_imz',0.,'GM_zmy',0.,'GM_gmr_wdh',0.,'GM_rmi_wdh',0.,'GM_imz_wdh',0.,'GM_zmy_wdh',0.,'GM_mix_gmr_clr',0.,'GM_mix_gmr_bgd',0.,'GM_mix_rmi_clr',0.,'GM_mix_rmi_bgd',0.,'GM_mix_imz_clr',0.,'GM_mix_imz_bgd',0.,'GM_mix_zmy_clr',0.,'GM_mix_zmy_bgd',0.,'GM_gmr_bgd',0.,'GM_gmr_wdh_bgd',0.,'GM_rmi_bgd',0.,'GM_rmi_wdh_bgd',0.,'GM_imz_bgd',0.,'GM_imz_wdh_bgd',0.,'GM_zmy_bgd',0.,'GM_zmy_wdh_bgd',0.,'GM_NN',0, 'Ntot',0)



         bcgmb=replicate(bcgmbt,num2)
         
     bcgmb.ra=-999.D
     bcgmb.dec=-999.D
    
     For j=0L,num4-1 do begin
         
         if( rj[inr1[j]+1]-1 ge rj[inr1[j]]) then begin
             print,j
             dd=rj[rj[inr1[j]]:rj[inr1[j]+1]-1]
         
             if(gal[inok[in1[dd[0]]]].used ne 1. and gal[inok[in1[dd[0]]]].gm_nn eq 2) then begin

                gal[in2[dd]].used=1    ;used=1 means being used by other BCG. isbcg=999 means it is a BCG
                gal[inok[in1[dd[0]]]].isbcg=999.
                   
                bcgmb[dd].bcgid=gal[inok[in1[dd[0]]]].objid
                bcgmb[dd].bcgra=gal[inok[in1[dd[0]]]].ra
                bcgmb[dd].bcgdec=gal[inok[in1[dd[0]]]].dec
                bcgmb[dd].bcg_gmr=gal[inok[in1[dd[0]]]].gmr
                bcgmb[dd].bcg_rmi=gal[inok[in1[dd[0]]]].rmi
                bcgmb[dd].bcg_imz=gal[inok[in1[dd[0]]]].imz
                bcgmb[dd].bcg_zmy=gal[inok[in1[dd[0]]]].zmy

                bcgmb[dd].bcg_gmr_err=gal[inok[in1[dd[0]]]].gmr_err
                bcgmb[dd].bcg_rmi_err=gal[inok[in1[dd[0]]]].rmi_err
                bcgmb[dd].bcg_imz_err=gal[inok[in1[dd[0]]]].imz_err
                bcgmb[dd].bcg_zmy_err=gal[inok[in1[dd[0]]]].zmy_err

                bcgmb[dd].bcg_photoz=gal[inok[in1[dd[0]]]].photoz
                bcgmb[dd].bcg_photoz_err=gal[inok[in1[dd[0]]]].photoz_err
                bcgmb[dd].bcgomag=gal[inok[in1[dd[0]]]].omag
                bcgmb[dd].bcgamag=gal[inok[in1[dd[0]]]].amag
                bcgmb[dd].GM_NN=gal[inok[in1[dd[0]]]].GM_NN

                bcgmb[dd].bcgmag_lh=gal[in2[dd]].bcgmag_lh
                bcgmb[dd].bcglh=gal[in2[dd]].bcglh 
         
                bcgmb[dd].lh=gal[in2[dd]].lh
             
                bcgmb[dd].NFW_lh=gal[in2[dd]].nfw_lh
                bcgmb[dd].ngals=gal[in2[dd]].ngals
                bcgmb[dd].rcenter=dist[dd]*angdist_lambda(gal[inok[in1[dd[0]]]].photoz)
              
                bcgmb[dd].GM_gmr=gal[in2[dd]].GM_gmr
                bcgmb[dd].GM_rmi=gal[in2[dd]].GM_rmi
                bcgmb[dd].GM_imz=gal[in2[dd]].GM_imz
                bcgmb[dd].GM_zmy=gal[in2[dd]].GM_zmy
                
                bcgmb[dd].GM_gmr_wdh=gal[in2[dd]].GM_gmr_wdh
                bcgmb[dd].GM_rmi_wdh=gal[in2[dd]].GM_rmi_wdh
                bcgmb[dd].GM_imz_wdh=gal[in2[dd]].GM_imz_wdh
                bcgmb[dd].GM_zmy_wdh=gal[in2[dd]].GM_zmy_wdh
                   
                bcgmb[dd].ra=gal[in2[dd]].ra
                bcgmb[dd].dec=gal[in2[dd]].dec
                bcgmb[dd].gmr=gal[in2[dd]].gmr
                bcgmb[dd].rmi=gal[in2[dd]].rmi
                bcgmb[dd].imz=gal[in2[dd]].imz
                bcgmb[dd].zmy=gal[in2[dd]].zmy
 ss=reverse(sort(alpha*gauss(gal[inok[in1[gg[0]]]].rmi,0.,mu,sigma)))
                bcgmb[dd].objid=gal[in2[dd]].objid
                bcgmb[dd].photoz=gal[in2[dd]].photoz
                bcgmb[dd].photoz_err=gal[in2[dd]].photoz_err
                bcgmb[dd].omag=gal[in2[dd]].omag
                bcgmb[dd].gmr_err=gal[in2[dd]].gmr_err
                bcgmb[dd].rmi_err=gal[in2[dd]].rmi_err
                bcgmb[dd].imz_err=gal[in2[dd]].imz_err
                bcgmb[dd].zmy_err=gal[in2[dd]].zmy_err
                bcgmb[dd].gr_ridge=gal[in2[dd]].gr_ridge
                bcgmb[dd].ri_ridge=gal[in2[dd]].ri_ridge
                bcgmb[dd].iz_ridge=gal[in2[dd]].iz_ridge
                bcgmb[dd].zy_ridge=gal[in2[dd]].zy_ridge   
       
                   
            endif else begin
                dup=dup+1
                print,'duplicate',dup
           
            endelse
        endif

    Endfor


           in1=0;clear mem
           in2=0;clear mem


           ok=where(bcgmb.ra ne -999.0 and bcgmb.dec ne -999.0)
           bcg1=bcgmb[ok]
           bcgmb=0;clear mem
           ;okd=rem_dup(bcg1.bcgid)
           galok=where(gal.isbcg eq 999)

          
        bcgt=create_struct('objid',long64(0),'omag',[0.,0.,0.,0.,0],'amag',[0.,0.,0.,0.,0],'ra',0.D,'dec',0.D,'photoz',0.,'photoz_err',0.,'Ngals',0,'lh',0.,'z',0.,'lim_i',0.,'gmr',0.,'rmi',0.,'imz',0.,'zmy',0.,'flag',0.,'gmr_err',0.,'rmi_err',0.,'imz_err',0.,'zmy_err',0.,'r200',0.,'ngals_r200',0,'ngals_500kpc',0,'used',0,'bcglh',0.0,'bcgmag_lh',0.,'NFW_lh',0.,'bcg_gr_lh',0.,'bcg_ri_lh',0.,'bcg_iz_lh',0.,'rcenter',0.,'concentration',0.,'gr_ridge',0.,'ri_ridge',0.,'iz_ridge',0.,'zy_ridge',0.,'GM_gmr',0.,'GM_rmi',0.,'GM_imz',0.,'GM_zmy',0.,'GM_gmr_wdh',0.,'GM_rmi_wdh',0.,'GM_imz_wdh',0.,'GM_zmy_wdh',0.,'GM_mix_gmr_clr',0.,'GM_mix_gmr_bgd',0.,'GM_mix_rmi_clr',0.,'GM_mix_rmi_bgd',0.,'GM_mix_imz_clr',0.,'GM_mix_imz_bgd',0.,'GM_mix_zmy_clr',0.,'GM_mix_zmy_bgd',0.,'GM_gmr_bgd',0.,'GM_gmr_wdh_bgd',0.,'GM_rmi_bgd',0.,'GM_rmi_wdh_bgd',0.,'GM_imz_bgd',0.,'GM_imz_wdh_bgd',0.,'GM_zmy_bgd',0.,'GM_zmy_wdh_bgd',0.,'central',0,'GM_NN',0,'GM_Ngals_Weighted',0.,'Ntot',0.)
       
      bcg=replicate(bcgt,n_elements(galok))
      h=0.7

      bcg.objid=gal[galok].objid
      bcg.omag=gal[galok].omag
      bcg.amag=gal[galok].amag
      bcg.ra=gal[galok].ra
      bcg.dec=gal[galok].dec
      bcg.photoz=gal[galok].photoz
      bcg.photoz_err=gal[galok].photoz_err
      bcg.z=gal[galok].z
      bcg.gm_ngals_weighted = gal[galok].gm_ngals_weighted
      bcg.ntot=gal[galok].ntot
      bcg.ngals=gal[galok].ngals
      bcg.lh=gal[galok].lh
     ; bcg.ngals_500kpc=gal[galok].ngals_500kpc  

      bcg.lim_i=gal[galok].lim_i
      bcg.gmr=gal[galok].gmr
      bcg.rmi=gal[galok].rmi
      bcg.imz=gal[galok].imz
      bcg.zmy=gal[galok].zmy

      bcg.photoz_err=gal[galok].photoz_err
      bcg.gmr_err=gal[galok].gmr_err
      bcg.rmi_err=gal[galok].rmi_err
      bcg.imz_err=gal[galok].imz_err
      bcg.zmy_err=gal[galok].zmy_err

      bcg.r200=0.156*gal[galok].ngals^0.6/h
      bcg.bcglh=gal[galok].bcglh
  
      bcg.NFW_lh=gal[galok].NFW_lh
 
      bcg.rcenter=0.
  
      bcg.gr_ridge=gal[galok].gr_ridge
      bcg.ri_ridge=gal[galok].ri_ridge
      bcg.iz_ridge=gal[galok].iz_ridge
      bcg.zy_ridge=gal[galok].zy_ridge

      bcg.GM_gmr=gal[galok].GM_gmr                 
      bcg.GM_rmi=gal[galok].GM_rmi  
      bcg.GM_imz=gal[galok].GM_imz 
      bcg.GM_zmy=gal[galok].GM_zmy 

      bcg.GM_gmr_wdh=gal[galok].GM_gmr_wdh                 
      bcg.GM_rmi_wdh=gal[galok].GM_rmi_wdh 
      bcg.GM_imz_wdh=gal[galok].GM_imz_wdh 
      bcg.GM_zmy_wdh=gal[galok].GM_zmy_wdh

      bcg.central=gal[galok].central
      bcg.gm_NN=gal[galok].GM_NN
      bcg.GM_mix_gmr_clr=gal[galok].GM_mix_gmr_clr
      bcg.GM_mix_gmr_bgd=gal[galok].GM_mix_gmr_bgd
      bcg.GM_mix_rmi_clr=gal[galok].GM_mix_rmi_clr
      bcg.GM_mix_rmi_bgd=gal[galok].GM_mix_rmi_bgd
      bcg.GM_mix_imz_clr=gal[galok].GM_mix_imz_clr
      bcg.GM_mix_imz_bgd=gal[galok].GM_mix_imz_bgd
      bcg.GM_mix_zmy_clr=gal[galok].GM_mix_zmy_clr
      bcg.GM_mix_zmy_bgd=gal[galok].GM_mix_zmy_bgd

      
      bcg.GM_gmr_bgd=gal[galok].GM_gmr_bgd
      bcg.GM_gmr_wdh_bgd=gal[galok].GM_gmr_wdh_bgd
      bcg.GM_rmi_bgd=gal[galok].GM_rmi_bgd
      bcg.GM_rmi_wdh_bgd=gal[galok].GM_rmi_wdh_bgd
      bcg.GM_imz_bgd=gal[galok].GM_imz_bgd
      bcg.GM_imz_wdh_bgd=gal[galok].GM_imz_wdh_bgd
      bcg.GM_zmy_bgd=gal[galok].GM_zmy_bgd
      bcg.GM_zmy_wdh_bgd=gal[galok].GM_zmy_wdh_bgd
     


      t1=systime(1)
      print,'---------run time:---------'
      ptime,t1-t0
      mwrfits,bcg,cat_dir+'des_mock_v'+ntostr(version,4)+'_BCG_gmbcg_v2.5_stripe_'+ntostr(patch)+'.fit',/create
      mwrfits,bcg1,cat_dir+'des_mock_v'+ntostr(version,4)+'_BCGMB_gmbcg_v2.5_stripe_'+ntostr(patch)+'.fit',/create
   

    return
   
end


